#include <TH1F.h>
#include <TH1D.h>
#include <TH2.h>
#include <TFile.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cctype>
#include <string>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <array>
#include <limits>

// Configuration constants
// MC geometrical factor (unchanged from the original analysis).
const double Amc = 2.0 * TMath::Pi() * TMath::Pi() * 1.38 * 1.38;
const double alpha = 2.6; // exponential to show flux_pow
const double egamma = 2.7; // exponential of prior power law
const int NATTEMPTS = 10; // max n. of unfolding iterations
int STARTING_DATA_BIN = 1;
const double STARTING_DATA_VAL = 30;//Min E value with N obs events > 0

const std::string Test_Stat = "chi2"; // "ks" (Kolmogorov-Smirnov) or "chi2" (Reduced Chi2)
double MIN_TS = (Test_Stat == "ks" ? 1e-4 : 1.);

const bool SMOOTHING = true;

// IMPORTANT: this loader reads the LAST column, effective livetime [s].
// month = 1..12: exactly one month; month = 0: full year, requiring all 12.
bool get_effective_livetime(const std::string &filename, int year, int month,
                            double &effectiveTime)
{
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "ERROR: cannot open livetime table: " << filename << std::endl;
        return false;
    }

    std::array<double, 13> times{};
    std::array<bool, 13> found{};
    std::string line;
    int lineNo = 0;
    while (std::getline(input, line)) {
        ++lineNo;
        const auto start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos || line[start] == '#') continue;

        std::istringstream iss(line);
        std::string period;
        double tFirst, tLast, liveRaw, excluded, liveEffective;
        if (!(iss >> period >> tFirst >> tLast >> liveRaw >> excluded >> liveEffective)) {
            std::cerr << "ERROR: invalid livetime row at line " << lineNo << std::endl;
            return false;
        }
        if (period.size() != 7 || period[2] != '/' ||
            !std::isdigit(static_cast<unsigned char>(period[0])) ||
            !std::isdigit(static_cast<unsigned char>(period[1])) ||
            !std::isdigit(static_cast<unsigned char>(period[3])) ||
            !std::isdigit(static_cast<unsigned char>(period[4])) ||
            !std::isdigit(static_cast<unsigned char>(period[5])) ||
            !std::isdigit(static_cast<unsigned char>(period[6]))) {
            std::cerr << "ERROR: invalid MM/YYYY at line " << lineNo << std::endl;
            return false;
        }
        const int rowMonth = std::stoi(period.substr(0, 2));
        const int rowYear = std::stoi(period.substr(3, 4));
        if (rowMonth < 1 || rowMonth > 12) {
            std::cerr << "ERROR: month out of range at line " << lineNo << std::endl;
            return false;
        }
        if (rowYear != year || (month != 0 && rowMonth != month)) continue;
        if (found[rowMonth]) {
            std::cerr << "ERROR: duplicate livetime row for " << period << std::endl;
            return false;
        }
        if (!std::isfinite(liveEffective) || liveEffective <= 0) {
            std::cerr << "ERROR: non-positive effective livetime for " << period << std::endl;
            return false;
        }
        found[rowMonth] = true;
        times[rowMonth] = liveEffective;
    }

    effectiveTime = 0.0;
    if (month != 0) {
        if (!found[month]) {
            std::cerr << "ERROR: no livetime entry for " << month << "/" << year << std::endl;
            return false;
        }
        effectiveTime = times[month];
    } else {
        for (int m = 1; m <= 12; ++m) {
            if (!found[m]) {
                std::cerr << "ERROR: missing livetime entry for " << m << "/" << year
                          << " (required for a complete annual histogram)" << std::endl;
                return false;
            }
            effectiveTime += times[m];
        }
    }
    return true;
}

// Global variables
std::vector<double> TRUGUESS;
std::vector<double> NORM_TRUE;

// Median energy calculation
double median_energy(double e_min, double e_max) {
    return TMath::Power((TMath::Power(e_max, 1.-alpha) - TMath::Power(e_min, 1.-alpha)) / 
                    ((1.-alpha) * (e_max - e_min)), -1.0/alpha);
}

void prior_prob(int nbins, Double_t *Ebins, std::vector<double>& p0) {
    double a_norm = (egamma - 1.) * TMath::Power(Ebins[10], egamma - 1.);
    double prior_sum = 0.;
    double sum_check = 0.;
    for (int ibin=0; ibin<nbins; ibin++){
        prior_sum += a_norm * (TMath::Power(Ebins[ibin+1], 1. - egamma) - TMath::Power(Ebins[ibin], 1. - egamma)) / (1. - egamma);
    }
    std::cout<<"prior_sum "<<prior_sum<<'\n';
    std::cout<<"Ebins[10] "<<Ebins[10]<<'\n';
    for (int ibin=0; ibin<nbins; ibin++){
        double prior_pow = a_norm * (TMath::Power(Ebins[ibin+1], 1. - egamma) - TMath::Power(Ebins[ibin], 1. - egamma)) / (1. - egamma);
        p0[ibin] = prior_pow / prior_sum;
        sum_check += p0[ibin];
    }
    if(sum_check<0.999 || sum_check>1.001){
        std::cerr << "Error: prior probability not normalized!" << std::endl;
    }
}

// Matrix operations helper class
class Matrix {
public:
    std::vector<std::vector<double>> data;
    int rows, cols;

    Matrix(int r, int c) : rows(r), cols(c) {
        data.resize(r, std::vector<double>(c, 0.0));
    }

    Matrix transpose() {
        Matrix result(cols, rows);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result.data[j][i] = data[i][j];
            }
        }
        return result;
    }

    std::vector<double> dot(const std::vector<double>& vec) {
        std::vector<double> result(rows, 0.0);
        for(int i = 0; i < rows; i++) {
            for(int j = 0; j < cols; j++) {
                result[i] += data[i][j] * vec[j];
            }
        }
        return result;
    }
};

// Normalize kernel function
void normalize_kernel(Matrix& kernel, std::vector<double>& efficiencies) {
    efficiencies.resize(kernel.rows);
    Matrix new_kernel(kernel.rows, kernel.cols);

    for(int i = 0; i < kernel.rows; i++) {
        double sum_kernel = 0.0;
        for(int j = 0; j < kernel.cols; j++) {
            sum_kernel += kernel.data[i][j];
        }
        efficiencies[i] = (NORM_TRUE[i] > 0.0 ? (sum_kernel /NORM_TRUE[i]) : 0.0);

        if(sum_kernel > 0.0) {
            for(int j = 0; j < kernel.cols; j++) {
                new_kernel.data[i][j] = (NORM_TRUE[i] > 0.0 ? (kernel.data[i][j] /NORM_TRUE[i]) : 0.0);
            }
        }
    }

    kernel = new_kernel;
}

// Smooth truth using ROOT histogram smoothing
std::vector<double> smooth_(const std::vector<double>& truth, const std::vector<double>& originaltruth) {
    int n = truth.size();
    assert(n == (int)originaltruth.size());
    
    TH1F h("smooth", "smooth", n, 0, n);
    
    // Fill histogram with truth/originaltruth ratios
    for(int i = 0; i < n; i++) {
        h.SetBinContent(i+1, originaltruth[i] != 0.0 ? truth[i] / originaltruth[i] : 0.0);
    }
    h.GetXaxis()->SetRange(STARTING_DATA_BIN, n);
    // Apply smoothing
    h.Smooth(1, "R");
    
    // Extract smoothed result multiplied by originaltruth
    std::vector<double> result(n);
    for(int i = 0; i < n; i++) {
        result[i] = h.GetBinContent(i+1) * originaltruth[i];
    }
    
    return result;
}

double compute_ks(const std::vector<double>& truth, const std::vector<double>& truthprevious){
    double ks{0.}, tru_sum{0.}, truprev_sum{0.};
    int ni = truth.size();
    std::vector<double> tru_cumsum, truprev_cumsum;
    tru_cumsum.resize(ni);
    truprev_cumsum.resize(ni);
    for(int i = 0; i < ni; i++) {
        tru_sum += truth[i];
        truprev_sum += truthprevious[i];   
    }
    for(int i = 0; i < ni; i++) {
        tru_cumsum[i] = truth[i];
        truprev_cumsum[i] = truthprevious[i];   
        if(i>0){
            tru_cumsum[i] += truth[i-1];
            truprev_cumsum[i] += truthprevious[i-1];
        }
        tru_cumsum[i] /= tru_sum;
        truprev_cumsum[i] /= truprev_sum;
        double diff_abs = TMath::Abs(tru_cumsum[i] - truprev_cumsum[i]);
        if(ks < diff_abs){ks = diff_abs;}
        
    }

    return ks;
}

double compute_chi2(const std::vector<double>& truth, const std::vector<double>& truthprevious){
    double chi2{0.}, tru_sum{0.}, truprev_sum{0.};
    int ni = truth.size();
    for(int i = 0; i < ni; i++) {
        tru_sum += truth[i];
        truprev_sum += truthprevious[i];   
    }
    for(int i = 0; i < ni; i++) {
        double sum_dist = ( (truth[i] + truthprevious[i]) > 0. ? (truth[i] + truthprevious[i]) : 1. );
        double diff_dist = (truth[i] * truprev_sum) - (truthprevious[i] * tru_sum);
        chi2 += diff_dist * diff_dist / sum_dist;
    }
    chi2 /= (tru_sum * truprev_sum * ni);



    /*for(size_t i = 0; i < truth.size(); i++) {
        double diff = truth[i] - truthprevious[i];
        if( truth[i] > 0. ){
            chi2 += diff * diff / truth[i];
        }
    }*/
    return chi2;
}


// Iterative Bayesian unfolding
void iterative_bayesian(const Matrix& kernel, const std::vector<double>& eff, const std::vector<double>& obs, 
                       const std::vector<double>& obserr, Matrix& invkernel,
                       std::vector<double>& truth, std::vector<double>& trutherr, bool verbose = true) {

    std::vector<double> obserr2(obserr.size());
    for(size_t i = 0; i < obserr.size(); i++) {
        obserr2[i] = obserr[i] * obserr[i];
    }

    truth = TRUGUESS;
    std::vector<double> originaltruth = TRUGUESS;
    std::vector<double> truthprevious;

    for(int attempt = 0; attempt < NATTEMPTS; attempt++) {
        std::vector<double> smoothtru_prob;
        if (SMOOTHING)
            smoothtru_prob = smooth_(truth, originaltruth);
        else
            smoothtru_prob = truth;//NO SMOOTHING
        double Ntrue{0.};
        for(std::vector<double>::iterator it = smoothtru_prob.begin(); it != smoothtru_prob.end(); ++it){
            Ntrue += *it;
        }
        for(size_t i = 0; i < smoothtru_prob.size(); i++){
            smoothtru_prob[i] /= Ntrue;
        }

        for(int i = 0; i < kernel.rows; i++) {
            double norm = 0.0;
            for(int j = 0; j < kernel.cols; j++) {
                invkernel.data[i][j] = smoothtru_prob[j] * kernel.data[i][j];
                norm += kernel.data[i][j] * smoothtru_prob[j];
            }

            if(norm > 0.0) {
                for(int j = 0; j < kernel.cols; j++) {
                    invkernel.data[i][j] /= norm;
                    if(eff[j] > 0.0 ) {
                        invkernel.data[i][j] /= eff[j];
                    } else{ invkernel.data[i][j] = 0.; }
                }
            } else if (norm < 0.0){
                for(int j = 0; j < kernel.cols; j++) {
                    invkernel.data[i][j] = 0.;
                }
            }
        }

        Matrix invkernel_t = invkernel.transpose();
        truth = invkernel_t.dot(obs);
        // Approximate statistical propagation at fixed unfolding weights.
        // For full uncertainty estimation use the original toy-MC option.
        trutherr.resize(truth.size());
        for (size_t j = 0; j < truth.size(); ++j) {
            double variance = 0.0;
            for (int i = 0; i < invkernel.rows; ++i) {
                const double w = invkernel.data[i][j];
                variance += w * w * obserr2[i];
            }
            trutherr[j] = std::sqrt(variance);
        }

        //Test statistic stopping criterion
        if(!truthprevious.empty()) {
            double ts{0.0};
            if (Test_Stat == "ks"){
                ts = compute_ks(truth, truthprevious);
            } else {
                ts = compute_chi2(truth, truthprevious);
            }
            if (verbose) std::cout << attempt << ") "<<Test_Stat<<" = " << ts << std::endl;
            if(ts < MIN_TS){
                attempt=NATTEMPTS;
                if (verbose) std::cout<< "Converged at iteration "<<attempt<< ": "<<Test_Stat<<" (" << ts <<") < " << MIN_TS << std::endl; 
            }
        }

        truthprevious = truth;
    }
    for(int i = 0; i < kernel.rows; i++) {
        for(int j = 0; j < kernel.cols; j++) {
            invkernel.data[i][j] *= eff[j];
        }
    }
}

// Compute standard deviation
std::vector<double> compute_std(const std::vector<std::vector<double>>& data) {
    if(data.empty()) return std::vector<double>();

    int n = data.size();
    int m = data[0].size();
    std::vector<double> means(m, 0.0);
    std::vector<double> stds(m, 0.0);

    for(int j = 0; j < m; j++) {
        for(int i = 0; i < n; i++) {
            means[j] += data[i][j];
        }
        means[j] /= n;
    }

    for(int j = 0; j < m; j++) {
        for(int i = 0; i < n; i++) {
            stds[j] += (data[i][j] - means[j]) * (data[i][j] - means[j]);
        }
        stds[j] = std::sqrt(stds[j] / n);
    }

    return stds;
}

// For ROOT: root -l -b -q 'unfolding_smooth_pHe_period.C(2016,1,"ROOT_FILES/PHe_skim_2016_01_...root","TXT_FILES/effective_livetime.txt",10000)'
// month=0: annual spectrum; month=1..12: monthly spectrum.
// dataFile is the EXACT output ROOT file from Load.
// nToys=0 is useful for a fast smoke test; default=10000 as original.
int unfolding_smooth_pHe_monthly(int year = 2023, int month = 0,
                                const char *dataFile = "",
                                const char *livetimeFile = "TXT_FILES/effective_livetime.txt",
                                int nToys = 10000) {
    if (year < 2000 || year > 2100 || month < 0 || month > 12 || nToys < 0) {
        std::cerr << "ERROR: invalid year/month/nToys" << std::endl;
        return 1;
    }
    if (!dataFile || !dataFile[0]) {
        std::cerr << "ERROR: provide the exact ROOT file produced by Load as third argument" << std::endl;
        return 1;
    }

    double livetime = 0.0;
    if (!get_effective_livetime(livetimeFile, year, month, livetime)) return 1;
    const double TotTime = livetime * Amc;
    const std::string yearStr = std::to_string(year);
    const std::string period = month == 0 ? yearStr :
        yearStr + "_" + (month < 10 ? "0" : "") + std::to_string(month);
    const std::string suffix =
        "_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC";
    gSystem->mkdir("ROOT_FILES", true);
    gSystem->mkdir("TXT_FILES", true);
    const std::string fout_name = "ROOT_FILES/unfold_results_pHe_" + period + suffix + ".root";
    const std::string fdat_name = "TXT_FILES/flux_spectrum_pHe_" + period + suffix + ".dat";

    // Keep response matrix, histograms and normalization from the original macro.
    const std::string response_file = "ROOT_FILES/PHe_MC_p_He_5PeV_5binperdecade_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_STKvert0e7_10TeV_EPOSLHC.root";
    const std::string response_histo = "h2Ntrig_wgt";
    const std::string data_file = dataFile;
    const std::string data_histo = "h1SelBGO_orb";
    const std::string ngen_file = "PHe_MC_FTFP_10TeV_EPOSLHC_h1Ngen_5binsPerDecade.root";
    const std::string ngen_histo = "h1Ngen";

    std::cout << "\n=== Unfolding " << period << " ==="
              << "\nData: " << data_file
              << "\nEffective livetime [s]: " << std::setprecision(12) << livetime
              << "\nExposure T_eff * Amc: " << TotTime
              << "\nToy MC samples: " << nToys
              << "\nOutput ROOT: " << fout_name
              << "\nOutput ASCII: " << fdat_name << std::endl;


    // Load normalized response matrix
    TFile* f = TFile::Open(response_file.c_str(), "READ");
    if(!f || f->IsZombie()) {
        std::cerr << "Error opening response file" << std::endl;
        return 1;
    }

    TH2* h = dynamic_cast<TH2*>(f->Get(response_histo.c_str()));
    if(!h) {
        std::cerr << "Error getting response histogram" << std::endl;
        return 1;
    }

    // Load data
    TFile* fdata = TFile::Open(data_file.c_str(), "READ");
    if(!fdata || fdata->IsZombie()) {
        std::cerr << "Error opening data file" << std::endl;
        return 1;
    }

    TH1* hdata = dynamic_cast<TH1*>(fdata->Get(data_histo.c_str()));
    if(!hdata) {
        std::cerr << "Error getting data histogram" << std::endl;
        return 1;
    }

    int nbin_obs = hdata->GetNbinsX();
    double Nobs = 0.0;
    for(int i = 0; i < nbin_obs; i++) {
        Nobs += hdata->GetBinContent(i+1);
    }
    std::cout << "Tot Nobs = " << Nobs << std::endl;
    if (!(Nobs > 0.0)) {
        std::cerr << "ERROR: no observed events in data histogram" << std::endl;
        return 1;
    }

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    if (hdata->GetNbinsX() != ny) {
        std::cerr << "ERROR: data bins (" << hdata->GetNbinsX()
                  << ") differ from response reconstructed-energy bins (" << ny << ")" << std::endl;
        return 1;
    }

    Matrix m(nx, ny);
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            m.data[i][j] = h->GetBinContent(i+1, j+1);
        }
    }

    // Get truth guess
    std::vector<Double_t> energy_bins(nx + 1);
    for (int i = 0; i < nx; ++i) energy_bins[i] = h->GetXaxis()->GetBinLowEdge(i + 1);
    energy_bins[nx] = h->GetXaxis()->GetBinUpEdge(nx);
    TRUGUESS.resize(nx);
    prior_prob(nx, energy_bins.data(), TRUGUESS);

    //uncomment if response-mat to be normalized
    TFile* fgen = TFile::Open(ngen_file.c_str(), "READ");
    if(!fgen || fgen->IsZombie()) {
        std::cerr << "Error opening ngen file" << std::endl;
        return 1;
    }

    TH1* hgen = dynamic_cast<TH1*>(fgen->Get(ngen_histo.c_str()));
    if(!hgen) {
        std::cerr << "Error getting response histogram" << std::endl;
        return 1;
    }
    if (hgen->GetNbinsX() != nx) {
        std::cerr << "ERROR: Ngen histogram has " << hgen->GetNbinsX()
                  << " bins, response truth axis has " << nx << std::endl;
        return 1;
    }
    NORM_TRUE.resize(nx);

    std::vector<double> eff(nx, 0.0);

    TH1D* truthguess = h->ProjectionX("hprior_guess");
    TH1D* hetru = h->ProjectionX("htru");
    TH1D* heff = dynamic_cast<TH1D*>(hetru->Clone("eff"));
    heff->SetDirectory(nullptr);
    heff->SetTitle("eff");
    TH1D* hacc = dynamic_cast<TH1D*>(hetru->Clone("acc"));
    hacc->SetDirectory(nullptr);
    hacc->SetTitle("acc");

    for (int i = 0; i < nx; ++i) {
        TRUGUESS[i] *= Nobs;
        truthguess->SetBinContent(i+1, TRUGUESS[i]);
        NORM_TRUE[i] = hgen->GetBinContent(i+1);
    }
    normalize_kernel(m, eff); // Normalize response by generated truth events.
    for (int i = 0; i < nx; ++i) {
        heff->SetBinContent(i+1, eff[i]);
        heff->SetBinError(i+1, 0.0); // MC efficiency uncertainty not computed here.
        hacc->SetBinContent(i+1, Amc * eff[i]);
        hacc->SetBinError(i+1, 0.0);
    }
    m = m.transpose();


    int ndata = hdata->GetNbinsX();
    STARTING_DATA_BIN = hdata->FindBin(STARTING_DATA_VAL+1.);
    std::vector<double> data(ndata);
    std::vector<double> dataerrs(ndata);

    for(int i = 0; i < ndata; i++) {
        data[i] = hdata->GetBinContent(i+1);
        dataerrs[i] = hdata->GetBinError(i+1);
    }

    // Perform unfolding
    std::vector<double> result, resulterr;
    Matrix invm(m.rows, m.cols);
    iterative_bayesian(m, eff, data, dataerrs, invm, result, resulterr);

    // Toy MC error estimate
    if (nToys > 0) {
        std::vector<double> lambdas(ndata);
        std::vector<double> coef(ndata);

        for(int i = 0; i < ndata; i++) {
            if(dataerrs[i] != 0.0) {
                lambdas[i] = (data[i] / dataerrs[i]) * (data[i] / dataerrs[i]);
                coef[i] = data[i] / lambdas[i];
            }

            if(data[i] == 0.0) {
                lambdas[i] = 0.0;
                coef[i] = 0.0;
            }
        }

        TRandom3 rng(1234);
        std::vector<std::vector<double>> toyresults;

        for(int attempt = 0; attempt < nToys; attempt++) {
            if (attempt % 500 == 0 || attempt + 1 == nToys)
                std::cout << "Running toy MC " << attempt + 1 << "/" << nToys << std::endl;

            std::vector<double> newdata(ndata);
            for(int i = 0; i < ndata; i++) {
                newdata[i] = rng.Poisson(lambdas[i]) * coef[i];
                if(std::isnan(newdata[i])) newdata[i] = 0.0;
            }

            std::vector<double> toyresult, toyresulterr;
            Matrix toy_invm(m.rows, m.cols);
            iterative_bayesian(m, eff, newdata, dataerrs, toy_invm, toyresult, toyresulterr, false);
            toyresults.push_back(toyresult);
        }

        resulterr = compute_std(toyresults);
    }

    // Calculate test observation
    std::vector<double> testobs = m.dot(result);

    // Create output histograms
    TH1D* hresult = (TH1D*)hetru->Clone();
    hresult->SetDirectory(0);
    hresult->SetName("result");
    hresult->SetTitle("result");

    // Create output histograms
    TH1D* hcts_unf = (TH1D*)hetru->Clone();
    hcts_unf->SetDirectory(0);
    hcts_unf->SetName("unfold_cts");
    hcts_unf->SetTitle("unfold_cts");

    TH1D* hflux = (TH1D*)hetru->Clone();
    hflux->SetDirectory(0);
    hflux->SetName("flux");
    hflux->SetTitle("flux");

    TH1D* hfluxpow = (TH1D*)hetru->Clone();
    hfluxpow->SetDirectory(0);
    hfluxpow->SetName("flux_pow");
    hfluxpow->SetTitle("flux_pow");

    std::ofstream outfile(fdat_name);
    if (!outfile) {
        std::cerr << "ERROR: cannot create output ASCII file: " << fdat_name << std::endl;
        return 1;
    }
    outfile << std::setprecision(10);
    //outfile << "# E  lowE  upE  Flux  StatErr\n";
    //TotalEnergyDAMPE[j]>>FluxDAMPE[j]>>Stat_errDAMPE[j]>>EnergyErrDAMPE[j]
    //outfile << std::setprecision(8) << std::scientific;

    for(int i = 0; i < nx; i++) {
        hresult->SetBinContent(i+1, result[i]);
        hresult->SetBinError(i+1, resulterr[i]);

        hcts_unf->SetBinContent(i+1, result[i]*hetru->GetBinContent(i+1));
        hcts_unf->SetBinError(i+1, resulterr[i]*hetru->GetBinContent(i+1));

        double lowE = hflux->GetXaxis()->GetBinLowEdge(i+1);
        double upE = hflux->GetXaxis()->GetBinUpEdge(i+1);

        double y = result[i] / ((upE - lowE) * TotTime);
        double yerr = resulterr[i] / ((upE - lowE) * TotTime);

        hflux->SetBinContent(i+1, y);
        hflux->SetBinError(i+1, yerr);

        double em = median_energy(lowE, upE);
        hfluxpow->SetBinContent(i+1, y * TMath::Power(em, alpha));
        hfluxpow->SetBinError(i+1, yerr * TMath::Power(em, alpha));

        //outfile << em << "  " << lowE << "  " << upE << "  " << y << "  " << yerr << "\n";
        //outfile << em << "\t" << y << "\t" << yerr << "\t0.0\t " << "\n";
        outfile << em << "\t" << y << "\t" << yerr << "\t0.0\t" << lowE << "\t" << upE << "\n";
    }

    hdata->SetLineColor(kRed);

    // Create testobs histogram
    TH1* htestobs = dynamic_cast<TH1*>(hdata->Clone("testobs"));
    htestobs->SetDirectory(0);
    htestobs->SetName("testobs");
    htestobs->SetTitle("testobs");

    for(int i = 0; i < ndata; i++) {
        htestobs->SetBinContent(i+1, testobs[i]);
    }
    htestobs->SetLineColor(kGreen+2);

    TH2* hunfold_matrix = dynamic_cast<TH2*>(h->Clone("unfold_matrix"));
    hunfold_matrix->SetDirectory(0);
    hunfold_matrix->SetName("unfold_matrix");
    hunfold_matrix->SetTitle("unfold_matrix");
    hunfold_matrix->GetXaxis()->SetTitle( h->GetYaxis()->GetTitle() );
    hunfold_matrix->GetYaxis()->SetTitle( h->GetXaxis()->GetTitle() );
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            hunfold_matrix->SetBinContent(j+1, i+1, invm.data[j][i]);
        }
    }


    // Save output
    TFile* fout = TFile::Open(fout_name.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "ERROR: cannot create output ROOT file: " << fout_name << std::endl;
        return 1;
    }
    fout->cd();
    h->Write();
    hunfold_matrix->Write();
    truthguess->Write();//prior distribution-guess
    htestobs->Write();//prediction of expected counts from unfolded flux
    hdata->Write();//counts before unfolding
    hresult->Write();//unfolded counts
    hflux->Write();
    hfluxpow->Write();
    hcts_unf->Write();
    heff->Write();
    hacc->Write();

    TH1D* hacc1 = (TH1D*)hetru->Clone("hacc1");
    //hacc1->SetDirectory(0);
    int maxibin=hacc1->GetNbinsX();
    for(int ibin=0; ibin<maxibin; ibin++){
        double cont=hacc1->GetBinContent(ibin+1);
        hacc1->SetBinContent(ibin+1,cont*Amc);
    }

    hetru->Write(); //efficiency
    hacc1->Write(); //acceptance

    fout->Close();

    std::cout << "Results saved to "<< fout_name << std::endl;
    // Cleanup
    delete f;
    delete fdata;
    delete fgen;
    delete fout;

    return 0;
}
