#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TFile.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TMath.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <cassert>
#include <fstream>
#include <iomanip>
#include <algorithm>

// Configuration constants
const bool TOY_MC_ERROR_ESTIMATE = true;
const int TOY_MC_ERROR_ESTIMATE_SAMPLES = 10000;
const double Amc = 2.0 * TMath::Pi() * TMath::Pi();
//const double livetime = 217488673.3;// 9 years
const double livetime = 217488673.3 + 23361925.7586;//10 years
const double TotTime = livetime * Amc;
const double alpha = 2.7; // exponential to show flux_pow
const double egamma = 2.7; // exponential of prior power law
const int NATTEMPTS = 5; // max n. of unfolding iterations
int STARTING_DATA_BIN = 1;
const double STARTING_DATA_VAL = 1e2;//Min E value with N obs events > 0

const std::string Test_Stat = "chi2"; // "ks" (Kolmogorov-Smirnov) or "chi2" (Reduced Chi2)
double MIN_TS = (Test_Stat == "ks" ? 1e-4 : 1.);

// Global variables
std::vector<double> TRUGUESS;
std::vector<double> NORM_TRUE;

// Median energy calculation
double median_energy(double e_min, double e_max) {
    return TMath::Power((TMath::Power(e_max, 1.-alpha) - TMath::Power(e_min, 1.-alpha)) / 
                    ((1.-alpha) * (e_max - e_min)), -1.0/alpha);
}

void prior_prob(int nbins, Double_t *Ebins, std::vector<double>& p0) {
    double a_norm = (egamma - 1.) * TMath::Power(Ebins[5], egamma - 1.);
    double prior_sum{0.}, sum_check{0.};
    for (int ibin=0; ibin<nbins; ibin++){
        prior_sum += a_norm * (TMath::Power(Ebins[ibin+1], 1. - egamma) - TMath::Power(Ebins[ibin], 1. - egamma)) / (1. - egamma);
    }
    std::cout<<"prior_sum "<<prior_sum<<'\n';
    std::cout<<"Ebins[5] "<<Ebins[5]<<'\n';
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
                       std::vector<double>& truth, std::vector<double>& trutherr) {

    std::vector<double> obserr2(obserr.size());
    for(size_t i = 0; i < obserr.size(); i++) {
        obserr2[i] = obserr[i] * obserr[i];
    }

    truth = TRUGUESS;
    std::vector<double> originaltruth = TRUGUESS;
    std::vector<double> truthprevious;

    for(int attempt = 0; attempt < NATTEMPTS; attempt++) {
        //std::vector<double> smoothtru_prob = truth;//NO SMOOTHING
        std::vector<double> smoothtru_prob = smooth_(truth, originaltruth);
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
        std::vector<double> trutherr2 = invkernel_t.dot(obserr2);

        trutherr.resize(truth.size());

        for(size_t i = 0; i < truth.size(); i++) {
            trutherr[i] = std::sqrt(trutherr2[i]);
        }

        //Test statistic stopping criterion
        if(!truthprevious.empty()) {
            double ts{0.0};
            if (Test_Stat == "ks"){
                ts = compute_ks(truth, truthprevious);
            } else {
                ts = compute_chi2(truth, truthprevious);
            }
            std::cout << attempt << ") "<<Test_Stat<<" = " << ts << std::endl;
            if(ts < MIN_TS){
                attempt=NATTEMPTS;
                std::cout<< "Converged at iteration "<<attempt<< ": "<<Test_Stat<<" (" << ts <<") < " << MIN_TS << std::endl; 
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

int unfolding_smooth_ne() {

    // Parse command line arguments
    std::string response_file = "OUTPUTS/Out_MC_FTFP_NEON_E2e7_v8_Feb26_20bins_norm34.root";
    std::string response_histo = "h2Ntrig_wgt";
    std::string data_file = "OUTPUTS/Out_DATA_120months_NEON_Mar26_20bins.root";
    std::string data_histo = "h1Nobs_nobkg";
    //std::string ngen_file = ""; //uncomment if response-mat to be normalized
    //std::string ngen_histo = ""; //uncomment if response-mat to be normalized


    // Load normalized response matrix
    TFile* f = TFile::Open(response_file.c_str(), "READ");
    if(!f || f->IsZombie()) {
        std::cerr << "Error opening response file" << std::endl;
        return 1;
    }

    TH2F* h = (TH2F*)f->Get(response_histo.c_str());
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

    TH1D* hdata = (TH1D*)fdata->Get(data_histo.c_str());
    if(!hdata) {
        std::cerr << "Error getting data histogram" << std::endl;
        return 1;
    }

    int nbin_obs = hdata->GetNbinsX();
    int Nobs{0};
    for(int i = 0; i < nbin_obs; i++) {
        Nobs += hdata->GetBinContent(i+1);
    }
    std::cout << "Tot Nobs = " << Nobs << std::endl;

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    Matrix m(nx, ny);
    for(int i = 0; i < nx; i++) {
        for(int j = 0; j < ny; j++) {
            m.data[i][j] = h->GetBinContent(i+1, j+1);
        }
    }

    // Get truth guess
    Double_t *energy_bins = (Double_t*)h->GetXaxis()->GetXbins()->GetArray();
    TRUGUESS.resize(nx);
    prior_prob(nx, energy_bins, TRUGUESS);//filling TRUGUESS with power law prop. to E^{-2.6}

    //uncomment if response-mat to be normalized
    /*TFile* fgen = TFile::Open(ngen_file.c_str(), "READ");
    if(!fgen || fgen->IsZombie()) {
        std::cerr << "Error opening ngen file" << std::endl;
        return 1;
    }

    TH1D* hgen = (TH1d*)fgen->Get(ngen_histo.c_str());
    if(!hgen) {
        std::cerr << "Error getting response histogram" << std::endl;
        return 1;
    }
    NORM_TRUE.resize(nx);*/

    std::vector<double> eff;
    eff.resize(nx);
    std::vector<double> acc;
    acc.resize(nx);
    
    TH1D* truthguess = h->ProjectionX("hprior_guess");
    TH1D* hetru = h->ProjectionX("htru");
    TH1D* heff = (TH1D*)hetru->Clone();
    heff->SetDirectory(0);
    heff->SetName("eff");
    heff->SetTitle("eff");
    TH1D* hacc = (TH1D*)hetru->Clone();
    hacc->SetDirectory(0);
    hacc->SetName("acc");
    hacc->SetTitle("acc");

    for(int i = 0; i < nx; i++) {
        TRUGUESS[i] *= Nobs;
        truthguess->SetBinContent(i+1, TRUGUESS[i]);
        //NORM_TRUE[i] = hgen->GetBinContent(i+1); //uncomment if response-mat to be normalized
        eff[i] = hetru->GetBinContent(i+1);
        acc[i] = Amc*hetru->GetBinContent(i+1);
        heff->SetBinContent(i+1, eff[i]);
        heff->SetBinContent(i+1, hetru->GetBinError(i+1));
        hacc->SetBinContent(i+1, Amc*hetru->GetBinError(i+1));
    }
    //normalize_kernel(m, eff); //uncomment if response-mat to be normalized
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
    if(TOY_MC_ERROR_ESTIMATE) {
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

        for(int attempt = 0; attempt < TOY_MC_ERROR_ESTIMATE_SAMPLES; attempt++) {
            std::cout << "\nRunning toy MC " << attempt << std::endl;

            std::vector<double> newdata(ndata);
            for(int i = 0; i < ndata; i++) {
                newdata[i] = rng.Poisson(lambdas[i]) * coef[i];
                if(std::isnan(newdata[i])) newdata[i] = 0.0;
            }

            std::vector<double> toyresult, toyresulterr;
            Matrix toy_invm(m.rows, m.cols);
            iterative_bayesian(m, eff, newdata, dataerrs, toy_invm, toyresult, toyresulterr);
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

    TH1D* hflux = (TH1D*)hetru->Clone();
    hflux->SetDirectory(0);
    hflux->SetName("flux");
    hflux->SetTitle("flux");

    TH1D* hfluxpow = (TH1D*)hetru->Clone();
    hfluxpow->SetDirectory(0);
    hfluxpow->SetName("flux_pow");
    hfluxpow->SetTitle("flux_pow");

    std::ofstream outfile("DATA_POINTS/flux_spectrum_neon_E2e7_20bins_nobkg_norm34_5iter_toy10000_10y.dat");
    outfile << "# E  lowE  upE  Flux  StatErr\n";
    //outfile << std::setprecision(8) << std::scientific;

    for(int i = 0; i < nx; i++) {
        hresult->SetBinContent(i+1, result[i]);
        hresult->SetBinError(i+1, resulterr[i]);

        double lowE = hflux->GetXaxis()->GetBinLowEdge(i+1);
        double upE = hflux->GetXaxis()->GetBinUpEdge(i+1);

        double y = result[i] / ((upE - lowE) * TotTime);
        double yerr = resulterr[i] / ((upE - lowE) * TotTime);

        hflux->SetBinContent(i+1, y);
        hflux->SetBinError(i+1, yerr);

        double em = median_energy(lowE, upE);
        hfluxpow->SetBinContent(i+1, y * TMath::Power(em, alpha));
        hfluxpow->SetBinError(i+1, yerr * TMath::Power(em, alpha));

        outfile << em << "  " << lowE << "  " << upE << "  " << y << "  " << yerr << "\n";
    }

    hdata->SetLineColor(kRed);

    // Create testobs histogram
    TH1D* htestobs = (TH1D*)hdata->Clone();
    htestobs->SetDirectory(0);
    htestobs->SetName("testobs");
    htestobs->SetTitle("testobs");

    for(int i = 0; i < ndata; i++) {
        htestobs->SetBinContent(i+1, testobs[i]);
    }
    htestobs->SetLineColor(kGreen+2);

    TH2D* hunfold_matrix = (TH2D*)h->Clone();
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
    TString fout_name = "OUTPUTS/unfold_results_neon_E2e7_20bins_nobkg_norm34_5iter_toy10000_10y.root";
    TFile* fout = TFile::Open(fout_name, "RECREATE");
    fout->cd();
    h->Write();
    hunfold_matrix->Write();
    truthguess->Write();//prior distribution-guess
    htestobs->Write();//prediction of expected counts from unfolded flux
    hdata->Write();//counts before unfolding
    hresult->Write();//unfolded counts
    hflux->Write();
    hfluxpow->Write();
    //heff->Write();
    //hacc->Write();

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
    delete fout;

    return 0;
}
