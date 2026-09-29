#include "TChain.h"
#include "TObject.h"
#include "TCut.h"
#include "TString.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
#include <string>
#include <iostream>
#include <cstring>
#include <fstream>
#include <dirent.h>
#include <vector>
#include <math.h>

//#include "../plot_data_functions.h"


using namespace std;


struct FitData : public TObject {
    TGraphAsymmErrors* g;                // grafico da usare
    std::vector<int> prior_group;        // per ogni punto: id del gruppo di prior (0,1,...)
    std::vector<double> sigma_tilde;     // opzionale se vuoi fornirle esternamente (puoi lasciarla vuota)
    std::vector<double> data_FluxR;      // per ogni punto: valore di FluxR (usato per SigmaTilde)
    std::vector<double> data_SysErrH;    // per ogni punto: Sys_errH corrispondente
    int nuisance_start = 5;              // indice del primo parametro di nuisance in par[]
};

struct SumTF1 { 
    SumTF1(const std::vector<TF1 *> & flist) : fFuncList(flist) {}
    double operator() (const double *x, const double *p) {
        double result = 0.;
        for (int i = 0; i < fFuncList.size(); ++i) {
            if(fFuncList[i]->EvalPar(x,p)>0){
                result += fFuncList[i]->EvalPar(x,p);
            }
        }
        return result; 
    } 
    std::vector<TF1*> fFuncList; 
};

struct graphAsymErr{  
    TString graph_title;
    std::vector<Double_t> gx;/*<E>*/
    std::vector<Double_t> gy;/*y*/
    std::vector<Double_t> gelx;/*<E>-Elo*/
    std::vector<Double_t> gely;/*yerrtot_lo*/
    std::vector<Double_t> gehx;/*Eup-<E>*/
    std::vector<Double_t> gehy;/*yerrtot_up*/
    Double_t gx_min;
    Double_t gx_max;
};


constexpr double prior_pow{2.6}; 
void compute_rel_abundances(int N_el, double epow, std::vector<TF1 *>element_fits, std::string func_name, TString fname);
Double_t log_pol(Double_t *x, Double_t *par);
void chek_sum_abund(int Nel, int ntot_Ebins, double rel_abund_int[Nel][30]);
int element_to_int(std::string el);
std::string int_to_element(int numb);
Color_t get_hcolor(int el);
double brokenmod(Double_t *x, Double_t *par);
Double_t brokenmod3(Double_t *x, Double_t *par);
Double_t brokenmod3_Wknee(Double_t *x, Double_t *par);
graphAsymErr read_file(std::string file_name, bool is_dampe);
graphAsymErr get_pow_data_withElw(graphAsymErr data, double pow, double prior_pow);
double compute_xlw(double e1, double e2, double pow);


int DAMPE_model(){
    double epow{2.6};
    // Create output file
    string output_file_name = "ROOT_FILES/DAMPE_composition_pHe_10GeV_10PeV_5Bin.root";
    TFile fout(output_file_name.c_str(), "RECREATE");
    fout.Close();

    string elements[2] = {"p", "He"};
    int N_el = sizeof(elements)/sizeof(elements[0]);
    std::vector<TF1 *> element_pol_fits, NeMgSi_pols;
    //string path_base = "/home/irene/Documents/DAMPE/scripts/allpart/allpart_with_composition_model/DAMPE_single_nuclei_fluxes/flux_pHeCOFe_PaperSoftening/";
    string path_base = "./TXT_FILES/";
    std::map<std::string, std::string> spectra_map;
            //spectra_map["DAMPE pHe"] = path_base+"dampe_pHe_2024_onlyStatErr.txt";
            spectra_map["DAMPE p"] = path_base+"DAMPE_p_2026_pHePaperDraft.txt";
            spectra_map["DAMPE He"] = path_base+"DAMPE_He_2026_pHePaperDraft.txt";
            

    

    double parameters_el[2][11] = { {4e-2, 2.8, 2.3, 3., 2.4, 1e3, 1e4, 1e5, 0.5, 1.0, 3.},
                                    {5e-3, 2.7, 2.5, 2.7, 2.7, 1e3, 2.5e4, 3.5e5, 0.5, 4.0, 7.}
    };

    double pars_LowLim[2][11] = { {1e-4,  2.6,  2.1,  2.6,  2.45, 9e2, 9e3, 1e5,  0.01,  0.01,  3.},
                                    {1e-4,  2.6,  2.1,  2.6,  2.55, 9e2, 2e4, 3e5,  0.01,  0.01, 3.}
    };

    double pars_UpLim[2][11] = { {3e4, 3., 2.6, 3.1, 2.6, 2e3, 4e4, 1.85e5,  10.0,  10.0,  10.},
                                    {3e4, 3., 2.6, 3.1, 2.75, 2e3, 4e4, 5e5,  10.0,  10.0,  10.}
    };

    TCanvas *c = new TCanvas("c", "c" ,1);
    c->SetLogx();
    c->SetLogy();
    c->SetGridx();
    c->SetGridy();
    TLegend *leg = new TLegend(0.5, 0.9, 0.9, 0.75); 
    leg->SetBorderSize(1);
    leg->SetTextFont(132);
    leg->SetLineColor(1);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    leg->SetFillColor(0);
    leg->SetFillStyle(4000);
    TMultiGraph *mg_single_spec = new TMultiGraph();
    for(int el = 0; el<N_el; el++){ 

        cout<<"\n\n------------ "<<elements[el]<<" data ------------\n";
        std::string map_first = "DAMPE "+elements[el];

        TGraphAsymmErrors *grae;

        graphAsymErr gr_dampe = read_file((spectra_map[map_first]).c_str(), true);
        gr_dampe.graph_title = map_first.c_str();

        graphAsymErr gr_pow_dampe;
        gr_pow_dampe = get_pow_data_withElw(gr_dampe, epow, prior_pow);

        int n_pts = gr_dampe.gx.size();        
        TVectorD ga_x(n_pts, &gr_pow_dampe.gx[0]);//qui
        TVectorD ga_y(n_pts, &gr_pow_dampe.gy[0]);//qui
        TVectorD vec_ex(60);
        vec_ex.Zero();
        TVectorD ga_elx(n_pts, &vec_ex[0]);
        TVectorD ga_ehx(n_pts, &vec_ex[0]);
        TVectorD ga_ely(n_pts, &gr_pow_dampe.gely[0]);//qui
        TVectorD ga_ehy(n_pts, &gr_pow_dampe.gehy[0]);//qui
        grae = new TGraphAsymmErrors(ga_x, ga_y, ga_elx, ga_ehx, ga_ely, ga_ehy);
        grae->SetName(gr_pow_dampe.graph_title);
        grae->SetTitle(gr_pow_dampe.graph_title);
        grae->SetFillStyle(0);
        grae->SetLineColor(get_hcolor(el+1));
        grae->SetMarkerColor(get_hcolor(el+1));
        grae->SetMarkerStyle(8);
        grae->SetMarkerSize(1);
        grae->GetXaxis()->SetLimits(1e2, 1e7);
        //grae->GetYaxis()->SetRangeUser(2e2, 3e4);
        mg_single_spec->Add(grae,"p");
        leg->AddEntry(grae, gr_pow_dampe.graph_title, "ep");

        //TString polinomial_name = ("smoothed BPL_"+elements[el]).c_str();
        //TF1 *polin = new TF1(polinomial_name, AsyErr_InterpolTF1(grae), 1e2, 1e7, 0);
        TString func_name = ("sbpl_"+elements[el]).c_str();
        TF1 *polin;
        TF1 *polin_noknee = new TF1(func_name,brokenmod3,1e2,1e7,12);
        polin_noknee->SetParameters(parameters_el[el][0], parameters_el[el][1] , parameters_el[el][2] , parameters_el[el][3] , parameters_el[el][4] , parameters_el[el][5], parameters_el[el][6], parameters_el[el][7], parameters_el[el][8], parameters_el[el][9], parameters_el[el][10]);
        for(int ipar=1; ipar<11; ipar++){
            polin_noknee->SetParLimits(ipar, pars_LowLim[el][ipar], pars_UpLim[el][ipar]);
        }
        polin_noknee->SetParNames("#Phi_{0}","#gamma_{0}", "#gamma_{1}", "#gamma_{2}", "#gamma_{3}", "E_{B1}", "E_{B2}", "E_{B3}", "s1", "s2", "s3");
        grae->Fit(func_name,"R"); 

        polin  = new TF1(func_name,brokenmod3_Wknee,1e2,1e7,12);
        double kneePHe{3e6};
        if(elements[el] == "He"){kneePHe = 6e6;}
        for(int ipar=1; ipar<11; ipar++){
            polin->FixParameter(ipar, polin_noknee->GetParameter(ipar));
            //polin->SetParLimits(ipar, pars_LowLim[el][ipar], pars_UpLim[el][ipar]);
        }
        polin->FixParameter(11, kneePHe);
        polin->SetParNames("#Phi_{0}","#gamma_{0}", "#gamma_{1}", "#gamma_{2}", "#gamma_{3}", "E_{B1}", "E_{B2}", "E_{B3}", "s1", "s2", "s3");
        grae->Fit(func_name,"R"); 
        
        leg->AddEntry(polin, func_name, "l");
        
        polin->SetLineColor(get_hcolor(el+1));
        TCanvas *cc = new TCanvas(("cc" + int_to_element(el+1)).c_str(), ("cc" + int_to_element(el+1)).c_str(),1);
        cc->SetLogx();
        cc->SetLogy();
        cc->SetGridx();
        cc->SetGridy();
        grae->Draw("AP");
        polin->Draw("same");
        cc->BuildLegend();

        element_pol_fits.push_back(polin);
        
        
   }
   //cout<<"\n HAWC Composition Model reproduced\n";
   TF1 *fsum = new TF1("Sum", SumTF1(element_pol_fits), 1e2, 1e7, 0);
   fsum->SetLineColor(kBlack);
   element_pol_fits.push_back(fsum);
   c->cd();
   leg->AddEntry(fsum, "SBPL sum", "l");
   mg_single_spec->GetXaxis()->SetLimits(1e2, 1e7);
   mg_single_spec->GetXaxis()->SetTitle("Energy [GeV]");
   mg_single_spec->GetYaxis()->SetRangeUser(2e2, 3e4);
   mg_single_spec->GetYaxis()->SetTitle(Form("E^{%.01f} #Phi [GeV^{%.01f} (m^{2} sr s)^{-1}]", epow, epow-1.));
   mg_single_spec->Draw("AP");
   fsum->Draw("same");
   for(auto &ipolin:element_pol_fits){
        ipolin->Draw("same");
   }
   leg->Draw();

   compute_rel_abundances(N_el, epow, element_pol_fits, "SBPL", output_file_name.c_str());
   

   return(0);
}
///////////////////////////////////////////////////////////

void compute_rel_abundances(int N_el, double epow, std::vector<TF1 *>element_fits, std::string func_name, TString fname){
   TFile fout(fname, "update");
   const int n_decades_p1 = 4;
   const int n_decades_p2 = 2;
   const int bins_per_decade_p1 = 5;
   const int bins_per_decade_p2 = 5;
   const int ntot_Ebins = (bins_per_decade_p1 * n_decades_p1) + (bins_per_decade_p2 * n_decades_p2);
   const double log_bin_width_p1 = 1. / bins_per_decade_p1;
   const double log_bin_width_p2 = 1. / bins_per_decade_p2;
   const double binning_factor_p1 = TMath::Power(10., log_bin_width_p1);
   const double binning_factor_p2 = TMath::Power(10., log_bin_width_p2);

   Double_t *energy_bins = new Double_t[ntot_Ebins + 1];
   for (int iBin = 0; iBin < ntot_Ebins+1; iBin++){
      if(iBin<(bins_per_decade_p1 * n_decades_p1)+1/*<=100TeV*/){
            energy_bins[iBin] = 10. * TMath::Power(binning_factor_p1, iBin);
      } else {
            energy_bins[iBin] = energy_bins[iBin-1] * binning_factor_p2;
      }
         
   }

   //double rel_abund_div[N_el][ntot_Ebins];
   double rel_abund_int[N_el][ntot_Ebins];

   TCanvas *c_abund = new TCanvas(("c_abund_"+func_name).c_str(), ("c_abund "+func_name).c_str(), 900, 700);
   c_abund->cd()->SetLogx();
   //c_abund->cd()->SetLogy();
   c_abund->cd()->SetGridx();
   c_abund->cd()->SetGridy();
   auto leg_abund = new TLegend(0.7, 0.9, 0.9, 0.75); 
   leg_abund->SetNColumns(2);
 
   TCanvas *c_flux = new TCanvas(("c_flux_"+func_name).c_str(), ("c_flux "+func_name).c_str(), 900, 700);
   c_flux->cd()->SetLogx();
   c_flux->cd()->SetLogy();
   c_flux->cd()->SetGridx();
   c_flux->cd()->SetGridy();
   auto leg_flux = new TLegend(0.3, 0.9, 0.9, 0.82); 
   leg_flux->SetNColumns(2);

   double sum_BinIntegr[ntot_Ebins], ibin_integr;
   for(int bin = 0; bin < ntot_Ebins; bin++){
      sum_BinIntegr[bin] = 0.;
      for(int el = 0; el < N_el; el++){
         if(energy_bins[bin]>=9.9){
            ibin_integr = element_fits.at(el)->Integral(energy_bins[bin], energy_bins[bin+1]);
         } else {
            ibin_integr = 0.;
         }
         sum_BinIntegr[bin] += ibin_integr;
         rel_abund_int[el][bin] = ibin_integr;
      }
   }
   
   for(int el = 0; el < N_el; el++){
        TF1 *fit = new TF1();
        element_fits.at(el)->Copy(*fit);
        TF1 *fit_dash = new TF1();
        element_fits.at(el)->Copy(*fit_dash);
        fit->SetLineColor(get_hcolor(el+1));
        fit->SetLineWidth(3);
        fit->SetRange(1e2, 1e6);
        fit->SetTitle((int_to_element(el+1)+" SBPL fit").c_str());
        fit_dash->SetLineColor(get_hcolor(el+1));
        fit_dash->SetLineWidth(1);
        fit_dash->SetLineStyle(9);
        fit_dash->GetXaxis()->SetTitleOffset(1.2);
        fit_dash->GetXaxis()->SetTitle("Energy [GeV]");
        fit_dash->GetYaxis()->SetTitle(Form("E^{%0.1f} #Phi [GeV^{%0.1f} (m^{2} sr s)^{-1}]", epow, epow-1.));
        fit_dash->GetYaxis()->SetRangeUser(2e2, 3e4);
        TH1D *h_ab = new TH1D(("h_ab_"+func_name+"_"+int_to_element(el+1)).c_str(),("normalized relative abundance: ("+func_name+"_"+int_to_element(el+1)+") / (sum of abundances)").c_str(), ntot_Ebins, energy_bins);
        h_ab->SetTitle(";Energy [GeV]; Rel. abundance");
        h_ab->SetDirectory(0);
        h_ab->SetStats(0);
        h_ab->SetLineColor(get_hcolor(el+1));
        h_ab->SetLineWidth(2);
        h_ab->GetXaxis()->SetTitleOffset(1.2);
        h_ab->GetYaxis()->SetRangeUser(0., .8);
        //cout<<(int_to_element(el+1))<<":\n";
        for(int bin = 0; bin < ntot_Ebins; bin++){
            if(sum_BinIntegr[bin] != 0.){
                rel_abund_int[el][bin] /= sum_BinIntegr[bin];
            }
            //cout<<rel_abund_int[el][bin]<<'\t';
            h_ab->SetBinContent(bin+1, rel_abund_int[el][bin]);
            //cout<<rel_abund_int[el][bin]<<'\n';
        }
        TH1D *h_ab_dash = (TH1D*)h_ab->Clone("h_ab_dash");
        h_ab_dash->SetDirectory(0);
        h_ab_dash->SetLineStyle(9);
        h_ab_dash->SetLineWidth(1);
        if(el == 0){
            c_abund->cd();
            h_ab_dash->Draw();
            h_ab->Draw("same");
            leg_abund->AddEntry(h_ab, (int_to_element(el+1)).c_str(), "l");
            c_flux->cd();
            fit_dash->Draw();
            fit->Draw("same");
            leg_flux->AddEntry(fit, (int_to_element(el+1)+" fit\t").c_str(), "l");
        } else {
            c_abund->cd();
            h_ab_dash->Draw("same");
            h_ab->Draw("same");
            
            leg_abund->AddEntry(h_ab, (int_to_element(el+1)).c_str(), "l");
            c_flux->cd();
            fit_dash->Draw("same");
            fit->Draw("same");
            leg_flux->AddEntry(fit, (int_to_element(el+1)+" fit").c_str(), "l");
        }
        h_ab->Write();
        //fit->Write();
    }
    c_abund->cd();
    leg_abund->Draw();
    c_flux->cd();
    element_fits.back()->Draw("same");
    TF1 *fit_dash_tot = new TF1();
    element_fits.back()->Copy(*fit_dash_tot);
    fit_dash_tot->SetLineStyle(9);
    fit_dash_tot->SetLineWidth(1);
    fit_dash_tot->Draw("same");
    element_fits.back()->SetRange(1e2, 1e6);
    element_fits.back()->SetLineWidth(3);
    element_fits.back()->Draw("same");
    leg_flux->AddEntry(element_fits.back(), "Sum of fits", "l");
    leg_flux->SetTextSize(0.03);
    leg_flux->Draw();

   fout.Close();

   chek_sum_abund(N_el, ntot_Ebins, rel_abund_int);

}

Double_t log_pol(Double_t *x, Double_t *par){
   /*2 parameters: 
         - par[0] is the A normalization const
         - par[1] is gamma1
         - par[2] is Z
         - par[3] is characteristic rigidity (cut off)*/
   Double_t fit_res;
   fit_res = par[0] + par[1] * TMath::Log10(x[0]) + par[2] * TMath::Power(TMath::Log10(x[0]), 2) + par[3] * TMath::Power(TMath::Log10(x[0]), 3) + par[4] * TMath::Power(TMath::Log10(x[0]), 4) + par[5] * TMath::Power(TMath::Log10(x[0]), 5) + par[6] * TMath::Power(TMath::Log10(x[0]), 6) + par[7] * TMath::Power(TMath::Log10(x[0]), 7) + par[8] * TMath::Power(TMath::Log10(x[0]), 8) + par[9] * TMath::Power(TMath::Log10(x[0]), 9) + par[10] * TMath::Power(TMath::Log10(x[0]), 10);
   return fit_res;
}

void chek_sum_abund(int Nel, int ntot_Ebins, double rel_abund_int[Nel][30]){
   double sum_in_bin;
   bool is_check_passed = true;
   for(int bin = 6; bin < ntot_Ebins; bin++){
      sum_in_bin = 0.;
      for(int el = 0; el < Nel; el++){
         sum_in_bin += rel_abund_int[el][bin];
         //if(bin<11){cout<<int_to_element(el+1)<<": "<<rel_abund_int[el][bin]<<'\n';}
      }
      if(fabs(sum_in_bin - 1.) > 0.00001){
         //cout<<"In bin n."<<bin+1<<" the sum of abundances is "<<sum_in_bin<<" != 1\n";
         printf("ATTENTION!!! In bin n. %i the sum of abundances is %.8lf != 1\n", bin+1, sum_in_bin);
         is_check_passed = false;
      }
   }
   if (is_check_passed) {
        printf("Chek on the abundances sum (for each bin) PASSED!\n");
   }
}

int element_to_int(std::string el){
   std::map<std::string, int>element_int;
      element_int["p"] = 1;
      element_int["He"] = 2;
      element_int["C"] = 3;
      element_int["O"] = 4;
      element_int["NeMgSi"] = 5;
      element_int["Fe"] = 6;
   return element_int[el];
}
std::string int_to_element(int numb){
   std::map<int, std::string>int_element;
      int_element[1] = "p";
      int_element[2] = "He";
      int_element[3] = "C";
      int_element[4] = "O";
      int_element[5] = "NeMgSi";
      int_element[6] = "Fe";
   return int_element[numb];
}

Color_t get_hcolor(int el){
   std::map<int, Color_t>colors;
      colors[1] = kOrange+8;//kRed;
      colors[2] = kOrange+1;//kAzure-3;//kAzure+1;
      colors[3] = kAzure+1;//kAzure+1;
      colors[4] = kBlue+1;//kAzure+1;
      colors[5] = kMagenta+1;//kPink+4
      colors[6] = kSpring-8;//kOrange-1;
   return colors[el];
   //darkness++;
}


Double_t brokenmod(Double_t *x, Double_t *par)
{
    const Double_t E  = x[0];
    if (E <= 0) return 0.0;

    // Parameters
    const Double_t J0 = par[0];

    // Gamma values before/after each break
    const Double_t g0 = par[1];
    const Double_t g1 = par[2];
    const Double_t g2 = par[3];
    //const Double_t g3 = par[4];

    // Break energies
    const Double_t E1 = par[4];
    const Double_t E2 = par[5];
    //const Double_t E3 = par[7];

    // Smoothness
    const Double_t s1  = par[6];
    const Double_t s2  = par[7];

    // Reference energy
    const Double_t E0 = par[8];

    // Main pre-break power law
    Double_t flux = J0 * TMath::Power(E / E0, -g0);

    // Array-friendly structure
    //const Double_t Ei[3]      = {E1, E2, E3};
    //const Double_t gammaL[3]  = {g0, g1, g2};
    //const Double_t gammaH[3]  = {g1, g2, g3};
    const Double_t Ei[2]     = {E1, E2};
    const Double_t gammaL[2] = {g0, g1};
    const Double_t gammaH[2] = {g1, g2};
    const Double_t s[2] = {s1, s2};


    for (int i = 0; i < 2; i++)
    {
        Double_t dgamma = gammaH[i] - gammaL[i];

        // Numerator factor: E-dependent part
        Double_t num = 1.0 + TMath::Power(E / Ei[i], 1.0 / s[i]);

        // Denominator factor: E0-dependent part (normalization)
        Double_t den = 1.0 + TMath::Power(E0 / Ei[i], 1.0 / s[i]);

        // Add this break
        flux *= TMath::Power(num / den, dgamma * s[i]);
    }

    return flux;
}

Double_t brokenmod3(Double_t *x, Double_t *par)
{
    const Double_t E  = x[0];
    if (E <= 0) return 0.0;

    // Parameters
    const Double_t J0 = par[0];

    // Gamma values before/after each break
    const Double_t g0 = par[1];
    const Double_t g1 = par[2];
    const Double_t g2 = par[3];
    const Double_t g3 = par[4];

    // Break energies
    const Double_t E1 = par[5];
    const Double_t E2 = par[6];
    const Double_t E3 = par[7];

    // Smoothness
    const Double_t s1  = par[8];
    const Double_t s2  = par[9];
    const Double_t s3  = par[10];

    // Reference energy
    const Double_t E0 = 1e2;

    // Main pre-break power law
    Double_t flux = J0 * TMath::Power(E / E0, -g0);

    // Array-friendly structure
    const Double_t Ei[3]      = {E1, E2, E3};
    const Double_t gammaL[3]  = {g0, g1, g2};
    const Double_t gammaH[3]  = {g1, g2, g3};
    const Double_t s[3] = {s1, s2, s3};

    for (int i = 0; i < 3; i++)
    {
        Double_t dgamma = gammaL[i] - gammaH[i];

        // Numerator factor: E-dependent part
        Double_t num = 1.0 + TMath::Power(E / Ei[i], s[i]);

        // Denominator factor: E0-dependent part (normalization)
        Double_t den = 1.0 + TMath::Power(E0 / Ei[i], s[i]);

        // Add this break
        flux *= TMath::Power(num / den, dgamma / s[i]);
    }
    flux *= TMath::Power(E, 2.6);

    return flux;
}

Double_t brokenmod3_Wknee(Double_t *x, Double_t *par)
{
    const Double_t E  = x[0];
    if (E <= 0) return 0.0;

    // Parameters
    const Double_t J0 = par[0];

    // Gamma values before/after each break
    const Double_t g0 = par[1];
    const Double_t g1 = par[2];
    const Double_t g2 = par[3];
    const Double_t g3 = par[4];

    // Break energies
    const Double_t E1 = par[5];
    const Double_t E2 = par[6];
    const Double_t E3 = par[7];

    // Smoothness
    const Double_t s1  = par[8];
    const Double_t s2  = par[9];
    const Double_t s3  = par[10];

    //  FIXED PHYSICAL BREAK (p+He knee)
    // ------------------------------
    const Double_t g4 = 3.1;
    
    const Double_t Ep  = 3e6;   // 3 PeV in GeV
    const Double_t EHe = 6e6;   // 6 PeV in GeV
    //const Double_t E4  = 3.57e6;//TMath::Power(Ep * Ep* Ep * EHe, 1./4.);   // mid-knee ~ sqrt(3 PeV * 6 PeV)
    const Double_t E4  = par[11];

    const Double_t s4 = 5.;//0.2;   // fixed smoothing for final knee
    // ------------------------------

    // Reference energy
    const Double_t E0 = 1e2;

    // Main pre-break power law
    Double_t flux = J0 * TMath::Power(E / E0, -g0);

    // Array-friendly structure
    const Double_t Ei[4]      = {E1, E2, E3, E4};
    const Double_t gammaL[4]  = {g0, g1, g2, g3};
    const Double_t gammaH[4]  = {g1, g2, g3, g4};
    const Double_t s[4] = {s1, s2, s3, s4};

    for (int i = 0; i < 4; i++)
    {
        Double_t dgamma = gammaL[i] - gammaH[i];

        // Numerator factor: E-dependent part
        Double_t num = 1.0 + TMath::Power(E / Ei[i], s[i]);

        // Denominator factor: E0-dependent part (normalization)
        Double_t den = 1.0 + TMath::Power(E0 / Ei[i], s[i]);

        // Add this break
        flux *= TMath::Power(num / den, dgamma / s[i]);
    }
    flux *= TMath::Power(E, 2.6);

    return flux;
}

graphAsymErr read_file(std::string file_name, bool is_dampe){
    ifstream file(file_name);
    std::string line, first, h_title, val;
    std::string col_names[11] = {"Qty", "<E>", "Elo", "Eup", "y", "ystat_lo", "ystat_up", "ysyst_lo", "ysyst_up", "yerrtot_lo", "yerrtot_up"};
    int nline = -1;
    double line_values[10];/*E, Elo, Eup, y, ystat_lo, ystat_up, ysyst_lo, ysyst_up, yerrtot_lo, yerrtot_up*/
    std::vector<Double_t> fx/*<E>*/, fy/*y*/, felx/*<E>-Elo*/, fely/*yerrtot_lo*/, fehx/*Eup-<E>*/, fehy/*yerrtot_up*/;
    Double_t fx_min{INFINITY}, fx_max{0.};
    while(std::getline(file, line)){
        nline++;
        istringstream iss(line);
        iss>>first;
        if (first.substr(0,1) == "#"){
            if (nline == 0){
                h_title=line;
                h_title.erase(0, 1);
                cout<<h_title<<":\t";
            } else if ((nline == 1 && is_dampe) || (nline == 2 && !is_dampe)){ //check on the order of columns
                for(int k = 0; k<10; k++){
                    iss>>val;
                    if(val != col_names[k]){
                        cout<<file_name<<'\n';
                        cout<<"ERROR: columns have a different order:\n"<<col_names[k]<<'\t'<<val<<'\n';
                    }
                }
            } 
            continue;
        } else {
            /*if(nline>38 && is_dampe) {
            fx.push_back(0.);
            fy.push_back(0.);
            felx.push_back(0.);
            fely.push_back(0.);
            fehx.push_back(0.);
            fehy.push_back(0.);
            } else{*/
            for(int k = 0; k<10; k++){
                iss>>line_values[k];
            }
            fx.push_back(line_values[0]);
            fy.push_back(line_values[3]);
            felx.push_back(line_values[0] - line_values[1]);
            fely.push_back(line_values[8]);
            fehx.push_back(line_values[2] - line_values[0]);
            fehy.push_back(line_values[9]);
            if(line_values[0] < fx_min){fx_min = line_values[0];}
            if(line_values[0] > fx_max){fx_max = line_values[0];}
            //}
            //cout<<nline<<": x = "<<fx.back()<<"\ty = "<<fy.back()<<'\n';
        }

    }
    if(is_dampe){cout<<nline-1<<" pts\n";} 
    else {cout<<nline-2<<" pts\n";}

    graphAsymErr graph_vecs;
    graph_vecs.graph_title = h_title.c_str();
    graph_vecs.gx = fx;
    graph_vecs.gy = fy;
    graph_vecs.gelx = felx;
    graph_vecs.gely = fely;
    graph_vecs.gehx = fehx;
    graph_vecs.gehy = fehy;
    graph_vecs.gx_min = fx_min;
    graph_vecs.gx_max = fx_max;

    return graph_vecs;

}

graphAsymErr get_pow_data_withElw(graphAsymErr data, double pow, double prior_pow){
    cout<<"\tUsing get_pow_data_withElw\n";
    int n_pts = data.gx.size();
    graphAsymErr data_pow;
    data_pow.graph_title = data.graph_title;
    data_pow.gx = data.gx;
    std::vector<Double_t> fy, ely, ehy;
    for(int i = 0; i< n_pts; i++){
        double erry, erry2, delta_log_x;
        double e1, e2, xlw; 
        e1 = data.gx.at(i) - data.gelx.at(i);
        e2 = data.gx.at(i) + data.gehx.at(i);
        xlw = compute_xlw(e1, e2, prior_pow);
        //cout<<"e1 = "<<e1<<"\te2 = "<<e2<<"\txlw = "<<xlw<<'\n';
        fy.push_back(data.gy.at(i) * TMath::Power(xlw, pow));
        erry = data.gely.at(i) * TMath::Power(xlw, pow);
        erry2 = data.gehy.at(i) * TMath::Power(xlw, pow);
        ely.push_back(erry);
        ehy.push_back(erry2);
        /*cout<<"data:\t";
        cout<<data.gx.at(i)<<'\t'<<data.gy.at(i)<<'\t'<<data.gelx.at(i)<<'\t'<<data.gehx.at(i)<<'\t'<<data.gely.at(i)<<'\t'<<data.gehy.at(i)<<'\n';
        cout<<"pow:\t";
        cout<<data.gx.at(i)<<'\t'<<fy.at(i)<<'\t'<<data.gelx.at(i)<<'\t'<<data.gehx.at(i)<<'\t'<<ely.at(i)<<'\t'<<ehy.at(i)<<'\n';
            */
        //cout<<"erry = "<<erry<<"\terry2 = "<<erry2<<"\t erry+erry2 = "<<erry+erry2<<"\t sqrt(erry*erry +erry2*erry2) = "<<std::sqrt(erry*erry + erry2*erry2)<<'\n';
    }
    data_pow.gy = fy;
    data_pow.gelx = data.gelx;
    data_pow.gehx = data.gehx;
    data_pow.gely = ely;
    data_pow.gehy = ehy;
    data_pow.gx_min = data.gx_min;
    data_pow.gx_max = data.gx_max;

    return data_pow;
}

double compute_xlw(double e1, double e2, double pow){
    double xlw = TMath::Power(((TMath::Power(e2, 1.-pow)) - (TMath::Power(e1, 1.-pow))) / ((1.-pow)*(e2-e1)),-1./pow);
    return xlw;
}
