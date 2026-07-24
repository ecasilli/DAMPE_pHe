
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <vector>

void fillSpectrum(
TChain *ch,
double perFileWeight,
bool applyGeoCorr,
double GeoCorr,
int noe,
int nbd,
Double_t *Ebin,
double *WBin,

TCut &bgo_acc,
TCut &cut06,
TCut &cut00,
TCut &cut01,
TCut &cut05,
TCut &SpCut,
TCut &Cut_PHe_MCCor6,

TH2D *hall,
TH2D *h06,
TH2D *h00,
TH2D *h01,
TH2D *h05,
TH2D *hSp,
TH2D *hFinal)
{

Long64_t nentries=ch->GetEntries();

TTreeFormula fMCE("fMCE","MC_EnergyT",ch);
TTreeFormula fBGO("fBGO","BGO_EnergyG_QuenchSatCorr_ML_ions_v3",ch);

TTreeFormula fBgoAcc("fBgoAcc",bgo_acc.GetTitle(),ch);
TTreeFormula fCut06("fCut06",cut06.GetTitle(),ch);
TTreeFormula fCut00("fCut00",cut00.GetTitle(),ch);
TTreeFormula fCut01("fCut01",cut01.GetTitle(),ch);
TTreeFormula fCut05("fCut05",cut05.GetTitle(),ch);
TTreeFormula fSpCut("fSpCut",SpCut.GetTitle(),ch);
TTreeFormula fChargeCut("fChargeCut",Cut_PHe_MCCor6.GetTitle(),ch);

double gc=1.0;
if(applyGeoCorr) gc=GeoCorr;

for(Long64_t ie=0;ie<nentries;ie++){

    if(ie%100000==0)
        cout<<"Processing "<<ie<<" / "<<nentries<<endl;

    ch->GetEntry(ie);

    double E=fMCE.EvalInstance();

    if(E<Ebin[0]) continue;
    if(E>=Ebin[noe]) continue;

    int j=(int)((TMath::Log10(E)-TMath::Log10((double)Ebin[0]))*nbd);

    if(j<0) j=0;
    if(j>=noe) j=noe-1;

    double w=WBin[j]*pow(E,-1.7)*perFileWeight*gc;

    double reco=fBGO.EvalInstance();

    hall->Fill(E,reco,w);

    if(!fBgoAcc.EvalInstance()) continue;

    if(!fCut06.EvalInstance()) continue;
    h06->Fill(E,reco,w);

    if(!fCut00.EvalInstance()) continue;
    h00->Fill(E,reco,w);

    if(!fCut01.EvalInstance()) continue;
    h01->Fill(E,reco,w);

    if(!fCut05.EvalInstance()) continue;
    h05->Fill(E,reco,w);

    if(!fSpCut.EvalInstance()) continue;
    hSp->Fill(E,reco,w);

    if(!fChargeCut.EvalInstance()) continue;

    hFinal->Fill(E,reco,w);
}

}

void Load_MC_He_p_5PeV_spectrum_fast(){


// =============================
// ----- Rilevamento automatico del sito ----

TString basePath;
TString hostname = gSystem->HostName();

if (hostname.Contains("cnaf")) {
    basePath = "/storage/gpfs_data/dampe/users/SKIM_2026_pHe";
} else if (hostname.Contains("le.infn.it")) {
    basePath = "/nfs/argo/dampe/SKIM_2026_pHe";
} else {
    cout << "ERROR: hostname non riconosciuto: " << hostname << endl;
    cout << "  Aggiungere il path corretto per questo sito." << endl;
    return;
}

cout << "Hostname: " << hostname << endl;
cout << "Base path: " << basePath << endl;

// =============================

const int nsetHe = 7;
const int nsetP  = 6;

std::vector<TString> filesP;
filesP.push_back(basePath + "/PROTON/Proton_10GeV_100GeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_100GeV_1TeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_1TeV_10TeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_10TeV_100TeV_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_100TeV_1PeV_EPOSLHC_FTFP_merged.root");
filesP.push_back(basePath + "/PROTON/Proton_1PeV_5PeV_EPOSLHC_FTFP_merged.root");

std::vector<TString> filesHe;
filesHe.push_back(basePath + "/HELIUM/He4_10GeV_100GeV_FTFP_BGO_Quenching_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_100GeV_1TeV_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_1TeV_10TeV_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_10TeV_100TeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_100TeV_500TeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_500TeV_1PeV_EPOSLHC_FTFP_merged.root");
filesHe.push_back(basePath + "/HELIUM/He4_1PeV_5PeV_EPOSLHC_FTFP_merged.root");

const int nbd  = 6;   // 6 bin per decade
const int ndec = 6;
Double_t arg1 = 1./double(nbd);
const int noe = nbd*ndec;
cout << "Number of energy bin " << noe << endl;
Double_t e0 = 10.;

Double_t Ebin[noe+1];
Ebin[0] = e0;
for (int i=1; i<noe+1; i++) Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);

// ---- Peso di normalizzazione per bin (era "weights[i]" ricalcolato ad ogni Draw) ----
// Dipende solo dagli estremi del bin: lo calcoliamo UNA volta sola qui,
// invece di ricostruirlo come TString/TCut 36 volte per ciascuno dei 13 file.
double WBin[noe];
for (int i=0; i<noe; i++){
    double E0 = Ebin[i], E1 = Ebin[i+1];
    WBin[i] = 1.7*TMath::Log(E1/E0) / (TMath::Power(E0,-1.7) - TMath::Power(E1,-1.7));
}

// =============================

TChain *sk_he[nsetHe];
TChain *sk_p[nsetP];
for (int i=0; i<nsetP;  i++) sk_p[i]  = new TChain("newtree");
for (int i=0; i<nsetHe; i++) sk_he[i] = new TChain("newtree");

// I pesi per-file sono costanti numeriche pure (nessun branch coinvolto):
// li valutiamo direttamente come double invece di passare per TCut/TTreeFormula.
double wP[nsetP] = {
    (1./1740700000.)*TMath::Log(10.),   // 10-100 GeV
    (1./519070000.)*TMath::Log(10.),    // 100GeV-1TeV
    (1./201032000.)*TMath::Log(10.),    // 1TeV-10TeV
    (1./138744900.)*TMath::Log(10.),    // 10TeV-100TeV
    (1./19441600.)*TMath::Log(10.),     // 100TeV-1PeV
    (1./5094200.)*TMath::Log(5.)        // 1PeV-5PeV
};

double wHe[nsetHe] = {
    (1./20380000.)*TMath::Log(10.),     // 10-100 GeV
    (1./19556000.)*TMath::Log(10.),     // 100GeV-1TeV
    (1./11346000.)*TMath::Log(10.),     // 1TeV-10TeV
    (1./17694000.)*TMath::Log(10.),     // 10TeV-100TeV
    (1./10305100.)*TMath::Log(5.),      // 100TeV-500TeV
    (1./10281395.)*TMath::Log(2.),      // 500TeV-1PeV
    (1./5064900.)*TMath::Log(5.)        // 1PeV-5PeV
};

double GeoCorr = (2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38);

// `````````````````````````````````````````````````````````````` PROTON
for (int i=0; i<nsetP; i++) sk_p[i]->Add(filesP[i]);

cout<<" "<<endl;
cout<<"  10 GeV - 100 GeV ...sk_p[0]-> "<<sk_p[0]->GetEntries()<<endl;
cout<<" 100 GeV - 1 TeV   ...sk_p[1]-> "<<sk_p[1]->GetEntries()<<endl;
cout<<"   1 TeV - 10 TeV  ...sk_p[2]-> "<<sk_p[2]->GetEntries()<<endl;
cout<<"  10 TeV - 100 TeV ...sk_p[3]-> "<<sk_p[3]->GetEntries()<<endl;
cout<<" 100 TeV - 1 PeV   ...sk_p[4]-> "<<sk_p[4]->GetEntries()<<endl;
cout<<"   1 PeV - 5 PeV   ...sk_p[5]-> "<<sk_p[5]->GetEntries()<<endl;
cout<<" "<<endl;

// `````````````````````````````````````````````````````````````` HELIUM
for (int i=0; i<nsetHe; i++) sk_he[i]->Add(filesHe[i]);

cout<<" "<<endl;
cout<<"  10 GeV - 100 GeV ...sk_he[0]-> "<<sk_he[0]->GetEntries()<<endl;
cout<<" 100 GeV - 1 TeV   ...sk_he[1]-> "<<sk_he[1]->GetEntries()<<endl;
cout<<"   1 TeV - 10 TeV  ...sk_he[2]-> "<<sk_he[2]->GetEntries()<<endl;
cout<<"  10 TeV - 100 TeV ...sk_he[3]-> "<<sk_he[3]->GetEntries()<<endl;
cout<<" 100 TeV - 500 TeV ...sk_he[4]-> "<<sk_he[4]->GetEntries()<<endl;
cout<<" 500 TeV - 1 PeV   ...sk_he[5]-> "<<sk_he[5]->GetEntries()<<endl;
cout<<"   1 PeV - 5 PeV   ...sk_he[6]-> "<<sk_he[6]->GetEntries()<<endl;
cout<<" "<<endl;

// ------------------------------------------------------------------------------
//    LISTA TAGLI (invariata rispetto all'originale)

TCut bgo_valid = "!( (BGO_slopeXZ_analy==0 && BGO_interceptXZ_analy==0) || (BGO_slopeYZ_analy==0 && BGO_interceptYZ_analy==0) )";
TCut bgo_acceptance =
    "fabs(BGO_interceptXZ_analy + 448.*BGO_slopeXZ_analy) < 280. && "
    "fabs(BGO_interceptYZ_analy + 448.*BGO_slopeYZ_analy) < 280. && "
    "fabs(BGO_interceptXZ_analy + 46.*BGO_slopeXZ_analy)  < 280. && "
    "fabs(BGO_interceptYZ_analy + 46.*BGO_slopeYZ_analy)  < 280.";
TCut bgo_acc = bgo_valid*bgo_acceptance;

TCut Trig_HEP = "BGO_HET>0.";
TCut cc204s   = "(BGO_EnergyG_QuenchSatCorr_ML_ions_v3>20.)";
TCut cut00    = cc204s*Trig_HEP;
TCut cut01    = "(PSD_ChargeY0>0.0 || PSD_ChargeY1>0.0) && (PSD_ChargeX0>0.0 || PSD_ChargeX1>0.0)";
TCut cut05    = "TMath::Abs(BGO_cbgomax[0]-BGO_cbgostk[0])<30. && TMath::Abs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.";
TCut cut06    = "fabs(STKtrack_to_PSD_topY)<400. && fabs(STKtrack_to_PSD_topX)<400.";
TCut SpCut    = "BGO_xtr > 12.";

// ---- selezione di carica (identica all'originale) ----
//TString PSDcharge = "( (((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)*PSD_ChargeY0)+(((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)*PSD_ChargeY1)+(((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)*PSD_ChargeX0)+(((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.)*PSD_ChargeX1) ) / ( ((TMath::Sign(1.,PSD_ChargeY0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeY1)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX0)+1.)/2.)+((TMath::Sign(1.,PSD_ChargeX1)+1.)/2.) )";
TString PSDcharge = "PSD_Global_Charge";

TString HeMPVf    = "(1.88586+(0.185264*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0867115*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0200701*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.000810826*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )";
TString HeWidthf  = "(-0.0824082+(0.201315*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0989824*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0218456*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00136192*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeGSigmaf = "(0.0396516)";

TString HeMcMPVf    = "(1.24779+(1.25798*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.698468*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.169728*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.0123542*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcWidthf  = "(-0.207322+(0.351724*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.173554*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0388396*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00242265*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcGSigmaf = "(0.0626831)";

TString HeMcMPVf_Cor    = "(1.94878+(0.0790806*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0187037*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.00120807*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00112927*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcWidthf_Cor  = "(0.00973057+(0.0360571*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.00109536*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.00185953*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.000683651*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString HeMcGSigmaf_Cor = "(0.0389501)";

TString HeMcFSig    = "sqrt("+HeMcWidthf+"**2+"+HeMcGSigmaf+"**2)";
TString HeMcFSigCor = "sqrt("+HeMcWidthf_Cor+"**2+"+HeMcGSigmaf_Cor+"**2)";
TString HeFSig      = "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

TString PMPVf    = "(0.993603+(-0.00402039*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(0.0408375*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(-0.0227868*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00494042*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )";
TString PWidthf  = "(-0.0611545+(0.140328*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0683705*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0143596*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.000527989*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PGSigmaf = "(5.5475e-08)";

TString PMcMPVf    = "(0.526726+(0.806415*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.488731*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.129689*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00965891*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcWidthf  = "(-0.151343+(0.289688*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.154242*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0383406*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.00270484*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcGSigmaf = "(0.000127013)";

TString PMcMPVf_Cor    = "(0.931721+(0.1114*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.0370493*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.000118024*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(0.00248452*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcWidthf_Cor  = "(-0.116835+(0.233937*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))+(-0.125023*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)+(0.0294152*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(-0.0020157*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4))";
TString PMcGSigmaf_Cor = "(9.39232e-05)";

TString PMcFSig    = "sqrt("+PMcWidthf+"**2+"+PMcGSigmaf+"**2)";
TString PMcFSigCor = "sqrt("+PMcWidthf_Cor+"**2+"+PMcGSigmaf_Cor+"**2)";
TString PFSig      = "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";

TString HeliumCor = "(("+PSDcharge+"-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+")";
TString ProtonCor = "(("+PSDcharge+"-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+")";

TString p_low_mcCor    = "("+PMcMPVf_Cor+"-"+ProtonCor+")<(2.5*"+PMcFSigCor+")";
TString he_high_mcCor6 = "("+HeliumCor+"-"+HeMcMPVf_Cor+")<(6*"+HeMcFSigCor+")";

TCut PLow_MCCor      = p_low_mcCor;
TCut HeHigh_MCCor6   = he_high_mcCor6;
TCut Cut_PHe_MCCor6  = PLow_MCCor*HeHigh_MCCor6;

TString STKcharge = "((((TMath::Sign(1.,STK_chargeY[0])+1.)/2.*STK_chargeY[0]+(TMath::Sign(1.,STK_chargeX[0])+1.)/2.*STK_chargeX[0])/((TMath::Sign(1.,STK_chargeY[0])+1.)/2.+(TMath::Sign(1.,STK_chargeX[0])+1.)/2.)))";

TString stk_charge_cut = "(" + STKcharge + " < 400.)";
TCut STKChargeSel = stk_charge_cut;

TString ChargeSelection =
    "((STK_vertexPrediction < 0.7) && (" + TString(Cut_PHe_MCCor6.GetTitle()) + ")) || "
    "((STK_vertexPrediction >= 0.7) && (" + STKcharge + " < 400.))";

TCut Cut_PHe_New = ChargeSelection;

// ---------------  TAGLI TOTALI CUMULATIVI (identici all'originale)
// ctot_06 = bgo_acc*cut06
// ctot_00 = ctot_06*cut00
// ctot_01 = ctot_00*cut01
// ctot_05 = ctot_01*cut05
// ctot_SpCut = ctot_05*SpCut
// ctot = ctot_SpCut*Cut_PHe_MCCor6

TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_unfolding_6binperdecade_2e5sigmaLow_6sigmaUp_noCut02_wPHe_PSDprogr_STKvert_fast.root","RECREATE");

TH2D *h2Ntrig_wgt_v3_all   = new TH2D("h2Ntrig_wgt_v3_all",   "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2D *h2Ntrig_wgt_v3_cut06 = new TH2D("h2Ntrig_wgt_v3_cut06", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2D *h2Ntrig_wgt_v3_cut00 = new TH2D("h2Ntrig_wgt_v3_cut00", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2D *h2Ntrig_wgt_v3_cut01 = new TH2D("h2Ntrig_wgt_v3_cut01", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2D *h2Ntrig_wgt_v3_cut05 = new TH2D("h2Ntrig_wgt_v3_cut05", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2D *h2Ntrig_wgt_v3_SpCut = new TH2D("h2Ntrig_wgt_v3_SpCut", "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);
TH2D *h2Ntrig_wgt_v3       = new TH2D("h2Ntrig_wgt_v3",       "Ntrig(Eo,Et)", noe, Ebin, noe, Ebin);



for(int i=0;i<nsetHe;i++){

    cout<<"Processing He sample "<<i<<endl;

    fillSpectrum(
        sk_he[i],
        wHe[i],
        true,
        GeoCorr,
        noe,
        nbd,
        Ebin,
        WBin,

        bgo_acc,
        cut06,
        cut00,
        cut01,
        cut05,
        SpCut,
        Cut_PHe_New,

        h2Ntrig_wgt_v3_all,
        h2Ntrig_wgt_v3_cut06,
        h2Ntrig_wgt_v3_cut00,
        h2Ntrig_wgt_v3_cut01,
        h2Ntrig_wgt_v3_cut05,
        h2Ntrig_wgt_v3_SpCut,
        h2Ntrig_wgt_v3
    );

    cout<<"Finished He "<<i<<endl;

}

for(int i=0;i<nsetP;i++){

    cout<<"Processing Proton sample "<<i<<endl;

    fillSpectrum(
        sk_p[i],
        wP[i],
        false,
        GeoCorr,
        noe,
        nbd,
        Ebin,
        WBin,

        bgo_acc,
        cut06,
        cut00,
        cut01,
        cut05,
        SpCut,
        Cut_PHe_New,

        h2Ntrig_wgt_v3_all,
        h2Ntrig_wgt_v3_cut06,
        h2Ntrig_wgt_v3_cut00,
        h2Ntrig_wgt_v3_cut01,
        h2Ntrig_wgt_v3_cut05,
        h2Ntrig_wgt_v3_SpCut,
        h2Ntrig_wgt_v3
    );

    cout<<"Finished Proton "<<i<<endl;

}

fout->cd();
h2Ntrig_wgt_v3_all->Write();
h2Ntrig_wgt_v3_cut06->Write();
h2Ntrig_wgt_v3_cut00->Write();
h2Ntrig_wgt_v3_cut01->Write();
h2Ntrig_wgt_v3_cut05->Write();
h2Ntrig_wgt_v3_SpCut->Write();
h2Ntrig_wgt_v3->Write();
fout->Close();

cout<<"End script."<<endl;
}