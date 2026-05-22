// run with root -l ChargeSelectionGraph.C+

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TStyle.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"

#include "TLatex.h"

#include <ctime>
#include <iomanip>

void DrawFitParameters(TF1* f, double x, double y, int color, const char* label, int npar)
{
    TLatex latex;

    latex.SetNDC();
    latex.SetTextSize(0.022);
    latex.SetTextFont(42);
    latex.SetTextColor(color);

    latex.DrawLatex(x, y, Form("%s", label));

    for(int i=0; i<npar; i++)
    {
        latex.DrawLatex(
            x,
            y-0.04*(i+1),
            Form("p%d = %.3e #pm %.3e",
                 i,
                 f->GetParameter(i),
                 f->GetParError(i))
        );
    }

    latex.DrawLatex(
        x,
        y-0.04*(npar+1),
        Form("#chi^{2}/NDF = %.1f / %d",
             f->GetChisquare(),
             f->GetNDF())
    );
}

void PrintTimestamp(const char* scriptName)
{
    time_t now = time(0);
    tm *ltm = localtime(&now);

    cout << "====================================" << endl;
    cout << "Running " << scriptName << endl;

    cout << "Date: "
         << 1900 + ltm->tm_year << "-"
         << setw(2) << setfill('0') << 1 + ltm->tm_mon << "-"
         << setw(2) << setfill('0') << ltm->tm_mday
         << "  "
         << setw(2) << setfill('0') << ltm->tm_hour << ":"
         << setw(2) << setfill('0') << ltm->tm_min << ":"
         << setw(2) << setfill('0') << ltm->tm_sec
         << endl;

    cout << "====================================" << endl;
}

std::string GetTimestamp()
{
    time_t now = time(0);

    char buf[64];

    strftime(buf, sizeof(buf),
             "%Y%m%d_%H%M%S",
             localtime(&now));

    return std::string(buf);
}

void ChargeSelectionGraph()
{
	PrintTimestamp("ChargeSelectionGraph.C");
	std::string timestamp = GetTimestamp();
	//std::string ParticleTag = "H";
	std::string ParticleTag = "He";

		gROOT->Reset();
		gStyle->SetOptFit(0);
		gStyle->SetOptStat(0);
		
		double EnergyI[20] = {0.};
		double EnergyF[20] = {0.};
		
		double MPVH[20] = {0.};
		double E_MPVH[20] = {0.};
		double SigmaLH[20] = {0.};
		double E_SigmaLH[20] = {0.};
		double SigmaGH[20] = {0.};
		double E_SigmaGH[20] = {0.};
		
		double MPVHMC[20] = {0.};
		double E_MPVHMC[20] = {0.};
		double SigmaLHMC[20] = {0.};
		double E_SigmaLHMC[20] = {0.};
		double SigmaGHMC[20] = {0.};
		double E_SigmaGHMC[20] = {0.};
	
		std::ifstream in;	
		in.open("LangausResults_He_corrected_v2.txt");
		//in.open("LangausResults_corrected_v2.txt");
		
		int i;
		
		double BGOEnergy[20] = {0.};
		double E_BGOEnergy[20] = {0.};
		
		double SigmaTotH[20] = {0.};
		double SigmaTotHMC[20] = {0.};
		double E_SigmaTotH[20] = {0.};
		double E_SigmaTotHMC[20] = {0.};
		
		
		for(int i=0; i<15; i++)
		{
			std::string dummy;
			std::getline(in, dummy);
			
			in >> EnergyI[i] >> EnergyF[i] >> MPVHMC[i] >> E_MPVHMC[i] >> SigmaLHMC[i] >> E_SigmaLHMC[i] >> SigmaGHMC[i] >> E_SigmaGHMC[i] >> MPVH[i] >> E_MPVH[i] >> SigmaLH[i] >> E_SigmaLH[i] >> SigmaGH[i] >> E_SigmaGH[i];
			
			//std::cout << EnergyI[i] << "\t" << EnergyF[i] << "\t" << MPVH[i] <<"\t" << SigmaGH[i] << "\t" << MPVHMC[i] << std::endl;
			
			BGOEnergy[i] = (EnergyI[i] + EnergyF[i])/2;
//			BGOEnergy[i] = pow(((pow(EnergyF[i],(2.7)) - pow(EnergyI[i],(2.7)))/((EnergyF[i] - EnergyI[i])*(1 - (2.7)))),(-1/2.7));
			E_BGOEnergy[i] = 0.01*(i+1);
//			E_BGOEnergy[i] = (EnergyF[i] - EnergyI[i])/2.;
			
			//std::cout << E_BGOEnergy[i] << std::endl;
			
			SigmaTotH[i] = sqrt((SigmaLH[i]*SigmaLH[i])+(SigmaGH[i]*SigmaGH[i]));
			SigmaTotHMC[i] = sqrt((SigmaLHMC[i]*SigmaLHMC[i])+(SigmaGHMC[i]*SigmaGHMC[i]));
			
			//std::cout << SigmaLHMC[i] << "\t" << SigmaGHMC[i] << "\t" << SigmaTotHMC[i] << std::endl;
			
			E_SigmaTotH[i] = sqrt((((SigmaLH[i]*SigmaLH[i])/((SigmaLH[i]*SigmaLH[i])+(SigmaGH[i]*SigmaGH[i])))*(E_SigmaLH[i]*E_SigmaLH[i]))+(((SigmaGH[i]*SigmaGH[i])/((SigmaLH[i]*SigmaLH[i])+(SigmaGH[i]*SigmaGH[i])))*(E_SigmaGH[i]*E_SigmaGH[i])));
			E_SigmaTotHMC[i] = sqrt((((SigmaLHMC[i]*SigmaLHMC[i])/((SigmaLHMC[i]*SigmaLHMC[i])+(SigmaGHMC[i]*SigmaGHMC[i])))*(E_SigmaLHMC[i]*E_SigmaLHMC[i]))+(((SigmaGHMC[i]*SigmaGHMC[i])/((SigmaLHMC[i]*SigmaLHMC[i])+(SigmaGHMC[i]*SigmaGHMC[i])))*(E_SigmaGHMC[i]*E_SigmaGHMC[i])));
			
			//std::cout << SigmaTotH[i] << "\t" << SigmaTotHMC[i] << std::endl;

		}	
		

		//gStyle->SetOptFit(222);		

cout << "\n * * * * * * *    MPV    * * * * * * * \n" << endl;
	
	//	TF1 *ChargeFitFuncH  = new TF1("ChargeFitFuncH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,41000.0);	//39820.0
		//TF1 *ChargeFitFuncH  = new TF1("ChargeFitFuncH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",25.0,100000.0);
		TF1 *ChargeFitFuncH  = new TF1("ChargeFitFuncH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,100000.0);
	//	TF1 *ChargeFitFuncH  = new TF1("ChargeFitFuncH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,15000.0);

		//ChargeFitFuncH->SetParameters(2.236, -0.6383, 0.622, -0.2396, 0.03738);	
		ChargeFitFuncH->SetLineColor(kBlue+1);
		TCanvas *C1 = new TCanvas("C1","C1",1000,800);
		
		TGraphErrors *ChargeSelectionMPVH = new TGraphErrors(14, BGOEnergy, MPVH, E_BGOEnergy, E_MPVH);
		
		C1->SetLogx();
		C1->SetTickx();
		C1->SetTicky();
		
		ChargeSelectionMPVH->SetMarkerStyle(20);
		ChargeSelectionMPVH->SetMarkerSize(1.2);
		ChargeSelectionMPVH->SetMarkerColor(kBlue+1);
		ChargeSelectionMPVH->SetLineColor(kBlue+1);
		ChargeSelectionMPVH->SetTitle(Form(" ;BGO Energy (GeV); MPV %s ", ParticleTag.c_str()));
		ChargeSelectionMPVH->SetFillStyle(0);  
		ChargeSelectionMPVH->Fit("ChargeFitFuncH","R");	
		//gStyle->SetOptFit(222);
		ChargeSelectionMPVH->Draw("AP");	

		//TF1 *ChargeFitFuncHMC  = new TF1("ChargeFitFuncHMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",25.0,100000.0);
		TF1 *ChargeFitFuncHMC  = new TF1("ChargeFitFuncHMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,100000.0);
		//TF1 *ChargeFitFuncHMC  = new TF1("ChargeFitFuncHMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,25120.0);
		//ChargeFitFuncHMC->SetParameters(2.236, -0.6383, 0.622, -0.2396, 0.03738);
		ChargeFitFuncHMC->SetLineColor(kRed+1);
			
	//	TCanvas *C2 = new TCanvas("C2","C2",1000,800);
		
		TGraphErrors *ChargeSelectionMPVHMC = new TGraphErrors(14, BGOEnergy, MPVHMC, E_BGOEnergy, E_MPVHMC);
		/*
		C2->SetLogx();
		C2->SetTickx();
		C2->SetTicky();
		*/
		ChargeSelectionMPVHMC->SetMarkerStyle(20);
		ChargeSelectionMPVHMC->SetMarkerSize(1.2);
		ChargeSelectionMPVHMC->SetMarkerColor(kRed+1);
		ChargeSelectionMPVHMC->SetLineColor(kRed+1);
		ChargeSelectionMPVHMC->SetTitle(Form(" ;BGO Energy (GeV); MPV %s ", ParticleTag.c_str()));
		ChargeSelectionMPVHMC->SetFillStyle(0);  
		ChargeSelectionMPVHMC->Fit("ChargeFitFuncHMC","R"); 
		//gStyle->SetOptFit(222);		
		ChargeSelectionMPVHMC->Draw("AP");

		cout << "TString PMPVf= \"(" << ChargeFitFuncH->GetParameter(0) << "+(" << ChargeFitFuncH->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))" << "+(" << ChargeFitFuncH->GetParameter(2) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)" << "+(" << ChargeFitFuncH->GetParameter(3) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(" << ChargeFitFuncH->GetParameter(4) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4) )\";" << endl;
    	cout << "TString PMcMPVf= \"(" << ChargeFitFuncHMC->GetParameter(0) << "+(" << ChargeFitFuncHMC->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))" << "+(" << ChargeFitFuncHMC->GetParameter(2) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)" << "+(" << ChargeFitFuncHMC->GetParameter(3) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(" << ChargeFitFuncHMC->GetParameter(4) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4)" << ")\";" << endl;
    
			
	TMultiGraph *mg = new TMultiGraph();	
	mg->SetTitle(Form(" ;BGO Energy (GeV);MPV %s ", ParticleTag.c_str()));
    	mg->Add(ChargeSelectionMPVH);
	mg->Add(ChargeSelectionMPVHMC);
	mg->Draw("ap");

mg->GetXaxis()->SetLimits(20, 33000.);
DrawFitParameters(ChargeFitFuncH,  0.15, 0.65, kBlue+1, "Data", 5);
DrawFitParameters(ChargeFitFuncHMC, 0.4, 0.65, kRed+1,  "MC",   5);

TLegend* leg = new TLegend(0.1352705,0.6745306,0.3857715,0.8748195); //0.1352705,0.6745306,0.3857715,0.8748195,NULL,"brNDC"
leg->SetFillColor(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.04);
leg->AddEntry(ChargeSelectionMPVHMC,"MC data","ep");
leg->AddEntry(ChargeSelectionMPVH,"Flight data","ep");

leg->Draw();

cout << "\n * * * * * * *    Width    * * * * * * * \n" << endl;

/*	
		TF1 *ChargeFitFuncSigmaH  = new TF1("ChargeFitFuncSigmaH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,25120.0);
		ChargeFitFuncSigmaH->SetParameters( -0.01816, 0.1312, 0.1498, -0.1056, 0.01933);
		
		TCanvas *C3 = new TCanvas("C3","C3",1000,800);
		
		TGraphErrors *SigmaH = new TGraphErrors(15, BGOEnergy, SigmaTotH, E_BGOEnergy, E_SigmaTotH);
		
		C3->SetLogx();
		
		SigmaH->SetMarkerStyle(20);
		SigmaH->SetMarkerSize(0.8);
		SigmaH->SetMarkerColor(2);
		SigmaH->SetLineColor(2);
		SigmaH->SetTitle(" ;BGO Energy (GeV); #sigma_{H} ");
		SigmaH->SetFillStyle(0);  
		SigmaH->Fit("ChargeFitFuncSigmaH","R"); 
	//	SigmaH->Draw("AP");
		

		TF1 *ChargeFitFuncSigmaHMC  = new TF1("ChargeFitFuncSigmaHMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,25120.0);
		ChargeFitFuncSigmaHMC->SetParameters( -0.01816, 0.1312, 0.1498, -0.1056, 0.01933);		
	
		TCanvas *C4 = new TCanvas("C4","C4",1000,800);
		
		TGraphErrors *SigmaHMC = new TGraphErrors(15, BGOEnergy, SigmaTotHMC, E_BGOEnergy, E_SigmaTotHMC);
		
		C4->SetLogx();
		
		SigmaHMC->SetMarkerStyle(20);
		SigmaHMC->SetMarkerSize(0.8);
		SigmaHMC->SetMarkerColor(2);
		SigmaHMC->SetLineColor(2);
		SigmaHMC->SetTitle(" ;BGO Energy (GeV); #sigma_{HMC} ");
		SigmaHMC->SetFillStyle(0);  
		SigmaHMC->Fit("ChargeFitFuncSigmaHMC","R"); 		
		SigmaHMC->Draw("AP");

*/		
		//TF1 *ChargeFitFuncWidthH  = new TF1("ChargeFitFuncWidthH","[0]+[1]*log10(x)",25.0,100000.0);
//		TF1 *ChargeFitFuncWidthH  = new TF1("ChargeFitFuncWidthH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3",25.,100000.0);
//		TF1 *ChargeFitFuncWidthH  = new TF1("ChargeFitFuncWidthH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,39820.0);	
		TF1 *ChargeFitFuncWidthH  = new TF1("ChargeFitFuncWidthH","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,100000.0);
		//ChargeFitFuncWidthH->SetParameters(-1.25302, 2.11415, -1.07456, 0.225429, -0.0132164);
		ChargeFitFuncWidthH->SetLineColor(kBlue+1);
		
		TCanvas *C5 = new TCanvas("C5","C5",1000,800);
		
		TGraphErrors *WidthH = new TGraphErrors(14, BGOEnergy, SigmaLH, E_BGOEnergy, E_SigmaLH);
		
		C5->SetLogx();
		C5->SetTickx();
		C5->SetTicky();

		//gStyle->SetOptFit(222);		
		WidthH->SetMarkerStyle(20);
		WidthH->SetMarkerSize(1.2);
		WidthH->SetMarkerColor(kBlue+1);
		WidthH->SetLineColor(kBlue+1);
		WidthH->SetTitle(Form(" ;BGO Energy (GeV); Width %s ", ParticleTag.c_str()));
		WidthH->SetFillStyle(0);  
		WidthH->Fit("ChargeFitFuncWidthH","R"); 
		WidthH->Draw("AP");
		
		
		//TF1 *ChargeFitFuncWidthHMC  = new TF1("ChargeFitFuncWidthHMC","[0]+[1]*log10(x)",25.0,100000.0);
		TF1 *ChargeFitFuncWidthHMC  = new TF1("ChargeFitFuncWidthHMC","[0]+[1]*log10(x)+[2]*log10(x)**2",63.10,100000.0);
//		TF1 *ChargeFitFuncWidthHMC  = new TF1("ChargeFitFuncWidthHMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",63.10,100000.0);
//		TF1 *ChargeFitFuncWidthHMC  = new TF1("ChargeFitFuncWidthHMC","[0]+[1]*log10(x)+[2]*log10(x)**2+[3]*log10(x)**3+[4]*log10(x)**4",25.0,25120.0);
		//ChargeFitFuncWidthHMC->SetParameters(-1.25302, 2.11415, -1.07456, 0.225429, -0.0132164);		
		ChargeFitFuncWidthHMC->SetLineColor(kRed+1);
		
//		TCanvas *C6 = new TCanvas("C6","C6",1000,800);
		
		TGraphErrors *WidthHMC = new TGraphErrors(14, BGOEnergy, SigmaLHMC, E_BGOEnergy, E_SigmaLHMC);
/*		
		C6->SetLogx();
		C6->SetTickx();
		C6->SetTicky();
*/
		
		//gStyle->SetOptFit(011);
		WidthHMC->SetMarkerStyle(20);
		WidthHMC->SetMarkerSize(1.2);
		WidthHMC->SetMarkerColor(kRed+1);
		WidthHMC->SetLineColor(kRed+1);
		WidthHMC->SetTitle(Form(" ;BGO Energy (GeV); Width %s ", ParticleTag.c_str()));
		WidthHMC->SetFillStyle(0);  
		WidthHMC->Fit("ChargeFitFuncWidthHMC","R"); 		
		WidthHMC->Draw("AP");

		cout << "TString PWidthf = \"(" << ChargeFitFuncWidthH->GetParameter(0) << "+(" << ChargeFitFuncWidthH->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))" << "+(" << ChargeFitFuncWidthH->GetParameter(2) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)" << "+(" << ChargeFitFuncWidthH->GetParameter(3) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(" << ChargeFitFuncWidthH->GetParameter(4) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4)" << ")\";" << endl;
    	cout << "TString PMcWidthf = \"(" << ChargeFitFuncWidthHMC->GetParameter(0) << "+(" << ChargeFitFuncWidthHMC->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))" << "+(" << ChargeFitFuncWidthHMC->GetParameter(2) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**2)" << "+(" << ChargeFitFuncWidthHMC->GetParameter(3) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**3)+(" << ChargeFitFuncWidthHMC->GetParameter(4) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3)**4)" << ")\";" << endl;
    

   	TMultiGraph *mg2 = new TMultiGraph();
	mg2->SetTitle(Form(" ;BGO Energy (GeV);Width %s ", ParticleTag.c_str()));
    	mg2->Add(WidthH);
	mg2->Add(WidthHMC);
	mg2->Draw("ap");

mg2->GetXaxis()->SetLimits(20, 33000.);
DrawFitParameters(ChargeFitFuncWidthH,  0.15, 0.65, kBlue+1, "Data", 5);
DrawFitParameters(ChargeFitFuncWidthHMC, 0.4, 0.65, kRed+1,  "MC",   5);

TLegend* leg1 = new TLegend(0.1352705,0.6745306,0.3857715,0.8748195); //0.1352705,0.6745306,0.3857715,0.8748195,NULL,"brNDC"
leg1->SetFillColor(0);
leg1->SetBorderSize(0);
leg1->SetTextSize(0.04);
leg1->AddEntry(WidthHMC,"MC data","ep");
leg1->AddEntry(WidthH,"Flight data","ep");

leg1->Draw();

cout << "\n * * * * * * *    Sigma    * * * * * * * \n" << endl;

		TF1 *FitFuncGausH  = new TF1("FitFuncGausH","[0]",25.,100000.0);
//		TF1 *FitFuncGausH  = new TF1("FitFuncGausH","[0]+[1]*x",25.0,100000.0);
//		TF1 *FitFuncGausH  = new TF1("FitFuncGausH","[0]+[1]*x",70.0,39820.0);

		//FitFuncGausH->SetParameters(0.5,1.0e-50);
//		FitFuncGausH->FixParameter(0,1e-06);
//		FitFuncGausH->FixParameter(1,0.0);

		FitFuncGausH->SetLineColor(kBlue+1);
			
		TCanvas *C7 = new TCanvas("C7","C7",1000,800);
		
		TGraphErrors *GausSigmaH = new TGraphErrors(14, BGOEnergy, SigmaGH, E_BGOEnergy, E_SigmaGH);
		
		C7->SetLogx();
		C7->SetTickx();
		C7->SetTicky();
		
		GausSigmaH->SetMarkerStyle(20);
		GausSigmaH->SetMarkerSize(1.2);
		GausSigmaH->SetMarkerColor(kBlue+1);
		GausSigmaH->SetLineColor(kBlue+1);
		GausSigmaH->SetTitle(Form(" ;BGO Energy (GeV);#sigma_{gaus} %s ", ParticleTag.c_str()));
		GausSigmaH->SetFillStyle(0); 
		GausSigmaH->Fit("FitFuncGausH","R");GausSigmaH->Fit("FitFuncGausH","R");
		GausSigmaH->Draw("AP");
		TF1 *FitFuncGausHMC  = new TF1("FitFuncGausHMC","[0]",63.10,15000.0);
//		TF1 *FitFuncGausHMC  = new TF1("FitFuncGausHMC","[0]+[1]*x",25.0,100000.0);
//		TF1 *FitFuncGausHMC  = new TF1("FitFuncGausHMC","[0]+[1]*x",25.0,25120.0);
//		TF1 *FitFuncGausHMC  = new TF1("FitFuncGausHMC","[0]+[1]*x",25.0,2600.0);

		//FitFuncGausHMC->SetParameters(0.5,1.0e-50);
		FitFuncGausHMC->SetLineColor(kRed+1);
		
//		TCanvas *C8 = new TCanvas("C8","C8",1000,800);
		
		TGraphErrors *GausSigmaHMC = new TGraphErrors(14, BGOEnergy, SigmaGHMC, E_BGOEnergy, E_SigmaGHMC);
/*		
		C8->SetLogx();
		C8->SetTickx();
		C8->SetTicky();
*/		
		GausSigmaHMC->SetMarkerStyle(20);
		GausSigmaHMC->SetMarkerSize(1.2);
		GausSigmaHMC->SetMarkerColor(kRed+1);
		GausSigmaHMC->SetLineColor(kRed+1);
		GausSigmaHMC->SetTitle(Form(" ;BGO Energy (GeV);#sigma_{gaus} %s ", ParticleTag.c_str()));
		GausSigmaHMC->SetFillStyle(0);   
		GausSigmaHMC->Fit("FitFuncGausHMC","R");GausSigmaHMC->Fit("FitFuncGausHMC","R");
		GausSigmaHMC->Draw("AP");

//		cout << "TString PGSigmaf = \"(" << FitFuncGausH->GetParameter(0) << "+(" << FitFuncGausH->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))" << ")\";" << endl;
//    	cout << "TString PMcGSigmaf = \"(" << FitFuncGausHMC->GetParameter(0) << "+(" << FitFuncGausHMC->GetParameter(1) << "*log10(BGO_EnergyG_QuenchSatCorr_ML_ions_v3))" << ")\";" << endl;
cout << "TString PGSigmaf = \"(" << FitFuncGausH->GetParameter(0) << ")\";" << endl;
cout << "TString PMcGSigmaf = \"(" << FitFuncGausHMC->GetParameter(0) << ")\";" << endl;
      
  
   	TMultiGraph *mg3 = new TMultiGraph();
	mg3->SetTitle(Form(" ;BGO Energy (GeV);#sigma_{gaus} %s ", ParticleTag.c_str()));
    	mg3->Add(GausSigmaH);
	mg3->Add(GausSigmaHMC);
	mg3->Draw("ap");

mg3->GetXaxis()->SetLimits(20, 33000.);
DrawFitParameters(FitFuncGausH,   0.45, 0.35, kBlue+1, "Data", 1);
DrawFitParameters(FitFuncGausHMC, 0.7, 0.35, kRed+1,  "MC",   1);

TLegend* leg2 = new TLegend(0.1352705,0.6745306,0.3857715,0.8748195); //0.1352705,0.6745306,0.3857715,0.8748195,NULL,"brNDC"
leg2->SetFillColor(0);
leg2->SetBorderSize(0);
leg2->SetTextSize(0.04);
leg2->AddEntry(GausSigmaHMC,"MC data","ep");
leg2->AddEntry(GausSigmaH,"Flight data","ep");

leg2->Draw();

/*
C1->SaveAs(Form("PLOTS/MPVPlot_%s_%s.pdf", ParticleTag.c_str(), timestamp.c_str()));
C5->SaveAs(Form("PLOTS/WidthPlot_%s_%s.pdf", ParticleTag.c_str(), timestamp.c_str()));
C7->SaveAs(Form("PLOTS/SigmaPlot_%s_%s.pdf", ParticleTag.c_str(), timestamp.c_str()));
*/
}

