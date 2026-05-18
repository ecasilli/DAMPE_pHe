

void draw_scatterplot_DATA_PSD_BGO(){
    gStyle->SetPalette(55);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.04);
    gStyle->SetPadBottomMargin(0.09);

    TFile aa("scatterplot_PHe_PSD_BGO_Orb120Month.root");
    TH1F *hDATA = (TH1F*)aa.Get("h2psdbgo"); hDATA->SetStats(0);
    hDATA->GetXaxis()->SetTitleOffset(1.2);
    hDATA->GetYaxis()->SetTitle("PSD charge");
    hDATA->GetYaxis()->SetTitleOffset(1.1);
    hDATA->GetYaxis()->SetRangeUser(0.8,5.6);
    hDATA->GetXaxis()->SetRangeUser(20.,900000.);
    TCanvas *cc = new TCanvas("cc","cc",780,780); cc->SetLogx(); //cc->SetLogz();
    hDATA->Draw("colz"); 

    /////////////////////////////////////////////////////////////////////////////////
    
    // HELIUM 05/2026
	TString HeMPVf= "(1.81869+(0.290607*log10(x))+(-0.146116*log10(x)**2)+(0.0338788*log10(x)**3)+(-0.00192039*log10(x)**4) )";
	TString HeWidthf = "(-0.404366+(0.742843*log10(x))+(-0.422432*log10(x)**2)+(0.10251*log10(x)**3)+(-0.00846052*log10(x)**4))";
	TString HeGSigmaf = "(0.060102+(-1.43351e-39*log10(x)))";

	TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

	TF1 *mpv_he = new TF1("mpv_he", HeMPVf, 20.,900000.); mpv_he->SetLineColor(6); mpv_he->SetLineWidth(3); mpv_he->Draw("same");

	TF1 *ch_high= new TF1("ch_high",HeMPVf+"+"+"6*"+HeFSig, 20.,900000.); ch_high->SetLineColor(8); ch_high->SetLineStyle(3); ch_high->SetLineWidth(3); ch_high->Draw("same");


	// PROTON 05/2026
	TString PMPVf= "(1.01823+(-0.0471954*log10(x))+(0.0688881*log10(x)**2)+(-0.0306054*log10(x)**3)+(0.0057345*log10(x)**4) )";
	TString PWidthf = "(0.0115231+(0.0107847*log10(x))+(0.0168705*log10(x)**2)+(-0.00982277*log10(x)**3)+(0.00196218*log10(x)**4))";
	TString PGSigmaf = "(4.15872e-08+(-1.46051e-40*log10(x)))";
	
	TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";

	TF1 *mpv_p = new TF1("mpv_p", PMPVf, 20.,900000.); mpv_p->SetLineColor(9); mpv_p->SetLineWidth(3); mpv_p->Draw("same");

	TF1 *ch_low= new TF1("ch_low",PMPVf+"-"+"2*"+PFSig, 20.,900000.); ch_low->SetLineColor(8); ch_low->SetLineStyle(3); ch_low->SetLineWidth(3); ch_low->Draw("same");

	cc->SaveAs("scatterplot_PHe_PSD_BGO_Orb120Month.png");

}
