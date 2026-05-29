

void draw_scatterplot_DATA_PSD_BGO(){
    gStyle->SetPalette(55);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.04);
    gStyle->SetPadBottomMargin(0.09);

    TFile aa("ROOT_FILES/scatterplot_PHe_PSD_BGO_Orb120Month.root");
    TH1F *hDATA = (TH1F*)aa.Get("h2psdbgo"); hDATA->SetStats(0);
    hDATA->GetXaxis()->SetTitleOffset(1.2);
    hDATA->GetYaxis()->SetTitle("PSD charge");
    hDATA->GetYaxis()->SetTitleOffset(1.1);
    hDATA->GetYaxis()->SetRangeUser(0.8,5.4);
    hDATA->GetXaxis()->SetRangeUser(20.,900000.);
    TCanvas *cc = new TCanvas("cc","cc",1);//780,780); 
    cc->SetLogx(); cc->SetLogz();
    hDATA->Draw("colz"); 

    /////////////////////////////////////////////////////////////////////////////////
    
    TString HeMPVf= "(1.88586+(0.185264*log10(x))+(-0.0867115*log10(x)**2)+(0.0200701*log10(x)**3)+(-0.000810826*log10(x)**4) )";
    TString HeWidthf = "(-0.0824082+(0.201315*log10(x))+(-0.0989824*log10(x)**2)+(0.0218456*log10(x)**3)+(-0.00136192*log10(x)**4))";
    TString HeGSigmaf = "(0.0396516)";

	TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

	TF1 *mpv_he = new TF1("mpv_he", HeMPVf, 20.,900000.); mpv_he->SetLineColor(6); mpv_he->SetLineWidth(3); mpv_he->Draw("same");

	TF1 *ch_high= new TF1("ch_high",HeMPVf+"+"+"6*"+HeFSig, 20.,900000.); ch_high->SetLineColor(1); ch_high->SetLineStyle(2); ch_high->SetLineWidth(3); ch_high->Draw("same");

    TString PMPVf= "(0.993603+(-0.00402039*log10(x))+(0.0408375*log10(x)**2)+(-0.0227868*log10(x)**3)+(0.00494042*log10(x)**4) )";
    TString PWidthf = "(-0.0611545+(0.140328*log10(x))+(-0.0683705*log10(x)**2)+(0.0143596*log10(x)**3)+(-0.000527989*log10(x)**4))";
    TString PGSigmaf = "(5.5475e-08)";
	
	TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";

	TF1 *mpv_p = new TF1("mpv_p", PMPVf, 20.,900000.); mpv_p->SetLineColor(9); mpv_p->SetLineWidth(3); mpv_p->Draw("same");

	TF1 *ch_low= new TF1("ch_low",PMPVf+"-"+"2.9*"+PFSig, 20.,900000.); ch_low->SetLineColor(1); ch_low->SetLineStyle(2); ch_low->SetLineWidth(3); ch_low->Draw("same");

	cc->SaveAs("PLOTS/scatterplot_PHe_PSD_BGO_Orb120Month_2e5sigmaLow_6sigmaUp_new3.pdf");

}
