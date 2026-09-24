

void draw_scatterplot_DATA_PSD_BGO(){
    gStyle->SetPalette(55);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadLeftMargin(0.09);
    gStyle->SetPadRightMargin(0.11);
    gStyle->SetPadTopMargin(0.04);
    gStyle->SetPadBottomMargin(0.09);

    TFile aa("ROOT_FILES/scatterplot_PHe_PSD_BGO_Orb120Month_23sett26.root");
    TH1F *hDATA = (TH1F*)aa.Get("h2psdbgo"); hDATA->SetStats(0);
    hDATA->GetXaxis()->SetTitleOffset(1.2);
    hDATA->GetYaxis()->SetTitle("PSD charge");
    hDATA->GetYaxis()->SetTitleOffset(1.1);
    hDATA->GetYaxis()->SetRangeUser(0.3,5.4);
    hDATA->GetXaxis()->SetRangeUser(20.,900000.);
    TCanvas *cc = new TCanvas("cc","cc",1);//780,780); 
    cc->SetLogx(); cc->SetLogz();
    hDATA->Draw("colz"); 

    /////////////////////////////////////////////////////////////////////////////////
    // * * * * * * * HELIUM
    
    TString HeMPVf= "(1.92403+(0.1541*log10(x))+(-0.082833*log10(x)**2)+(0.0208528*log10(x)**3)+(-0.0011579*log10(x)**4) )";
    TString HeWidthf = "(-0.0105333+(0.129859*log10(x))+(-0.0713011*log10(x)**2)+(0.0169532*log10(x)**3)+(-0.00117928*log10(x)**4))";
    TString HeGSigmaf = "(0.0227879)";


	TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";

	TF1 *mpv_he = new TF1("mpv_he", HeMPVf, 20.,900000.); mpv_he->SetLineColor(6); mpv_he->SetLineWidth(3); mpv_he->Draw("same");

	TF1 *ch_high_he= new TF1("ch_high_he",HeMPVf+"+"+"6*"+HeFSig, 20.,900000.); ch_high_he->SetLineColor(1); ch_high_he->SetLineStyle(2); ch_high_he->SetLineWidth(3); ch_high_he->Draw("same");
    TF1 *ch_low_he = new TF1("ch_low_he", HeMPVf+"-"+"2*"+HeFSig, 20.,900000.); ch_low_he->SetLineColor(1); ch_low_he->SetLineStyle(2); ch_low_he->SetLineWidth(3); ch_high_he->Draw("same");

    //////////////////////////////////////////////////////////////////////////////////
    // * * * * * * * PROTON

    TString PMPVf= "(0.988625+(0.0166979*log10(x))+(0.0136979*log10(x)**2)+(-0.0114983*log10(x)**3)+(0.00313886*log10(x)**4) )";
    TString PWidthf = "(-0.0152013+(0.0753179*log10(x))+(-0.0344289*log10(x)**2)+(0.00656681*log10(x)**3)+(0.000105438*log10(x)**4))";
    TString PGSigmaf = "(2.02488e-08)";
	
	TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";

	TF1 *mpv_p = new TF1("mpv_p", PMPVf, 20.,900000.); mpv_p->SetLineColor(9); mpv_p->SetLineWidth(3); mpv_p->Draw("same");

	TF1 *ch_high_p= new TF1("ch_high_p",PMPVf+"-"+"2.6*"+PFSig, 20.,900000.); ch_high_p->SetLineColor(1); ch_high_p->SetLineStyle(2); ch_high_p->SetLineWidth(3); ch_high_p->Draw("same");
    TF1 *ch_low_p = new TF1("ch_low_p", PMPVf+"+"+"4.0*"+PFSig, 20.,900000.); ch_low_p->SetLineColor(1); ch_low_p->SetLineStyle(2); ch_low_p->SetLineWidth(3); ch_low_p->Draw("same");

	//========================
    // Legenda in alto a sinistra
    //========================
    TLegend *leg = new TLegend(0.12, 0.55, 0.39, 0.9);
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);

    leg->AddEntry(ch_high_he, "Helium MPV + 6 #sigma", "l");
    leg->AddEntry(mpv_he,  "Helium MPV",          "l");
    leg->AddEntry(ch_low_he, "Helium MPV - 2 #sigma", "l");

    leg->AddEntry(ch_high_p,  "Proton MPV + 4.0 #sigma", "l");
    leg->AddEntry(mpv_p,   "Proton MPV",          "l");
    leg->AddEntry(ch_low_p,  "Proton MPV - 2.6 #sigma", "l");

    leg->Draw();

    cc->SaveAs("PLOTS/scatterplot_PHe_PSDprogressive_BGO_Orb120Month_23sett26.pdf");

}
