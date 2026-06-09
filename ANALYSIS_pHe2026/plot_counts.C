
void plot_counts() {

    // === Stile globale ===
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // ---------------------------------------------------------------
    // Funzione lambda per caricare un clone TH1D da un file ROOT
    // ---------------------------------------------------------------
    auto LoadHist = [](const char* filename, const char* histname, const char* newname) -> TH1D* {
        TFile *f = TFile::Open(filename, "READ");
        if (!f || f->IsZombie()) {
            ::Error("LoadHist", "Impossibile aprire %s", filename);
            return nullptr;
        }
        TH1D *h = (TH1D*) f->Get(histname);
        if (!h) {
            ::Error("LoadHist", "%s non trovata in %s", histname, filename);
            return nullptr;
        }
        TH1D *hclone = (TH1D*) h->Clone(newname);
        hclone->SetDirectory(0);  // sopravvive alla chiusura del file
        f->Close();
        return hclone;
    };

    // ---------------------------------------------------------------
    // Funzione lambda per clonare un istogramma e dividerlo per il livetime
    // restituisce un nuovo TH1D con contenuto = counts/livetime
    // ---------------------------------------------------------------
    auto MakeRate = [](TH1D* h, double livetime, const char* newname) -> TH1D* {
        if (!h || livetime <= 0.) {
            ::Error("MakeRate", "Istogramma nullo o livetime non valido");
            return nullptr;
        }
        TH1D *hrate = (TH1D*) h->Clone(newname);
        hrate->SetDirectory(0);
        hrate->Scale(1.0 / livetime);   // bin content e errori scalati coerentemente
        return hrate;
    };

    // === Stile comune ===
    auto SetStyle = [](TH1D *h, int color, int mstyle, double msize ) { //lw = 2) {
        if (!h) return;
        h->SetLineColor(color);
        h->SetMarkerColor(color);
        h->SetMarkerStyle(mstyle);
        h->SetMarkerSize(msize);
        h->SetStats(0);
    };

    // === Titolo comune ===
    auto SetTitle = [](TH1D *h, const char* titleX, const char* titleY ) { 
        if (!h) return;
        h->SetTitle(" ");
        h->GetXaxis()->SetTitle(titleX);
        h->GetYaxis()->SetTitle(titleY);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetTitleOffset(1.4);
        h->GetXaxis()->SetTitleSize(0.03);
        h->GetYaxis()->SetTitleSize(0.03);
        h->GetXaxis()->SetLabelSize(0.03);
        h->GetYaxis()->SetLabelSize(0.03);
        //h->GetYaxis()->SetRangeUser(0.,0.30);
    };

    // ---------------------------------------------------------------
    // Funzione lambda per creare e configurare una TCanvas standard
    // ---------------------------------------------------------------
    auto MakeCanvas = [](const char* name, const char* title) -> TCanvas* {
        TCanvas *c = new TCanvas(name, title, 900, 650);
        c->SetLogx();
        c->SetLogy();
        c->SetGridx();
        c->SetGridy();
        c->SetLeftMargin(0.1);
        c->SetRightMargin(0.06);
        c->SetBottomMargin(0.1);
        c->SetTopMargin(0.07);
        return c;
    };

    // ---------------------------------------------------------------
    // === Livetime nei 3 casi (in secondi) ===
    // ---------------------------------------------------------------
    const double lt_All = 242576599.4; // ~10 anni in secondi 
    const double lt_6ys = 2183*86400*0.7667; // ~ 6 anni in secondi 
    const double lt_4ys = 242576599-(2183*86400*0.7667);   // ~ 4 anni in secondi 
    //const double lt_9ys = 217488673.3;  // ~9 years


   // === Canvas ===
    TCanvas *c1 = MakeCanvas("c1", "Counts before unfolding");
    TCanvas *c2 = MakeCanvas("c2", "Counts after unfolding");
    TCanvas *c3 = MakeCanvas("c3", "Rate before unfolding");
    TCanvas *c4 = MakeCanvas("c4", "Rate after unfolding");

    // === Caricamento istogrammi ===
    const char *fAll = "ROOT_FILES/unfold_results_pHe_2026_Orb120Month_except25low_MLionsv3_2e5sigmaLow_6sigmaUp_new.root";
    //const char *fAll = "ROOT_FILES/unfold_results_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.root";
    const char *f6ys = "ROOT_FILES/unfold_results_pHe_2026_72months_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.root";
    const char *f4ys = "ROOT_FILES/unfold_results_pHe_2026_rem48months_MLionsv3_2e5sigmaLow_6sigmaUp_new_smooth.root";

    TH1D *hAllbef = LoadHist(fAll, "h1SelBGO_orb_v3", "hAllbef");
    TH1D *hAllaft = LoadHist(fAll, "unfold_cts",      "hAllaft");
    TH1D *h6ysbef = LoadHist(f6ys, "h1SelBGO_orb_v3", "h6ysbef");
    TH1D *h6ysaft = LoadHist(f6ys, "unfold_cts",      "h6ysaft");
    TH1D *h4ysbef = LoadHist(f4ys, "h1SelBGO_orb_v3", "h4ysbef");
    TH1D *h4ysaft = LoadHist(f4ys, "unfold_cts",      "h4ysaft");
    
    SetStyle(hAllbef, kBlack,    20, 0.7);
    SetStyle(hAllaft, kBlack,    20, 0.7);
    SetStyle(h6ysbef, kOrange+1, 20, 0.7);
    SetStyle(h6ysaft, kOrange+1, 20, 0.7);
    SetStyle(h4ysbef, kGreen+1,  20, 0.7);
    SetStyle(h4ysaft, kGreen+1,  20, 0.7);

    // === Istogrammi di rate (counts / livetime) ===
    TH1D *hAllbef_rate = MakeRate(hAllbef, lt_All, "hAllbef_rate");
    TH1D *hAllaft_rate = MakeRate(hAllaft, lt_All, "hAllaft_rate");
    TH1D *h6ysbef_rate = MakeRate(h6ysbef, lt_6ys, "h6ysbef_rate");
    TH1D *h6ysaft_rate = MakeRate(h6ysaft, lt_6ys, "h6ysaft_rate");
    TH1D *h4ysbef_rate = MakeRate(h4ysbef, lt_4ys, "h4ysbef_rate");
    TH1D *h4ysaft_rate = MakeRate(h4ysaft, lt_4ys, "h4ysaft_rate");

    SetStyle(hAllbef_rate, kBlack,    20, 0.7);
    SetStyle(hAllaft_rate, kBlack,    20, 0.7);
    SetStyle(h6ysbef_rate, kOrange+1, 20, 0.7);
    SetStyle(h6ysaft_rate, kOrange+1, 20, 0.7);
    SetStyle(h4ysbef_rate, kGreen+1,  20, 0.7);
    SetStyle(h4ysaft_rate, kGreen+1,  20, 0.7);

    // === Titoli ===
    SetTitle(hAllbef, "BGO deposited energy (GeV)", "Counts");
    SetTitle(hAllaft, "Primary energy (GeV)", "Counts after unfolding");

    SetTitle(hAllbef_rate, "BGO deposited energy (GeV)", "Rate (s^{-1})");
    SetTitle(hAllaft_rate, "Primary energy (GeV)", "Rate after unfolding (s^{-1})");
    hAllbef_rate->GetYaxis()->SetRangeUser(1.5e-9, 0.7);
    hAllaft_rate->GetYaxis()->SetRangeUser(7e-11, 0.7);
/*
    // === Legenda ===
    TLegend *legend = new TLegend(0.64, 0.63, 0.88, 0.86);
    legend->SetTextSize(0.028);
    legend->SetLineWidth(1);
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);
    legend->SetFillColor(0);

    if (hAllbef) legend->AddEntry(hAllbef, "10 years (2016-2025)", "pe");
    if (h6ysbef) legend->AddEntry(h6ysbef, "6 years (2016-2021)", "pe");
    if (h4ysbef) legend->AddEntry(h4ysbef, "4 years (2022-2025)", "pe");
*/
    // === Disegno — canvas 1 (counts before) ===
    c1->cd();
    if (hAllbef) hAllbef->Draw("E");
    //if (h6ysbef) h6ysbef->Draw("E same");
    //if (h4ysbef) h4ysbef->Draw("E same");
    //legend->Draw();

    // === Disegno — canvas 2 (counts after) ===
    c2->cd();
    if (hAllaft) hAllaft->Draw("E");
    //if (h6ysaft) h6ysaft->Draw("E same");
    //if (h4ysaft) h4ysaft->Draw("E same");
    //legend->Draw();

    // === Disegno — canvas 3 (rate before) ===
    c3->cd();
    if (hAllbef_rate) hAllbef_rate->Draw("E");
    //if (h6ysbef_rate) h6ysbef_rate->Draw("E same");
    //if (h4ysbef_rate) h4ysbef_rate->Draw("E same");
    //legend->Draw();

    // === Disegno — canvas 4 (rate after) ===
    c4->cd();
    if (hAllaft_rate) hAllaft_rate->Draw("E");
    //if (h6ysaft_rate) h6ysaft_rate->Draw("E same");
    //if (h4ysaft_rate) h4ysaft_rate->Draw("E same");
    //legend->Draw();

    // === Salvataggio in *.png e *.pdf ==
    c1->SaveAs("PLOTS/counts_before_unfolding_pHe_Orb120Month_except25low.pdf");
    c1->SaveAs("PLOTS/counts_before_unfolding_pHe_Orb120Month_except25low.png");

    c2->SaveAs("PLOTS/counts_after_unfolding_pHe_Orb120Month_except25low.pdf");
    c2->SaveAs("PLOTS/counts_after_unfolding_pHe_Orb120Month_except25low.png");

    c3->SaveAs("PLOTS/rate_before_unfolding_pHe_Orb120Month_except25low.pdf");
    c3->SaveAs("PLOTS/rate_before_unfolding_pHe_Orb120Month_except25low.png");

    c4->SaveAs("PLOTS/rate_after_unfolding_pHe_Orb120Month_except25low.pdf");
    c4->SaveAs("PLOTS/rate_after_unfolding_pHe_Orb120Month_except25low.png");


}