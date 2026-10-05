
void plot_acceptances_EPOSLHC_comp() {

    const int col1 = TColor::GetColor("#3190da");  // kP10Blue 
    const int col2 = TColor::GetColor("#ffa90e");  // kP10Yellow
    const int col3 = TColor::GetColor("#bd1f01");  // kP10Red
    const int col4 = TColor::GetColor("#94a4a2");  // kP10Gray
    const int col5 = TColor::GetColor("#832db6");  // kP10Violet
    const int col6 = TColor::GetColor("#a96b59");  // kP10Brown
    const int col7 = TColor::GetColor("#92dadd");  // kP10Cyan

    // === Stile globale ===
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    // === Canvas ===
    TCanvas *c1 = new TCanvas("c1", "Partial Acceptances", 900, 650);
    c1->SetLogx();
    c1->SetGridx();
    c1->SetGridy();
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    // ---------------------------------------------------------------
    // Funzione lambda per caricare un clone di hacc1 da un file ROOT
    // ---------------------------------------------------------------
    auto LoadHist = [](const char* filename, const char* newname) -> TH1D* {
        TFile *f = TFile::Open(filename, "READ");
        if (!f || f->IsZombie()) {
            ::Error("LoadHist", "Impossibile aprire %s", filename);
            return nullptr;
        }
        TH1D *h = (TH1D*) f->Get("hacc1");
        if (!h) {
            ::Error("LoadHist", "hacc1 non trovata in %s", filename);
            return nullptr;
        }
        TH1D *hclone = (TH1D*) h->Clone(newname);
        hclone->SetDirectory(0);  // sopravvive alla chiusura del file
        f->Close();
        return hclone;
    };

    // === Caricamento istogrammi ===
    const char *base = "ROOT_FILES/unfold_results_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_24sett26_wPHe_kernel_5bin";
    const char *base1 = "ROOT_FILES/unfold_results_pHe_2026_Orb120Month_3sigmaLow_6sigmaUp_PSDprogr_STKcharge450_comb_vert0e7_10TeV_EPOSLHC";

    TH1D *hcutChSel  = LoadHist(Form("%s.root", base ), "hcutChSel"); 
    TH1D *hcutChSel1 = LoadHist(Form("%s.root", base1), "hcutChSel"); 

    // === Stile comune ===
    auto SetStyle = [](TH1D *h, int color, int lw = 2) {
        if (!h) return;
        h->SetLineColor(color);
        h->SetLineWidth(lw);
        h->SetStats(0);
    };

    SetStyle(hcutChSel , 2);
    SetStyle(hcutChSel1, 9);

    // === Titoli ===
    if (hcutChSel) {
        hcutChSel->SetTitle(" ");
        hcutChSel->GetXaxis()->SetTitle("MC true energy (GeV)");
        hcutChSel->GetYaxis()->SetTitle("Acceptance (m^{2} sr)");
        hcutChSel->GetYaxis()->SetTitleOffset(1.4);
        hcutChSel->GetXaxis()->SetTitleSize(0.045);
        hcutChSel->GetYaxis()->SetTitleSize(0.045);
        hcutChSel->GetXaxis()->SetLabelSize(0.040);
        hcutChSel->GetYaxis()->SetLabelSize(0.040);
        hcutChSel->GetYaxis()->SetRangeUser(0.,0.145);
    }

    // === Disegno ===
    if (hcutChSel ) hcutChSel->Draw("HIST ");
    if (hcutChSel1) hcutChSel1->Draw("HIST same");

    // === Legenda ===
    TLegend *legend = new TLegend(0.2, 0.72, 0.55, 0.85);
    legend->SetTextSize(0.03);
    legend->SetLineColor(0);
    //legend->SetLineWidth(1);
    //legend->SetBorderSize(1);
    legend->SetFillStyle(1001);
    legend->SetFillColor(0);

    if (hcutChSel ) legend->AddEntry(hcutChSel, "10-100TeV FTFP",   "l");
    if (hcutChSel1) legend->AddEntry(hcutChSel1, "10-100TeV EPOSLHC_FTFP",   "l");

    legend->Draw();

    c1->RedrawAxis();

    // === Salvataggio ===
    c1->SaveAs("PLOTS/partial_acceptances_smooth_3sLow_6sUP_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC_final.pdf");
    c1->SaveAs("PLOTS/partial_acceptances_smooth_3sLow_6sUP_PSDprogr_STKch450_comb_vert0e7_10TeV_EPOSLHC_final.png");

    ::Info("plot_acceptances", "Plot salvato con successo.");
}

