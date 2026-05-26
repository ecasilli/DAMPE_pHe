
void plot_acceptances() {

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
    const char *base = "unfold_results_pHe_2026_MLionsv3_2e5sigmaLow_6sigmaUp_new";

    TH1D *hcut00    = LoadHist(Form("%s_cut00.root",  base), "hcut00");
    TH1D *hcut01    = LoadHist(Form("%s_cut01.root",  base), "hcut01");
    TH1D *hcut02    = LoadHist(Form("%s_cut02.root",  base), "hcut02");
    TH1D *hcut05    = LoadHist(Form("%s_cut05.root",  base), "hcut05");
    TH1D *hcut06    = LoadHist(Form("%s_cut06.root",  base), "hcut06");
    TH1D *hcutSpCut = LoadHist(Form("%s_SpCut.root",  base), "hcutSpCut");
    TH1D *hcutChSel = LoadHist(Form("%s.root", base), "hcutChSel"); // nuovo

    // === Stile comune ===
    auto SetStyle = [](TH1D *h, int color, int lw = 2) {
        if (!h) return;
        h->SetLineColor(color);
        h->SetLineWidth(lw);
        h->SetStats(0);
    };

    SetStyle(hcut00,    kCyan+1);
    //SetStyle(hcut01,    col2);
    SetStyle(hcut02,    46);
    SetStyle(hcut05,    col2);
    SetStyle(hcut06,    8);
    SetStyle(hcutSpCut, col4);
    SetStyle(hcutChSel, 9);

    // === Titoli ===
    if (hcut00) {
        hcut00->SetTitle(" ");
        hcut00->GetXaxis()->SetTitle("MC true energy (GeV)");
        hcut00->GetYaxis()->SetTitle("Acceptance (m^{2} sr)");
        hcut00->GetYaxis()->SetTitleOffset(1.4);
        hcut00->GetXaxis()->SetTitleSize(0.045);
        hcut00->GetYaxis()->SetTitleSize(0.045);
        hcut00->GetXaxis()->SetLabelSize(0.040);
        hcut00->GetYaxis()->SetLabelSize(0.040);
        hcut00->GetYaxis()->SetRangeUser(0.,0.30);
    }

    // === Disegno ===
    if (hcut00)    hcut00->Draw("HIST");
    //if (hcut01)    hcut01->Draw("HIST same");
    if (hcut02)    hcut02->Draw("HIST same");
    if (hcut05)    hcut05->Draw("HIST same");
    if (hcut06)    hcut06->Draw("HIST same");
    if (hcutSpCut) hcutSpCut->Draw("HIST same");
    if (hcutChSel) hcutChSel->Draw("HIST same");

    // === Legenda ===
    TLegend *legend = new TLegend(0.18, 0.55, 0.58, 0.85);
    legend->SetTextSize(0.028);
    legend->SetLineWidth(1);
    legend->SetBorderSize(1);
    legend->SetFillStyle(1001);
    legend->SetFillColor(0);

    if (hcut00)    legend->AddEntry(hcut00,    "HET selection",                                        "l");
    //if (hcut01)    legend->AddEntry(hcut01,    "PSD charge on Y and X",                                "l");
    if (hcut02)    legend->AddEntry(hcut02,    "E_{BGO,0}+E_{BGO,1} < E_{BGO,2}+E_{BGO,3}",           "l");
    if (hcut05)    legend->AddEntry(hcut05,    "BGO-STK match",                                        "l");
    if (hcut06)    legend->AddEntry(hcut06,    "STK track in PSD volume",                              "l");
    if (hcutSpCut) legend->AddEntry(hcutSpCut, "Electron cut",                                        "l");
    if (hcutChSel) legend->AddEntry(hcutChSel, "Charge selection",                                    "l");

    legend->Draw();

    c1->RedrawAxis();

    // === Salvataggio ===
    c1->SaveAs("PLOTS/partial_acceptances_perCut_20260526_2.pdf");
    c1->SaveAs("PLOTS/partial_acceptances_perCut_20260526_2.png");

    ::Info("plot_acceptances", "Plot salvato con successo.");
}

