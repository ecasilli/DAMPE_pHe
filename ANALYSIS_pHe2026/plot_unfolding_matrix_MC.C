
void plot_unfolding_matrix_MC() {

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
        c->SetLogz();
        c->SetLeftMargin(0.1);
        c->SetRightMargin(0.12);
        c->SetBottomMargin(0.1);
        c->SetTopMargin(0.07);
        return c;
    };


   // === Canvas ===
    TCanvas *c1 = MakeCanvas("c1", "Response matrix");
    TCanvas *c2 = MakeCanvas("c2", "Unfolded matrix");

    // === Caricamento istogrammi ===
    const char *fAll = "ROOT_FILES/unfold_results_pHe_2026_Orb120Month_except25low_MLionsv3_2e5sigmaLow_6sigmaUp_new.root";

    TH1D *hresponse = LoadHist(fAll, "h2Ntrig_wgt_v3", "hresponse");
    TH1D *hunfolded = LoadHist(fAll, "unfold_matrix",  "hunfolded");

    // === Titoli ===
    SetTitle(hresponse, "MC true energy (GeV)", "BGO deposited energy (GeV)");
    SetTitle(hunfolded, "MC true energy (GeV)", "BGO deposited energy (GeV)");

    // === Disegno — canvas 1 (counts before) ===
    c1->cd();
    if (hresponse) hresponse->Draw("COLZ");

    // === Disegno — canvas 2 (counts after) ===
    c2->cd();
    if (hunfolded) hunfolded->Draw("COLZ");
    

    // === Salvataggio in *.png e *.pdf ==
    c1->SaveAs("PLOTS/response_matrix_pHe_2026.pdf");
    c1->SaveAs("PLOTS/response_matrix_pHe_2026.png");

    c2->SaveAs("PLOTS/unfolded_matrix_pHe_2026.pdf");
    c2->SaveAs("PLOTS/unfolded_matrix_pHe_2026.png");


}