#include <vector>
#include <iostream>

void plot_STK_charge_per_BGO_energy(TString mode = "DATA") {

    // ---- Switch DATA / MC ----
    bool isMC;
    if (mode == "DATA") {
        isMC = false;
    } else if (mode == "MC") {
        isMC = true;
    } else {
        std::cerr << "Errore: argomento '" << mode
                   << "' non valido. Usa \"DATA\" oppure \"MC\"." << std::endl;
        return;
    }

    std::cout << "Analyzing " << mode << std::endl; 

    // ---- Parametri generali ----
    const char* fileName = isMC
        ? "PHe_MC_p_He_5PeV_STK_charge_200bins.root"
        : "PHe_STK_Orb120Month_240bins.root";
    const double xLineADC  = 400.0;   // STK charge cut (ADC counts)

    int histColor   = isMC ? kRed+1  : kBlue+1;   // colore istogramma
    int lineColor   = isMC ? kBlue+1 : kRed+1;    // colore linea verticale
    const char* tag = isMC ? "MC" : "DATA";       // prefisso per il salvataggio

    gStyle->SetOptStat(0); // niente box delle statistiche
    gStyle->SetOptTitle(1);

    TFile* f = TFile::Open(fileName, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Errore: impossibile aprire il file " << fileName << std::endl;
        return;
    }

    const int nHist = 5;
    const char* hNames[nHist]  = {"h01", "h02", "h03", "h04", "h05"};
    const char* hTitles[nHist] = {
        "10 GeV < E_{BGO} < 100 GeV",
        "100 GeV < E_{BGO} < 1 TeV",
        "1 TeV < E_{BGO} < 10 TeV",
        "10 TeV < E_{BGO} < 100 TeV",
        "100 TeV < E_{BGO} < 1 PeV"
    };

    std::vector<TCanvas*> canvases;

    for (int i = 0; i < nHist; ++i) {
 
        TH1F* h = (TH1F*)f->Get(hNames[i]);
        if (!h) {
            std::cerr << "Warning: " << hNames[i] << " doesn't exist, skip." << std::endl;
            continue;
        }
 
        TString cname = Form("c_%s", hNames[i]);
        TCanvas* c = new TCanvas(cname, hTitles[i], 50 + i*30, 50 + i*30, 800, 600);
        //c->SetGrid();
        canvases.push_back(c);
 
        // ---- Stile dell'istogramma ----
        h->SetTitle("");
        h->SetLineColor(histColor);
        h->SetMarkerColor(histColor);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.6);
        h->GetXaxis()->SetTitle("STK signal (ADC counts)");
        h->GetYaxis()->SetTitle("Number of events");
        h->GetXaxis()->SetTitleSize(0.035);
        h->GetYaxis()->SetTitleSize(0.035);
        h->GetXaxis()->SetLabelSize(0.035);
        h->GetYaxis()->SetLabelSize(0.035);
        h->GetXaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetTitleOffset(1.2);
        h->GetYaxis()->SetMaxDigits(3);
 
        h->Draw("PE"); 
        c->Update(); 
 
        // ---- Linea verticale a xLineADC ----
        double yMin = gPad->GetUymin();
        double yMax = gPad->GetUymax();
 
        TLine* line = new TLine(xLineADC, yMin, xLineADC, yMax);
        line->SetLineColor(lineColor);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");

        // ---- Titolo dentro l'area del plot ----
        TPaveText* pt = new TPaveText(0.5, 0.8, 0.85, 0.87, "NDC");
        pt->SetBorderSize(0);
        pt->SetFillStyle(0);      // trasparente, niente riquadro
        pt->SetTextFont(42);
        pt->SetTextSize(0.04);
        pt->SetTextAlign(12);     // allineato a sinistra
        pt->AddText(hTitles[i]);
        pt->Draw("same");
 
        // ---- Legenda: solo la linea, il titolo è già nel riquadro sopra il plot ----
        TLegend* leg = new TLegend(0.58, 0.78, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.035);
        leg->AddEntry(line, Form("ADC = %.0f", xLineADC), "l");
        //leg->Draw();
 
        c->Modified();
        c->Update();
 
        // ---- Salvataggio su file ----
        c->SaveAs(Form("PLOTS/%s_%s.png", tag, hNames[i]));
        c->SaveAs(Form("PLOTS/%s_%s.pdf", tag, hNames[i]));
    }




}