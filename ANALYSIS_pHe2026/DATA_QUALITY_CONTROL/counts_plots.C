/*
 * plot_counts_dampe.C
 * ===================
 * Legge i file DAMPE pHe organizzati per anno e range di energia,
 * e produce:
 *
 *   1) plot_total.pdf  — tutti gli anni, bin = 1 settimana, 6 colori per energia
 *   2) plot_YYYY.pdf   — una canvas per anno (2016-2025), bin = 1 giorno
 *
 * NOTA: il branch "time" è Int_t e conta secondi dall'epoch DAMPE.
 *       Imposta TIME_OFFSET qui sotto (secondi Unix dell'epoch DAMPE).
 *
 * Utilizzo:
 *   nohup root -l -b -q plot_counts_dampe.C > plot_counts_dampe.log 2>&1 &
 */

#include <iostream>

// ── USER CONFIG ───────────────────────────────────────────────────────────────
const int START_YEAR = 2016;
const int END_YEAR   = 2025;
const int WEEK       = 7 * 86400;
const int DAY        = 86400;

// Offset dell'epoch DAMPE in secondi Unix.
// "time" nel TTree = secondi dal 01/01/2013 00:00 UTC (epoch interna DAMPE).
// Calcolato come: unix(01/01/2016 00:00 UTC) - 94608350 = 1356998050
// dove 94608350 è il primo valore di "time" nel file del 01/01/2016.
// => time=0 corrisponde a ~01/01/2013 00:00 UTC
const Int_t TIME_OFFSET = 1356998050;

struct EnergyBin {
    const char* label;
    const char* suffix;
    Color_t     color;
};

// Matplotlib tab10 palette
const int mplBlue   = 2001;
const int mplOrange = 2002;
const int mplGreen  = 2003;
const int mplRed    = 2004;
const int mplPurple = 2005;
const int mplBrown  = 2006;

const int NEBINS = 6;
EnergyBin ebins[NEBINS] = {
    { "20-100 GeV",  "002_010", mplBlue   },
    { "100-250 GeV", "010_025", mplOrange },
    { "250-500 GeV", "025_050", mplGreen  },
    { "0.5-1 TeV",   "050_100", mplRed    },
    { "1-5 TeV",     "100_500", mplPurple },
    { ">5 TeV",      "500_000", mplBrown  },
};
// ─────────────────────────────────────────────────────────────────────────────

TString getBasePath() {
    TString host = gSystem->HostName();
    if (host.Contains("cnaf"))     return "/storage/gpfs_data/dampe/users";
    if (host.Contains("le.infn"))  return "/nfs/argo/dampe";
    std::cout << "WARNING: hostname non riconosciuto: " << host << std::endl;
    return ".";
}

void addFile(TChain* ch, const TString& base, int year, const char* suffix) {
    TString fname;
    fname = base + TString::Format("/SKIM_2026_pHe/FLIGHT/skim_flight_%s_%d_merged.root", suffix, year);
    ch->Add(fname);
}

void setTimeAxis(TH1D* h, bool shortFmt = false) {
    h->GetXaxis()->SetTimeDisplay(1);
    h->GetXaxis()->SetTimeOffset(0, "gmt");
    h->GetXaxis()->SetTimeFormat(shortFmt ? "%m/%Y" : "%d/%m/%Y");
    h->GetXaxis()->SetNdivisions(506);
    h->GetXaxis()->SetLabelSize(0.038);
    h->GetXaxis()->SetLabelOffset(0.01);
}

TLegend* makeLegend(TH1D* hh[], int n) {
    TLegend* leg = new TLegend(0.88, 0.50, 0.98, 0.88);
    leg->SetTextSize(0.032);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    for (int i = 0; i < n; i++)
        leg->AddEntry(hh[i], ebins[i].label, "l");
    return leg;
}

inline Double_t toUnix(Int_t t_dampe) {
    return (Double_t)(t_dampe) + (Double_t)TIME_OFFSET;
}

Double_t yearStartUnix(int year) {
    // TDatime::Convert(kTRUE) restituisce secondi Unix in UTC
    TDatime d(year, 1, 1, 0, 0, 0);
    return (Double_t)d.Convert(kTRUE);
}

// ── MAIN ──────────────────────────────────────────────────────────────────────
void counts_plots() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetPadGridY(1);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    new TColor(mplBlue,   31./255.,119./255.,180./255.);
    new TColor(mplOrange,255./255.,127./255., 14./255.);
    new TColor(mplGreen,  44./255.,160./255., 44./255.);
    new TColor(mplRed,   214./255., 39./255., 40./255.);
    new TColor(mplPurple,148./255.,103./255.,189./255.);
    new TColor(mplBrown, 140./255., 86./255., 75./255.);

    TString base = getBasePath();
    std::cout << "BasePath: " << base << std::endl;
    std::cout << "TIME_OFFSET: " << TIME_OFFSET << "\n" << std::endl;

    Double_t t_global_min = yearStartUnix(START_YEAR);
    Double_t t_global_max = yearStartUnix(END_YEAR + 1);
    int n_weeks = (int)((t_global_max - t_global_min) / WEEK) + 1;

    TH1D* h_tot[NEBINS];
    for (int ie = 0; ie < NEBINS; ie++) {
        h_tot[ie] = new TH1D(
            Form("h_tot_%d", ie), "",
            n_weeks, t_global_min, t_global_min + n_weeks * WEEK);
        h_tot[ie]->SetLineColor(ebins[ie].color);
        h_tot[ie]->SetLineWidth(2);
    }

    for (int yr = START_YEAR; yr <= END_YEAR; yr++) {
        std::cout << "\n=== Anno " << yr << " ===" << std::endl;

        Double_t t_yr_min = yearStartUnix(yr);
        Double_t t_yr_max = yearStartUnix(yr + 1);
        int n_days = (int)((t_yr_max - t_yr_min) / DAY);

        TH1D* h_yr[NEBINS];
        for (int ie = 0; ie < NEBINS; ie++) {
            h_yr[ie] = new TH1D(
                Form("h_yr_%d_%d", yr, ie), "",
                n_days, t_yr_min, t_yr_max);
            h_yr[ie]->SetLineColor(ebins[ie].color);
            h_yr[ie]->SetLineWidth(2);
        }

        for (int ie = 0; ie < NEBINS; ie++) {
            TChain* ch = new TChain("newtree");
            addFile(ch, base, yr, ebins[ie].suffix);
            Long64_t n = ch->GetEntries();
            std::cout << "  " << ebins[ie].label << ": " << n << " eventi" << std::endl;

            Int_t t_val = 0;
            ch->SetBranchStatus("*", 0);
            ch->SetBranchStatus("time", 1);
            ch->SetBranchAddress("time", &t_val);

            for (Long64_t ev = 0; ev < n; ev++) {
                ch->GetEntry(ev);
                Double_t t_unix = toUnix(t_val);
                h_yr[ie]->Fill(t_unix);
                h_tot[ie]->Fill(t_unix);
            }
            delete ch;
        }

        // Canvas annuale
        TCanvas* c = new TCanvas(Form("c_%d", yr), Form("Anno %d", yr), 1500, 550);
        c->SetLeftMargin(0.08);
        c->SetLogy();
        c->SetRightMargin(0.13);
        c->SetBottomMargin(0.15);

        Double_t ymax = 0;
        for (int ie = 0; ie < NEBINS; ie++)
            if (h_yr[ie]->GetMaximum() > ymax) ymax = h_yr[ie]->GetMaximum();

        TH1D* hframe = (TH1D*)h_yr[0]->Clone(Form("frame_%d", yr));
        hframe->Reset();
        hframe->SetMaximum(ymax * 1.15);
        hframe->SetMinimum(0.5);
        hframe->SetTitle(Form("DAMPE pHe - %d; ;Counts / day", yr));
        setTimeAxis(hframe, false);
        hframe->Draw("AXIS");

        for (int ie = 0; ie < NEBINS; ie++)
            h_yr[ie]->Draw("HIST SAME");

        makeLegend(h_yr, NEBINS)->Draw();
        c->Update();
        c->SaveAs(Form("plot_%d.pdf", yr));
        c->SaveAs(Form("plot_%d.png", yr));
        std::cout << "  Salvato: plot_" << yr << ".pdf" << std::endl;

        for (int ie = 0; ie < NEBINS; ie++) delete h_yr[ie];
        delete hframe;
        delete c;
    }

    // Canvas totale
    std::cout << "\n=== Plot totale ===" << std::endl;
    TCanvas* ctot = new TCanvas("ctot", "Tutti gli anni", 1700, 550);
    ctot->SetLeftMargin(0.07);
    ctot->SetLogy();
    ctot->SetRightMargin(0.13);
    ctot->SetBottomMargin(0.15);

    Double_t ymax_tot = 0;
    for (int ie = 0; ie < NEBINS; ie++)
        if (h_tot[ie]->GetMaximum() > ymax_tot) ymax_tot = h_tot[ie]->GetMaximum();

    TH1D* hframe_tot = (TH1D*)h_tot[0]->Clone("frame_tot");
    hframe_tot->Reset();
    hframe_tot->SetMaximum(ymax_tot * 1.15);
    hframe_tot->SetMinimum(0.5);
    hframe_tot->SetTitle("DAMPE pHe - 2016-2025; ;Counts / week");
    setTimeAxis(hframe_tot, true);
    hframe_tot->Draw("AXIS");

    for (int ie = 0; ie < NEBINS; ie++)
        h_tot[ie]->Draw("HIST SAME");

    makeLegend(h_tot, NEBINS)->Draw();
    ctot->Update();
    ctot->SaveAs("plot_total_20260608.pdf");
    ctot->SaveAs("plot_total_20260608.png");
    std::cout << "Salvato: plot_total.pdf" << std::endl;

    std::cout << "\nFatto!" << std::endl;
}