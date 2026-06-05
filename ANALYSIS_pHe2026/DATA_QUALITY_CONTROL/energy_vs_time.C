/*
 * plot_energy_vs_time.C
 * =====================
 * Energia media giornaliera per range di energia, DAMPE pHe 2016-2025.
 *
 * Produce:
 *   plot_energy_total_all.png      — tutti gli anni, tutti i range sovrapposti
 *   plot_energy_total_YYYY.png     — per anno, tutti i range sovrapposti
 *   plot_energy_SUFFIX_all.png     — per range, tutti gli anni
 *   plot_energy_SUFFIX_YYYY.png    — per range, per anno
 *
 * Utilizzo:
 *   nohup root -l -b -q plot_energy_vs_time.C > plot_energy_vs_time.log 2>&1 &
 */

#include <iostream>

// ── USER CONFIG ───────────────────────────────────────────────────────────────
const int START_YEAR   = 2016;
const int END_YEAR     = 2025;
const int DAY          = 86400;
const Int_t TIME_OFFSET = 1356998050;

struct EnergyBin {
    const char* label;
    const char* suffix;
    Color_t     color;
    Double_t    emin;   // GeV - range Y del TProfile
    Double_t    emax;   // GeV
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
    if (year == 2025 && TString(suffix) == "002_010")
        fname = base + TString::Format(
            "/SKIM_2026_pHe/FLIGHT/skim_flight_%s_%d_partially_merged.root", suffix, year);
    else
        fname = base + TString::Format(
            "/SKIM_2026_pHe/FLIGHT/skim_flight_%s_%d_merged.root", suffix, year);
    ch->Add(fname);
}

Double_t toUnix(Int_t t)      { return (Double_t)t + (Double_t)TIME_OFFSET; }
Double_t yearStartUnix(int yr) { TDatime d(yr,1,1,0,0,0); return (Double_t)d.Convert(kTRUE); }

void setTimeAxis(TAxis* ax, bool shortFmt = false) {
    ax->SetTimeDisplay(1);
    ax->SetTimeOffset(0, "gmt");
    ax->SetTimeFormat(shortFmt ? "%m/%Y" : "%d/%m/%Y");
    ax->SetNdivisions(506);
    ax->SetLabelSize(0.038);
    ax->SetLabelOffset(0.01);
}

// Disegna e salva una canvas con i TProfile elencati in idx[]
void saveCanvas(const char* fname, const char* title,
                TProfile** hh,
                int* idx, int nidx,
                bool shortFmt) {

    TCanvas* c = new TCanvas(fname, title, 1400, 550);
    c->SetLeftMargin(0.10);
    c->SetRightMargin(0.13);
    c->SetBottomMargin(0.15);
    c->SetGrid(0, 1);
    c->SetLogy();

    bool first = true;
    for (int k = 0; k < nidx; k++) {
        TProfile* h = hh[idx[k]];
        if (first) {
            h->SetTitle(Form("%s;Data;Energia media [GeV]", title));
            setTimeAxis(h->GetXaxis(), shortFmt);
            h->GetYaxis()->SetLabelSize(0.038);
            h->GetYaxis()->SetTitleOffset(1.2);
            h->SetMinimum(10.);
            h->SetMaximum(1e5);
            h->Draw("E0");
            first = false;
        } else {
            h->Draw("E0 SAME");
        }
    }

    TLegend* leg = new TLegend(0.88, 0.50, 0.98, 0.88);
    leg->SetTextSize(0.032);
    leg->SetBorderSize(1);
    leg->SetFillStyle(1001);
    for (int k = 0; k < nidx; k++)
        leg->AddEntry(hh[idx[k]], ebins[idx[k]].label, "p");
    leg->Draw();

    c->Update();
    c->SaveAs(fname);
    std::cout << "  Salvato: " << fname << std::endl;
    delete c;
}

// ── MAIN ──────────────────────────────────────────────────────────────────────
void energy_vs_time() {

    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    new TColor(mplBlue,   31./255.,119./255.,180./255.);
    new TColor(mplOrange,255./255.,127./255., 14./255.);
    new TColor(mplGreen,  44./255.,160./255., 44./255.);
    new TColor(mplRed,   214./255., 39./255., 40./255.);
    new TColor(mplPurple,148./255.,103./255.,189./255.);
    new TColor(mplBrown, 140./255., 86./255., 75./255.);

    TString base = getBasePath();
    std::cout << "BasePath: " << base << "\n" << std::endl;

    const int NYEARS = END_YEAR - START_YEAR + 1;

    Double_t t_global_min = yearStartUnix(START_YEAR);
    Double_t t_global_max = yearStartUnix(END_YEAR + 1);
    int n_days_total = (int)((t_global_max - t_global_min) / DAY) + 1;

    // ── Alloca TProfile globali ───────────────────────────────────────────────
    // TProfile(name, title, nbins, xmin, xmax, ymin, ymax)
    // ymin/ymax dicono a ROOT di ignorare Fill() fuori da quel range Y
    TProfile* p_tot[NEBINS];
    for (int ie = 0; ie < NEBINS; ie++) {
        p_tot[ie] = new TProfile(
            Form("p_tot_%d", ie), "",
            n_days_total, t_global_min, t_global_max,
            ebins[ie].emin, ebins[ie].emax);
        p_tot[ie]->SetLineColor(ebins[ie].color);
        p_tot[ie]->SetMarkerColor(ebins[ie].color);
        p_tot[ie]->SetMarkerStyle(7);
        p_tot[ie]->SetMarkerSize(0.5);
    }

    // ── Alloca TProfile per anno ──────────────────────────────────────────────
    TProfile* p_yr[NYEARS][NEBINS];
    for (int iy = 0; iy < NYEARS; iy++) {
        int yr = START_YEAR + iy;
        Double_t t0 = yearStartUnix(yr);
        Double_t t1 = yearStartUnix(yr + 1);
        int nd = (int)((t1 - t0) / DAY);
        for (int ie = 0; ie < NEBINS; ie++) {
            p_yr[iy][ie] = new TProfile(
                Form("p_yr_%d_%d", yr, ie), "",
                nd, t0, t1,
                ebins[ie].emin, ebins[ie].emax);
            p_yr[iy][ie]->SetLineColor(ebins[ie].color);
            p_yr[iy][ie]->SetMarkerColor(ebins[ie].color);
            p_yr[iy][ie]->SetMarkerStyle(7);
            p_yr[iy][ie]->SetMarkerSize(0.5);
        }
    }

    // ── Riempimento ───────────────────────────────────────────────────────────
    for (int iy = 0; iy < NYEARS; iy++) {
        int yr = START_YEAR + iy;
        std::cout << "=== Anno " << yr << " ===" << std::endl;

        for (int ie = 0; ie < NEBINS; ie++) {
            TChain* ch = new TChain("newtree");
            addFile(ch, base, yr, ebins[ie].suffix);
            Long64_t n = ch->GetEntries();
            std::cout << "  " << ebins[ie].label << ": " << n << " eventi" << std::endl;

            Int_t    t_val = 0;
            Double_t e_val = 0.;
            ch->SetBranchStatus("*", 0);
            ch->SetBranchStatus("time", 1);
            ch->SetBranchStatus("BGO_EnergyG_SatCorr_ML_ions_v3", 1);
            ch->SetBranchAddress("time", &t_val);
            ch->SetBranchAddress("BGO_EnergyG_SatCorr_ML_ions_v3", &e_val);

            for (Long64_t ev = 0; ev < n; ev++) {
                ch->GetEntry(ev);
                Double_t t_unix = toUnix(t_val);
                p_yr[iy][ie]->Fill(t_unix, e_val);
                p_tot[ie]->Fill(t_unix, e_val);
            }
            delete ch;
        }
    }

    std::cout << "\n=== Salvo i plot ===" << std::endl;
    int all_idx[NEBINS] = {0, 1, 2, 3, 4, 5};

    // 1) Tutti gli anni, tutti i range sovrapposti
    saveCanvas("MEAN_ENERGY_PLOTS/plot_energy_total_all.png",
               "DAMPE pHe - Energia media giornaliera (2016-2025)",
               p_tot, all_idx, NEBINS, true);

    // 2) Per anno, tutti i range sovrapposti
    for (int iy = 0; iy < NYEARS; iy++) {
        int yr = START_YEAR + iy;
        saveCanvas(Form("MEAN_ENERGY_PLOTS/plot_energy_total_%d.png", yr),
                   Form("DAMPE pHe - Energia media giornaliera %d", yr),
                   p_yr[iy], all_idx, NEBINS, false);
    }

    // 3) Per range, tutti gli anni
    for (int ie = 0; ie < NEBINS; ie++) {
        int one[1] = {ie};
        saveCanvas(Form("MEAN_ENERGY_PLOTS/plot_energy_%s_all.png", ebins[ie].suffix),
                   Form("DAMPE pHe - %s (2016-2025)", ebins[ie].label),
                   p_tot, one, 1, true);
    }

    // 4) Per range, per anno
    for (int ie = 0; ie < NEBINS; ie++) {
        int one[1] = {ie};
        for (int iy = 0; iy < NYEARS; iy++) {
            int yr = START_YEAR + iy;
            saveCanvas(Form("MEAN_ENERGY_PLOTS/plot_energy_%s_%d.png", ebins[ie].suffix, yr),
                       Form("DAMPE pHe - %s - %d", ebins[ie].label, yr),
                       p_yr[iy], one, 1, false);
        }
    }

    std::cout << "\nFatto!" << std::endl;
}