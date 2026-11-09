
#include <iostream>
#include <vector>
#include <cmath>

#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TMath.h"

using namespace std;

// Info su MC sample
struct MCsample {

    TString tag;
    TString filename;

    double Emin;   // GeV
    double Emax;   // GeV

    double Ngen;
};

// Controlla che due istogrammi abbiano lo stesso binning
bool SameBinning(TH1 *h1, TH1 *h2) {
    if (h1->GetNbinsX() != h2->GetNbinsX())
        return false;

    int nbins = h1->GetNbinsX();

    for (int ibin = 1; ibin <= nbins + 1; ibin++) {

        double x1 = h1->GetXaxis()->GetBinLowEdge(ibin);
        double x2 = h2->GetXaxis()->GetBinLowEdge(ibin);
        double scale = std::max(1.0, std::max(fabs(x1), fabs(x2)));

        if (fabs(x1 - x2) > 1.e-10 * scale)
            return false;
    }

    return true;
}


// Legge l'istogramma da un file e applica il peso del sample
TH1* LoadAndWeightSample(const MCsample &sample, const TString &histname){
    cout << endl;
    cout << "============================================" << endl;
    cout << sample.tag << endl;
    cout << "File: " << sample.filename << endl;

    TFile *fin = TFile::Open(sample.filename, "READ");

    if (!fin || fin->IsZombie()) {
        cerr << "ERROR: cannot open file " << sample.filename << endl;
        return nullptr;
    }

    TH1 *hin = dynamic_cast<TH1*>(fin->Get(histname));

    if (!hin) {
        cerr << "ERROR: histogram " << histname << " not found in " << sample.filename << endl;
        fin->Close();
        delete fin;
        return nullptr;
    }

    // Numero di bin realmente contenuti nell'istogramma
    double Nhist = hin->Integral(0, hin->GetNbinsX()+1);
    cout << "E range       = " << sample.Emin << " - " << sample.Emax << " GeV" << endl;
    cout << "Ngen expected = " << sample.Ngen << endl;
    cout << "Entries       = " << hin->GetEntries() << endl;
    cout << "Integral      = " << Nhist << endl;


    // Controllo Ngen
    double diff = Nhist - sample.Ngen;
    double relDiff = diff / sample.Ngen;
    cout << "Difference    = " << diff << "  (" << 100. * relDiff << " %)" << endl;

    if (fabs(relDiff) > 1.e-5) 
        cout << "WARNING: histogram integral differs " << "from expected Ngen!" << endl;

    // Clono l'istogramma
    TString newname = "hWeighted_" + sample.tag;
    TH1 *h = dynamic_cast<TH1*>(hin->Clone(newname));
    h->SetDirectory(nullptr);

    if (h->GetSumw2N()==0) h->Sumw2();

    // Peso del sample
    // w = ln(Emax/Emin) / Ngen
    double weight = log(sample.Emax / sample.Emin) / sample.Ngen;

    cout << "MC weight     = " << weight << endl;

    // Applico il peso a tutto l'istogramma
    // i bin fuori dal range del sample sono zero, quindi non contribuiscono
    h->Scale(weight);

    fin->Close();
    delete fin;

    return h;

}


// Costruisce la somma di una lista di sample
TH1* BuildTotal(const vector<MCsample> &samples, const TString &histname, const TString &totalName, vector<TH1*> &weightedSamples) {
    TH1 *hTotal = nullptr;


    for (const auto &sample : samples) {
        TH1 *h = LoadAndWeightSample(sample, histname);
        if (!h)
            continue;

        weightedSamples.push_back(h);

        // primo sample: creo l'istogramma totale con lo stesso binning
        if (!hTotal) {
            hTotal = dynamic_cast<TH1*>(h->Clone(totalName));
            hTotal->SetDirectory(nullptr);
            hTotal->Reset("ICES");

            if (hTotal->GetSumw2N() == 0)
                hTotal->Sumw2();
        }

        // controllo il binning
        if (!SameBinning(hTotal, h)) {
            cerr << "ERROR: incompatible binning for " << sample.tag << endl;
            continue;
        }


        // Sommo il sample pesato
        hTotal->Add(h);
    }

    return hTotal;
}


// stampa i contributi al bin che contiene una certa energia
void PrintContributions(double energy, const vector<MCsample> &samples, const vector<TH1*> &weightedSamples) {
    if (weightedSamples.empty())
        return;

    int ibin = weightedSamples[0]->GetXaxis()->FindBin(energy);

    double low = weightedSamples[0]->GetXaxis()->GetBinLowEdge(ibin);
    double high = weightedSamples[0]->GetXaxis()->GetBinUpEdge(ibin);

    cout << endl;
    cout << "============================================" << endl;
    cout << "Contributions around E = " << energy << " GeV" << endl;
    cout << "Bin " << ibin << " = [" << low << ", " << high << "] GeV" << endl;

    double sum = 0.;

    for (unsigned int i = 0; i < weightedSamples.size(); i++) {
        double content = weightedSamples[i]->GetBinContent(ibin);

        if (content == 0.)
            continue;

        cout << samples[i].tag << " : " << content << endl;
        sum += content;
    }

    cout << "SUM = " << sum << endl;
}


// MAIN
void p0_5PeV_5bins_FLUKA()
{
    // ========================================================
    // Nome dell'istogramma presente negli skim
    // ========================================================

    TString histname =
        "h_energy_truth_5bins_cut_0";


    // ========================================================
    // DIRECTORY
    // ========================================================

    TString dirP =
        "/nfs/argo/dampe/SKIM_2026/PROTON/FLUKA/";

    TString dirHe =
        "/nfs/argo/dampe/SKIM_2026/HELIUM/FLUKA/";


    // ========================================================
    // HELIUM
    //
    // Energie espresse tutte in GeV
    // ========================================================

    vector<MCsample> samplesHe = {

        {
            "He_10GeV_100GeV",
            dirHe +
            "He4-v6r0p18_10GeV_100GeV_FLUKA_partial.root",
            1.e1,
            1.e2,
            474272213.
            //198586517.
        },

        {
            "He_100GeV_1TeV",
            dirHe +
            "He4-v6r0p18_100GeV_1TeV_FLUKA_merged.root",
            1.e2,
            1.e3,
            66009860.
        },

        {
            "He_1TeV_10TeV",
            dirHe +
            "He4-v6r0p18_1TeV_10TeV_FLUKA_merged.root",
            1.e3,
            1.e4,
            77281442.
        },

        {
            "He_10TeV_100TeV",
            dirHe +
            "He4-v6r0p18_10TeV_100TeV_FLUKA_merged.root",
            1.e4,
            1.e5,
            19993301.
        },

        {
            "He_100TeV_500TeV",
            dirHe +
            "He4-v6r0p10_100TeV_500TeV_FLUKA_merged.root",
            1.e5,
            5.e5,
            20132200.
        },

        {
            "He_500TeV_1PeV",
            dirHe +
            "He4-v6r0p18_500TeV_1PeV_FLUKA_merged.root",
            5.e5,
            1.e6,
            1246818.
        } ,

        {
            "He_1PeV_5PeV",
            dirHe +
            "He4-v6r0p20_1PeV_5PeV_FLUKA_merged.root",
            1.e6,
            5.e6,
            2263134.
        }
    };


    // ========================================================
    // PROTON
    // ========================================================

    vector<MCsample> samplesP = {

        {
            "P_10GeV_100GeV",
            dirP +
            "Proton_10GeV_100GeV_FLUKA_merged.root",
            1.e1,
            1.e2,
            1055454555.
        },

        {
            "P_100GeV_1TeV",
            dirP +
            "Proton_100GeV_1TeV_FLUKA_merged.root",
            1.e2,
            1.e3,
            579900871.
        },

        {
            "P_1TeV_10TeV",
            dirP +
            "Proton_1TeV_10TeV_FLUKA_merged.root",
            1.e3,
            1.e4,
            274526589.
        },
        
        {
            "P_10TeV_100TeV",
            dirP +
            "Proton_10TeV_100TeV_FLUKA_merged.root",
            1.e4,
            1.e5,
            208226000.
        },

        {
            "P_100TeV_1PeV",
            dirP +
            "Proton_100TeV_1PeV_FLUKA_merged.root",
            1.e5,
            1.e6,
            9999400.
        } ,

        {
            "P_1PeV_5PeV",
            dirP +
            "Proton_1PeV_5PeV_FLUKA_merged.root",
            1.e6,
            5.e6,
            2564400.
        }
    };


    // ========================================================
    // Costruzione istogrammi He e P
    // ========================================================

    vector<TH1*> weightedHe;
    vector<TH1*> weightedP;


    cout << endl;
    cout << "############################################" << endl;
    cout << "                  HELIUM                    " << endl;
    cout << "############################################" << endl;

    TH1 *hHe = BuildTotal( samplesHe, histname, "h1NgenHe", weightedHe);

    cout << endl;
    cout << "############################################" << endl;
    cout << "                  PROTON                    " << endl;
    cout << "############################################" << endl;


    TH1 *hP = BuildTotal( samplesP, histname, "h1NgenP", weightedP);

    if (!hHe || !hP) {
        cerr << "ERROR building total histograms" << endl;
        return;
    }


    // ========================================================
    // P + He
    // ========================================================

    if (!SameBinning(hHe, hP)) {
        cerr << "ERROR: P and He have different binning" << endl;
        return;
    }

    TH1 *hPHe = dynamic_cast<TH1*>(hP->Clone("h1Ngen"));

    hPHe->SetDirectory(nullptr);
    hPHe->Add(hHe);


    // ========================================================
    // DIAGNOSTICA IMPORTANTE
    //
    // Stampiamo proprio il bin problematico:
    // 500 TeV = 5e5 GeV
    // ========================================================

    PrintContributions(5.e5, samplesHe, weightedHe);
    // anche 5 PeV
    PrintContributions(5.e6, samplesHe, weightedHe);


    // ========================================================
    // DRAW
    // ========================================================

    TCanvas *can = new TCanvas("can", "MC weights", 1000, 700);

    can->SetLogx();
    can->SetGridx();
    can->SetGridy();

    hPHe->SetTitle("P + He MC generation weights;True energy [GeV];Weighted MC generation");
    hPHe->Draw("HIST E");


    // ========================================================
    // OUTPUT
    // ========================================================

    TFile *fout = new TFile("PHe_MC_FLUKA_5PeV_h1Ngen_5binsPerDecade.root", "RECREATE");
    fout->cd();

    // Totali
    hHe->Write();
    hP->Write();
    hPHe->Write();

    // Salvo anche i singoli sample pesati:
    // utilissimo per controllare 500 TeV
    for (auto h : weightedHe)
        h->Write();

    for (auto h : weightedP)
        h->Write();

    fout->Close();

    cout << endl;
    cout << "Output written successfully." << endl;
}


