/*
 * diagnose_tree.C
 * ===============
 * Ispeziona il primo file disponibile:
 *   - lista tutti i branch del TTree "newtree"
 *   - stampa tipo e range del branch "time" (o quello trovato)
 *   - stampa i primi 10 valori
 *
 * Utilizzo:
 *   root -l -b -q diagnose_tree.C
 */
/*
void diagnose_tree() {

    TString base;
    TString host = gSystem->HostName();
    if (host.Contains("cnaf"))    base = "/storage/gpfs_data/dampe/users";
    else if (host.Contains("le.infn")) base = "/nfs/argo/dampe";
    else { std::cout << "Hostname non riconosciuto: " << host << std::endl; return; }

    // Apri il primo file disponibile
    TString fpath = base + "/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2016_merged.root";
    TFile* f = TFile::Open(fpath, "READ");
    if (!f || f->IsZombie()) {
        std::cout << "File non trovato: " << fpath << std::endl;
        return;
    }
    std::cout << "File aperto: " << fpath << "\n" << std::endl;

    TTree* tree = (TTree*)f->Get("newtree");
    if (!tree) { std::cout << "TTree 'newtree' non trovato!" << std::endl; return; }

    std::cout << "Entries totali: " << tree->GetEntries() << "\n" << std::endl;

    // ── 1) Lista tutti i branch con tipo ─────────────────────────────────────
    std::cout << "=== Branch del TTree ===" << std::endl;
    TObjArray* branches = tree->GetListOfBranches();
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch* br = (TBranch*)branches->At(i);
        TLeaf*   lf = br->GetLeaf(br->GetName());
        std::cout << "  " << br->GetName()
                  << "  \t[tipo: " << (lf ? lf->GetTypeName() : "?") << "]"
                  << std::endl;
    }

    // ── 2) Cerca branch con "time" nel nome (case-insensitive) ───────────────
    std::cout << "\n=== Branch contenenti 'time' nel nome ===" << std::endl;
    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch* br = (TBranch*)branches->At(i);
        TString bname = br->GetName();
        if (bname.Contains("time", TString::kIgnoreCase) ||
            bname.Contains("Time", TString::kIgnoreCase) ||
            bname.Contains("TIME", TString::kIgnoreCase)) {
            TLeaf* lf = br->GetLeaf(br->GetName());
            std::cout << "  TROVATO: " << bname
                      << "  [tipo: " << (lf ? lf->GetTypeName() : "?") << "]"
                      << std::endl;
        }
    }

    // ── 3) Prova a leggere "time" come Double_t e come Float_t ───────────────
    std::cout << "\n=== Test lettura branch 'time' ===" << std::endl;

    // Prova Double_t
    Double_t t_double = 0;
    Int_t    res_d = tree->SetBranchAddress("time", &t_double);
    if (res_d >= 0) {
        std::cout << "SetBranchAddress con Double_t: OK" << std::endl;
        tree->GetEntry(0);
        std::cout << "  Primo valore (Double_t): " << t_double << std::endl;
        if (t_double > 1e9 && t_double < 2e9)
            std::cout << "  --> Sembra un Unix timestamp valido!" << std::endl;
        else
            std::cout << "  --> ATTENZIONE: valore fuori range Unix timestamp (1e9-2e9)" << std::endl;
    } else {
        std::cout << "SetBranchAddress con Double_t: FALLITO (res=" << res_d << ")" << std::endl;
    }

    // Reset e prova Float_t
    tree->ResetBranchAddresses();
    Float_t t_float = 0;
    Int_t res_f = tree->SetBranchAddress("time", &t_float);
    if (res_f >= 0) {
        tree->GetEntry(0);
        std::cout << "  Primo valore (Float_t):  " << t_float << std::endl;
    }

    // Reset e prova Int_t
    tree->ResetBranchAddresses();
    Int_t t_int = 0;
    Int_t res_i = tree->SetBranchAddress("time", &t_int);
    if (res_i >= 0) {
        tree->GetEntry(0);
        std::cout << "  Primo valore (Int_t):    " << t_int << std::endl;
    }

    // ── 4) Stampa i primi 10 valori con tree->Draw (metodo infallibile) ──────
    std::cout << "\n=== Primi 10 valori via tree->Draw ===" << std::endl;
    tree->Scan("time", "", "", 10, 0);

    // ── 5) GetMinimum / GetMaximum ────────────────────────────────────────────
    std::cout << "\n=== Range del branch 'time' ===" << std::endl;
    std::cout << "  Min: " << tree->GetMinimum("time") << std::endl;
    std::cout << "  Max: " << tree->GetMaximum("time") << std::endl;

    f->Close();
    std::cout << "\nDiagnosi completata." << std::endl;
}
*/
void diagnose_tree() {
    TString base;
    TString host = gSystem->HostName();
    if (host.Contains("cnaf"))         base = "/storage/gpfs_data/dampe/users";
    else if (host.Contains("le.infn")) base = "/nfs/argo/dampe";
    else { std::cout << "Hostname non riconosciuto" << std::endl; return; }

    TFile* f = TFile::Open(base + "/SKIM_2026_pHe/FLIGHT/skim_flight_002_010_2016_merged.root");
    TTree* tree = (TTree*)f->Get("newtree");

    const char* bname = "BGO_EnergyG_SatCorr_ML_ions_v3";

    // Tipo del branch
    TBranch* br = tree->GetBranch(bname);
    TLeaf*   lf = br ? br->GetLeaf(bname) : nullptr;
    std::cout << "\nBranch '" << bname << "':" << std::endl;
    std::cout << "  Tipo: " << (lf ? lf->GetTypeName() : "NON TROVATO") << std::endl;

    // Primi 10 valori via Scan (infallibile, non dipende dal tipo)
    std::cout << "\nPrimi 10 valori via Scan:" << std::endl;
    tree->Scan(Form("time:%s", bname), "", "", 10, 0);

    // Min/Max
    std::cout << "\nMin: " << tree->GetMinimum(bname) << std::endl;
    std::cout << "Max: " << tree->GetMaximum(bname) << std::endl;

    f->Close();
}