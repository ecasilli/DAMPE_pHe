
void confucio_mask(){
	TFile aa("ROOT_FILES/PHe_MC_p_He_5PeV_unfolding_6binperdecade_2e5sigmaLow_6sigmaUp_new.root");
	TH1D *h2trig_corr = (TH1D*) aa.Get("h2Ntrig_wgt_v3");
	//TH1D *h2trig_corr2 = (TH1D*) aa.Get("h2Ntrig_wgt_v3");


	int nbd = 6; // 6 bins per decade
	int ndec = 6;

	float arg1 = 1./float(nbd);
	int noe = nbd*ndec;
	cout << "Number of energy bin " << noe << endl;
	float e0 = 10.;

	float Ebin[noe+1];
	Ebin[0] = e0;
	for (int i=1; i<noe+1; i++) {
	    Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);
	}


	for(int ybin=1; ybin<noe+1; ybin++){
		for (int xbin=1; xbin<noe+1; xbin++){
			if (ybin<(-9.12+0.8*xbin)){
				h2trig_corr->SetBinContent(xbin,ybin,0);
				//h2trig_corr2->SetBinContent(xbin,ybin,0);
				
				
			}
			if (ybin<3) h2trig_corr->SetBinContent(xbin,ybin,0);
		}
	}

	TCanvas *c1 = new TCanvas("c1","c1",1); c1->SetLogx();  c1->SetLogy();  c1->SetLogz(); 
	h2trig_corr->Draw("colz");

	TFile *fout = new TFile("ROOT_FILES/PHe_MC_p_He_5PeV_unfolding_6binperdecade_2e5sigmaLow_6sigmaUp_new_CORRECTED3.root", "RECREATE");
	fout->cd();
	h2trig_corr->Write();	
	//h2trig_corr2->Write();	
	fout->Close();

}