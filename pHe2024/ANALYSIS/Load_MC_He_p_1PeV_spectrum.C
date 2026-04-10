{
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <string>
//
//

// =============================

const int nsetHe = 6; //number of different energy intervals used to produce the MC data
const int nsetP = 5;

//== the denominator is the total event number before any selection for each energy bin (the generated event number) 
Double_t nwP[nsetP] = {
		  (1./2068200000.)*log(100.),           //1-100 GeV//      
		  (1./1056770000.)*log(10.),          //100GeV-1TeV//
                  (1./252596200.)*log(10.),             //1TeV-10TeV//
		  (1./117399700.)*log(10.),                 //10TeV-100TeV
		  (1./15121700.)*log(10.)                 //100TeV-1PeV
			 };
//
Double_t nwHe[nsetHe] = {(1./204200000.)*log(10.),          //10-100 GeV//
                    (1./101500000.)*log(10.),          //100GeV-1TeV//
                    (1./99545000.)*log(10.),          //1TeV -10TeV//
                    (1./99446452.)*log(10.),          //10TeV-100TeV//
			(1./18673400.)*log(5.),			//100 TeV - 500TeV
			(1./6102260.)*log(2.)			//500TeV - 1 PeV
		   };


//
TChain *skim= new TChain("newtree");
TChain *sk_he[nsetHe];
TChain *sk_p[nsetP];
for (int i=0; i< nsetP; i++) { sk_p[i] = new TChain("newtree"); };
for (int i=0; i< nsetHe; i++) { sk_he[i] = new TChain("newtree"); };


// -p-
// 1 GeV - 100 GeV --- 14 files
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part1.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part2.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part3.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part4.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part5.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part6.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part7.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part8.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part9.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part10.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part11.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part12.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part13.root");
sk_p[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_1GeV_100GeV_FTFP_HP_part14.root");

// 100 GeV - 1 TeV --- 177 files 
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part9.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part8.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part89.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part88.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part87.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part86.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part85.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part84.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part83.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part82.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part81.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part80.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part7.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part79.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part78.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part77.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part76.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part75.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part74.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part73.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part72.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part71.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part70.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part6.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part69.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part68.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part67.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part66.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part65.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part64.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part63.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part62.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part61.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part60.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part5.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part59.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part58.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part57.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part56.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part55.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part54.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part53.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part52.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part51.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part50.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part4.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part49.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part48.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part47.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part46.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part45.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part44.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part43.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part42.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part41.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part40.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part3.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part39.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part38.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part37.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part36.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part35.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part34.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part33.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part32.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part31.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part30.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part2.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part29.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part28.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part27.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part26.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part25.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part24.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part23.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part22.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part21.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part20.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part1.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part19.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part18.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part17.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part16.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part15.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part14.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part13.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part12.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part11.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_part10.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part87.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part86.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part85.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part84.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part83.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part82.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part81.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part80.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part79.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part78.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part77.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part76.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part75.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part74.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part73.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part72.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part71.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part70.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part69.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part68.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part67.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part66.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part65.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part64.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part63.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part62.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part61.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part60.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part59.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part58.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part57.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part56.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part55.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part54.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part53.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part52.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part51.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part50.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part49.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part48.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part47.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part46.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part45.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part44.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part43.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part42.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part41.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part40.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part39.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part38.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part37.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part36.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part35.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part34.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part33.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part32.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part31.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part30.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part29.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part28.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part27.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part26.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part25.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part24.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part23.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part22.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part21.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part20.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part19.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part18.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part17.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part16.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part15.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part14.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part13.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part12.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part11.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part10.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part09.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part08.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part07.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part06.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part05.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part04.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part03.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part02.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part01.root");
sk_p[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/PROTON_FTFP_v6r0p10_SATcor/skim_MC_allProton-v6r0p10_100GeV_1TeV_FTFP_HP_p2_part00.root");

// 1 TeV - 10 TeV 
for (int i=0; i<10; i++){
    sk_p[2]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[2]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP_part%i.root",i));
}
for (int i=0; i<10; i++){
    sk_p[2]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP-p2_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[2]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP-p2_part%i.root",i));
}


// 10 TeV - 100 TeV 
for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP_part%i.root",i));
}

for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p1_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p1_part%i.root",i));
}

for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p3_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p3_part%i.root",i));
}

for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p4_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p4_part%i.root",i));
}

for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p6_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p6_part%i.root",i));
}

for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p7_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p7_part%i.root",i));
}

for (int i=0; i<10; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p8_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_p[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p8_part%i.root",i));
}

// 100 TeV - 1 PeV EU
for (int i=0; i<10; i++){
sk_p[4]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p12p13_100TeV_1PeV_EPOSLHC_FTFP_part0%i.root",i));
}
for (int i=10; i<90; i++){
sk_p[4]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p12p13_100TeV_1PeV_EPOSLHC_FTFP_part%i.root",i));
}

 cout<<" "<<endl;
 cout<<" ...sk_p[0]-> "<<sk_p[0]->GetEntries()<<endl;
 cout<<" ...sk_p[1]-> "<<sk_p[1]->GetEntries()<<endl;
 cout<<" ...sk_p[2]-> "<<sk_p[2]->GetEntries()<<endl;
 cout<<" ...sk_p[3]-> "<<sk_p[3]->GetEntries()<<endl;
 cout<<" ...sk_p[4]-> "<<sk_p[4]->GetEntries()<<endl;
 cout<<" "<<endl;

// N_Gen
Double_t ngen_P[nsetP]={0.0,0.0,0.0,0.0,0.0};
Double_t ngen_P_w[nsetP]={0.0,0.0,0.0,0.0,0.0}; 
Int_t nfiles0=14, nfiles1=177, nfiles2=188,nfiles3=683, nfiles4=90;
TFile *Pfile0[nfiles0],*Pfile1[nfiles1],*Pfile2[nfiles2],*Pfile3[nfiles3],*Pfile4[nfiles4];
TH1 *hP0[nfiles0], *hP1[nfiles1],*hP2[nfiles2],*hP3[nfiles3],*hP4[nfiles4];

TH1*hPTot2;
TH1*hPTot3;
TH1*hPTot4;
/*
for (int i=0; i<10; i++) {
 Pfile2[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP_part0%i.root",i));
 Pfile2[i]->GetObject("h_energy_truth_cut_0",hP2[i]);
}
for (int i=10; i<89; i++) {
 Pfile2[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP_part%i.root",i));
 Pfile2[i]->GetObject("h_energy_truth_cut_0",hP2[i]);
}
for (int i=89; i<99; i++) {
 Pfile2[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP-p2_part0%i.root",i-89));
 Pfile2[i]->GetObject("h_energy_truth_cut_0",hP2[i]);
}
for (int i=99; i<178; i++) {
 Pfile2[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_1TeV_10TeV_FTFP-p2_part%i.root",i-89));
 Pfile2[i]->GetObject("h_energy_truth_cut_0",hP2[i]);
}

for (int i=0; i<178; i++) {
ngen_P[2] += hP2[i]->GetEntries();
}
//
cout<<"Number of entries 1 TeV - 10 TeV: "<< setprecision(15)<<ngen_P[2]<<endl;
//****************************

for (int i=0; i<10; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP_part0%i.root",i));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=10; i<89; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP_part%i.root",i));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=89; i<99; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p1_part0%i.root",i-89));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=99; i<178; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p1_part%i.root",i-89));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=178; i<188; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p3_part0%i.root",i-178));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=188; i<267; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p3_part%i.root",i-178));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=267; i<277; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p4_part0%i.root",i-267));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=277; i<356; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p4_part%i.root",i-267));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=356; i<366; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p6_part0%i.root",i-356));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=366; i<445; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p6_part%i.root",i-356));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=445; i<455; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p7_part0%i.root",i-445));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=455; i<534; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p7_part%i.root",i-445));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=534; i<544; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p8_part0%i.root",i-534));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}
for (int i=544; i<623; i++) {
 Pfile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p10_10TeV_100TeV_FTFP-p8_part%i.root",i-534));
 Pfile3[i]->GetObject("h_energy_truth_cut_0",hP3[i]);
}


for (int i=0; i<623; i++) {
ngen_P[3] += hP3[i]->GetEntries();
}
//
cout<<"Number of entries 10 TeV - 100 TeV: "<< setprecision(15)<<ngen_P[3]<<endl;

//****************************
for (int i=0; i<10; i++) {
 Pfile4[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p12p13_100TeV_1PeV_EPOSLHC_FTFP_part0%i.root",i));
 Pfile4[i]->GetObject("h_energy_truth_cut_0",hP4[i]);
}
for (int i=10; i<90; i++) {
 Pfile4[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_PROTON/skim_MC_allProton-v6r0p12p13_100TeV_1PeV_EPOSLHC_FTFP_part%i.root",i));
 Pfile4[i]->GetObject("h_energy_truth_cut_0",hP4[i]);
}


for (int i=0; i<nfiles4; i++) {
ngen_P[4] += hP4[i]->GetEntries();
}
//
cout<<"Number of entries EPOS P: "<< setprecision(15)<<ngen_P[4]<<endl;
// -------------------------------------------------------------
*/
//``````````````````````````````````````````````````````````````````````````
// HE 10 GeV - 100 GeV
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part1.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part2.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part3.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part4.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part5.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part6.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part7.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part8.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part9.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part10.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part11.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part12.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part13.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part14.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part15.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part16.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part17.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part18.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part19.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part20.root");  //
sk_he[0]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_10GeV_100GeV_FTFP_part21.root");  //
//
// HE 100 GeV - 1 TeV
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part1.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part2.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part3.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part4.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part5.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part6.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part7.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part8.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part9.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part10.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part11.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part12.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part13.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part14.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part15.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part16.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part17.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part18.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part19.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part20.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part21.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part22.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part23.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part24.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part25.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part26.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part27.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part28.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part29.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part30.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part31.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part32.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part33.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part34.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part35.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part36.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part37.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part38.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part39.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part40.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part41.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part42.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part43.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part44.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part45.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part46.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part47.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part48.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part49.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part50.root");  //
sk_he[1]->Add("/nfs/argo/dampe/SkimMC_v6r0p10/HELIUM_FTFP_v6r0p10/allHe4_v6r0p10_100GeV_1TeV_FTFP_part51.root");  // 
//
// HE 1 TeV - 10 TeV with BGO saturation correction
// 1 TeV - 10 TeV 
for (int i=0; i<10; i++){
    sk_he[2]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_1TeV_10TeV_FTFP_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_he[2]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_1TeV_10TeV_FTFP_part%i.root",i));
}
//
// HE 10 TeV - 100 TeV with BGO saturation correction
for (int i=0; i<10; i++){
    sk_he[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_10TeV_100TeV_FTFP_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_he[3]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_10TeV_100TeV_FTFP_part%i.root",i));
}
//--------------------------------------------------------------
// HE 100 TeV - 500 TeV with BGO saturation correction
for (int i=0; i<10; i++){
    sk_he[4]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_100TeV_500TeV_EPOSLHC_part0%i.root",i));
}
for (int i=10; i<89; i++){
    sk_he[4]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_100TeV_500TeV_EPOSLHC_part%i.root",i));
}
// --------------------------------------------------------------------------
//
// 500 TeV - 1 PeV PMO
for (int i=0; i<10; i++){
sk_he[5]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP_SATcor_part0%i.root",i));
}
for (int i=10; i<90; i++){
sk_he[5]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP_SATcor_part%i.root",i));
}
// 500 TeV - 1 PeV EU
for (int i=0; i<10; i++){
sk_he[5]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP-p1_SATcor_part0%i.root",i));
}
for (int i=10; i<90; i++){
sk_he[5]->Add(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP-p1_SATcor_part%i.root",i));
}

//
// --------------------------------------------------------------------------
//

cout<<" "<<endl;
cout<<" ...sk_he[0]-> "<<sk_he[0]->GetEntries()<<endl;
cout<<" ...sk_he[1]-> "<<sk_he[1]->GetEntries()<<endl;
cout<<" ...sk_he[2]-> "<<sk_he[2]->GetEntries()<<endl;
cout<<" ...sk_he[3]-> "<<sk_he[3]->GetEntries()<<endl;
cout<<" ...sk_he[4]-> "<<sk_he[4]->GetEntries()<<endl;
cout<<" ...sk_he[5]-> "<<sk_he[5]->GetEntries()<<endl;
cout<<" "<<endl;
//
// -----------------------------------------------------------------------------------------
// ------------------------------------------------------------------- Evaluate Ngen:
//
Double_t ngen_He[6]={0.0,0.0,0.0,0.0,0.0,0.0};
Double_t ngen_He_w[6]={0.0,0.0,0.0,0.0,0.0,0.0};
TFile *Hefile0[21],*Hefile1[51],*Hefile2[89],*Hefile3[89],*Hefile4[89], *Hefile5[180];
TH1       *he0[21],    *he1[51],    *he2[89],    *he3[89],    *he4[89], *he5[180];
TH1     *he0_w[21],  *he1_w[51],  *he2_w[89],  *he3_w[89],  *he4_w[89], *he5_w[180];
//
for (int i=0; i<10; i++) {
 Hefile2[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_1TeV_10TeV_FTFP_part0%i.root",i));
 Hefile2[i]->GetObject("h_energy_truth_cut_0",he2[i]);
}
for (int i=10; i<89; i++) {
 Hefile2[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_1TeV_10TeV_FTFP_part%i.root",i));
 Hefile2[i]->GetObject("h_energy_truth_cut_0",he2[i]);
}
//
for (int i=0; i<89; i++) {
 ngen_He[2] += he2[i]->GetEntries();
}
//
cout<<"Number of entries He 1-10 TeV: "<<setprecision(15)<<ngen_He[2]<<endl;
//*********************************
for (int i=0; i<10; i++) {
 Hefile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_10TeV_100TeV_FTFP_part0%i.root",i));
 Hefile3[i]->GetObject("h_energy_truth_cut_0",he3[i]);
}
for (int i=10; i<89; i++) {
 Hefile3[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_10TeV_100TeV_FTFP_part%i.root",i));
 Hefile3[i]->GetObject("h_energy_truth_cut_0",he3[i]);
}
//
for (int i=0; i<89; i++) {
 ngen_He[3] += he3[i]->GetEntries();
}
//
cout<<"Number of entries He 10-100 TeV: "<<setprecision(15)<<ngen_He[3]<<endl;
//*********************************
for (int i=0; i<10; i++) {
 Hefile4[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_100TeV_500TeV_EPOSLHC_part0%i.root",i));
 Hefile4[i]->GetObject("h_energy_truth_cut_0",he4[i]);
}
for (int i=10; i<89; i++) {
 Hefile4[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4-v6r0p10_100TeV_500TeV_EPOSLHC_part%i.root",i));
 Hefile4[i]->GetObject("h_energy_truth_cut_0",he4[i]);
}
//
for (int i=0; i<89; i++) {
 ngen_He[4] += he4[i]->GetEntries();
}
//
cout<<"Number of entries He 100-500 TeV: "<<setprecision(15)<<ngen_He[4]<<endl;
// -------------------------------------------------------------
//

for (int i=0; i<10; i++) {
 Hefile5[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP_SATcor_part0%i.root",i));
 Hefile5[i]->GetObject("h_energy_truth_cut_0",he5[i]);
}
for (int i=10; i<90; i++) {
 Hefile5[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP_SATcor_part%i.root",i));
 Hefile5[i]->GetObject("h_energy_truth_cut_0",he5[i]);
}
for (int i=90; i<100; i++) {
 Hefile5[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP-p1_SATcor_part0%i.root",i-90));
 Hefile5[i]->GetObject("h_energy_truth_cut_0",he5[i]);
}
for (int i=100; i<180; i++) {
 Hefile5[i]= new TFile(Form("/nfs/argo/dampe/v6r0p13/BGO_SAT_COR/MC_HELIUM/skim_MC_allHe4_v6r0p13_500TeV_1PeV_FTFP-p1_SATcor_part%i.root",i-90));
 Hefile5[i]->GetObject("h_energy_truth_cut_0",he5[i]);
}

//
for (int i=0; i<180; i++) {
 ngen_He[5] += he5[i]->GetEntries();
}
//
cout<<"Number of entries PeV He: "<<setprecision(15)<<ngen_He[5]<<endl;
// -------------------------------------------------------------
//
// -------------------------------------------------------------------------------
//


//== here we have a normalization factor (SAME NUMBERS, as before)


TCut wP[nsetP] = {"(1./2068200000.)*log(100.)",           //1-100 GeV//      
                  "(1./1056770000.)*log(10.)",          //100GeV-1TeV//
                  "(1./252596200.)*log(10.)" ,            //1TeV-10TeV//
		  "(1./117399700.)*log(10.)",               //10TeV-100TeV
		  "(1./15121700.)*log(10.)"               //100TeV-1PeV
                      };
//
TCut wHe[nsetHe] = {"(1./204200000.)*log(10.)",          //10-100 GeV//
                    "(1./101500000.)*log(10.)",          //100GeV-1TeV//
                    "(1./99545000.)*log(10.)",          //1TeV -10TeV//
                    "(1./99446452.)*log(10.)",         //10TeV-100TeV//
		"(1./18673400.)*log(5.)",       //100TeV-500TeV EPOS//
		"(1./6102260.)*log(2.)"			//500TeV - 1 PeV 
		};

TCut wHeN1[nsetHe];
TCut wPN1[nsetP];

TCut wcomP = "2*3.1415*3.1415*1.38*1.38";                                                 // !!!!!!!!!!!!!!!!!!!!!!
TCut wcomHe = "2*3.14159*3.14159";  
TCut GeoCorr = "(2.*3.14159*3.14159)/(2.*3.14159*3.14159*1.38*1.38)"; 

TCut wspectAmsPs =  "(0.4544*(1./45.)**(-2.849)*(MC_EnergyT)**(-1.849)*(1.+(MC_EnergyT/336.)**(5.542))**(0.024))";
TCut wspectAmsHes =  "((0.0948/2.)*(1./(2*45.))**(-2.780)*(MC_EnergyT)**(-1.780)*(1.+(MC_EnergyT/(2*245.))**(4.407))**(0.027))";

for (int i=0; i< nsetP; i++) { wPN1[i] = wP[i]*wspectAmsPs; };
for (int i=0; i< nsetHe; i++) { wHeN1[i] = wHe[i]*wspectAmsHes; };


TCut Trig_HEP="BGO_HET>0.";

TCut cc204s = "(BGO_EnergyG>20.)";

TCut cTrig="(BGO_EnergyG<1.e+5)"*Trig_HEP; 
//
// --------------------------------------------------------------------------
// *** Selection CUTS ***
//
//
TCut cut00  = cc204s*Trig_HEP;
//
TCut cut01="PSD_psdY_Emax_layer>0.0 && PSD_psdX_Emax_layer>0.0";              // un segnale su entrambe i due piani di PSD

// New Cut order:
// >new-cut02:
//TCut cut02="((nbgohit01*ebgotot01+nbgohit02*ebgotot02)*(ebgiotot07+ebgotot06)/((ebgotot01+ebgotot02)*(nbgohit07*ebgotot07+nbgohit06*ebgotot06))<1)";
TCut cut02="BGO_EneLay[0]+BGO_EneLay[1]<BGO_EneLay[2]+BGO_EneLay[3]";

//TCut cut03="Etrack01>0.0 && Etrack23>0.0";                         // un segnale su entrambe i due piani di PSD
TCut cut03ups="(10.*(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_EnergyY_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_EnergyY_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_psdY_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_psdY_pathlength[1]))>0."; // un segnale su entrambe i due piani di PSD
TCut cut03dws="(10.*(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_EnergyX_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_EnergyX_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_psdX_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_psdX_pathlength[1]))>0."; // un segnale su entrambe i due piani di PSD
TCut cut03=cut03ups*cut03dws; 

//TCut cut04="(fabs(sqrt(Etrack01)-sqrt(epsdmax01/facslope))<5. && fabs(sqrt(Etrack23)-sqrt(epsdmax23/facslope))<5.)";
TCut cut04ups="abs(sqrt(10.*(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_EnergyY_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_EnergyY_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyY_corr[0])+1)*PSD_psdY_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyY_corr[1])+1)*PSD_psdY_pathlength[1]))-sqrt(PSD_psdY_Emax_layer*STK_theta_correction))<5.";
TCut cut04dws="abs(sqrt(10.*(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_EnergyX_corr[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_EnergyX_corr[1])/(0.5*(TMath::Sign(1.,PSD_EnergyX_corr[0])+1)*PSD_psdX_pathlength[0]+0.5*(TMath::Sign(1.,PSD_EnergyX_corr[1])+1)*PSD_psdX_pathlength[1]))-sqrt(PSD_psdX_Emax_layer*STK_theta_correction))<5.";
TCut cut04=cut04ups*cut04dws;


//TCut cut05="(fabs(cbgomax01-cbgostk01)<30.) && (fabs(cbgomax02-cbgostk02)<30.)";
TCut cut05="TMath::Abs(BGO_cbgomax[0]-BGO_cbgostk[0])<30. && TMath::Abs(BGO_cbgomax[1]-BGO_cbgostk[1])<30.";


 TCut ctot=cut00*cut01*cut02*cut03*cut04*cut05;




TCut SpCut= "log10(BGO_nhits_layers[12]*BGO_EneLay[12]/BGO_EnergyG+BGO_nhits_layers[13]*BGO_EneLay[13]/BGO_EnergyG)>-0.0001*(BGO_nbars)**2.0-1.34";

 
 
TCut cutEt = "MC_EnergyT>20.";

//quenching correction:
//	1-0.026*(Ene_true/10(GeV/n))^(-0.32)
//	Ene_true/10(GeV/n) == MC_EnergyT/40 

//1-0.026*(MC_EnergyT/40)^(-0.32)

//************** NEW CHARGE SELECTION 06/04/2022 ***********************//

TString HeMcMPVf="((12.69)+(-6.543)*log10(BGO_EnergyG)+(3.776)*log10(BGO_EnergyG)**2+(-0.9158)*log10(BGO_EnergyG)**3+(0.1075)*log10(BGO_EnergyG)**4)";
TString HeMcMPVf_Cor="((10.95)+(-4.175)*log10(BGO_EnergyG)+(2.653)*log10(BGO_EnergyG)**2+(-0.6897)*log10(BGO_EnergyG)**3+(0.07261)*log10(BGO_EnergyG)**4)";
TString HeMPVf="((7.852)+(1.035)*log10(BGO_EnergyG)+(-0.4836)*log10(BGO_EnergyG)**2+(0.1124)*log10(BGO_EnergyG)**3+(-0.001648)*log10(BGO_EnergyG)**4)";
TString HeMcWidthf="((2.04)+(-3.504)*log10(BGO_EnergyG)+(2.481)*log10(BGO_EnergyG)**2+(-0.7159)*log10(BGO_EnergyG)**3+(0.08252)*log10(BGO_EnergyG)**4)";
TString HeMcWidthf_Cor="((0.6497)+(-0.9196)*log10(BGO_EnergyG)+(0.7154)*log10(BGO_EnergyG)**2+(-0.1817)*log10(BGO_EnergyG)**3+(0.01912)*log10(BGO_EnergyG)**4)";
TString HeWidthf="((0.2505)+(0.02188)*log10(BGO_EnergyG)+(0.1554)*log10(BGO_EnergyG)**2+(-0.06573)*log10(BGO_EnergyG)**3+(0.01144)*log10(BGO_EnergyG)**4)";


TString HeMcGSigmaf="(0.5807)+(-1.46e-42)*log10(BGO_EnergyG)";
TString HeMcGSigmaf_Cor="(0.5533)+(-2.406e-22)*log10(BGO_EnergyG)";
TString HeGSigmaf="(0.4789)+(-4.384e-33)*log10(BGO_EnergyG)";	


TString HeMcFSig= "sqrt("+HeMcWidthf+"**2+"+HeMcGSigmaf+"**2)";
TString HeMcFSigCor= "sqrt("+HeMcWidthf_Cor+"**2+"+HeMcGSigmaf_Cor+"**2)";
TString HeFSig= "sqrt("+HeWidthf+"**2+"+HeGSigmaf+"**2)";


TString PMcMPVf="((3.945)+(-3.557)*log10(BGO_EnergyG)+(2.546)*log10(BGO_EnergyG)**2+(-0.7748)*log10(BGO_EnergyG)**3+(0.09923)*log10(BGO_EnergyG)**4)";
TString PMcMPVf_Cor="((1.948)+(-0.1462)*log10(BGO_EnergyG)+(0.3469)*log10(BGO_EnergyG)**2+(-0.1806)*log10(BGO_EnergyG)**3+(0.03454)*log10(BGO_EnergyG)**4)";
TString PMPVf="((2.071)+(-0.3114)*log10(BGO_EnergyG)+(0.4153)*log10(BGO_EnergyG)**2+(-0.1876)*log10(BGO_EnergyG)**3+(0.03373)*log10(BGO_EnergyG)**4)";
TString PMcWidthf="((0.8941)+(-1.455)*log10(BGO_EnergyG)+(1.07)*log10(BGO_EnergyG)**2+(-0.3192)*log10(BGO_EnergyG)**3+(0.03884)*log10(BGO_EnergyG)**4)";
TString PMcWidthf_Cor="((-0.1775)+(0.3985)*log10(BGO_EnergyG)+(-0.09441)*log10(BGO_EnergyG)**2+(-0.02121)*log10(BGO_EnergyG)**3+(0.01006)*log10(BGO_EnergyG)**4)";
TString PWidthf="((0.03606)+(0.09163)*log10(BGO_EnergyG)+(0.05757)*log10(BGO_EnergyG)**2+(-0.04983)*log10(BGO_EnergyG)**3+(0.01135)*log10(BGO_EnergyG)**4)";

TString PMcGSigmaf="(1.056E-09)+(-1.461E-40)*log10(BGO_EnergyG)";
TString PMcGSigmaf_Cor="(8.949E-10)+(-1.461E-40)*log10(BGO_EnergyG)";
TString PGSigmaf="(6.945E-09)+(-1.461E-40)*log10(BGO_EnergyG)";		


//********************** END NEW **********************************//

TString PMcFSig= "sqrt("+PMcWidthf+"**2+"+PMcGSigmaf+"**2)";
TString PMcFSigCor= "sqrt("+PMcWidthf_Cor+"**2+"+PMcGSigmaf_Cor+"**2)";
TString PFSig= "sqrt("+PWidthf+"**2+"+PGSigmaf+"**2)";


//Questo è il valore di PSD energy corretto da usare per i nuovi istogrammi di carica
TString HeliumCor ="(((PSD_Etrack01+PSD_Etrack23)/2-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+")";
TString ProtonCor ="(((PSD_Etrack01+PSD_Etrack23)/2-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+")";

TString p_low_mcCor= "("+PMcMPVf_Cor+"-"+ProtonCor +")<(2*"+PMcFSigCor+")";
TString he_high_mcCor6= "("+HeliumCor+"-"+HeMcMPVf_Cor+")<(6*"+HeMcFSigCor+")";


TString p_low_mc= "("+PMcMPVf+"-"+"(PSD_Etrack01+PSD_Etrack23)/2)<(2*"+PMcFSig+")";
TString he_high_mc6= "((PSD_Etrack01+PSD_Etrack23)/2-"+HeMcMPVf+")<(6*"+HeMcFSig+")";

TString p_low= "("+PMPVf+"-"+"(PSD_Etrack01+PSD_Etrack23)/2)<(2*"+PFSig+")*(((PSD_Etrack01+PSD_Etrack23)/2)>0.)";
TString he_high6= "((PSD_Etrack01+PSD_Etrack23)/2-"+HeMPVf+")<(6*"+HeFSig+")";



TCut PLow_MC = p_low_mc;
TCut HeHigh_MC6= he_high_mc6;

TCut PLow_MCCor = p_low_mcCor;
TCut HeHigh_MCCor6= he_high_mcCor6;

TCut PLow = p_low;
TCut HeHigh6= he_high6;




TCut Cut_PHe_MC6 =PLow_MC*HeHigh_MC6;
TCut Cut_PHe_MCCor6 =PLow_MCCor*HeHigh_MCCor6;
TCut Cut_PHe6 =PLow*HeHigh6;


TCut ctotCont=cc204s*cut00*cut01*cut02*cut03*cut04*cut05*SpCut;
TCut ctotCont6=cut00*cut01*cut02*cut03*cut04*cut05*SpCut*Cut_PHe6;
TCut ctotContMC6=cut00*cut01*cut02*cut03*cut04*cut05*SpCut*Cut_PHe_MC6;
TCut ctotContMCCor6=cut00*cut01*cut02*cut03*cut04*cut05*SpCut*Cut_PHe_MCCor6;


 int nbd = 48;
 int ndec = 6;

 float arg1 = 1./float(nbd);
 int noe = nbd*ndec;
 cout << "Number of energy bin " << noe << endl;
 float e0 = 1.;

 float Ebin[noe+1];
 Ebin[0] = e0;
 for (int i=1; i<noe+1; i++) {
    Ebin[i] = Ebin[i-1]*TMath::Power(10., arg1);
 }


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/*
TCut cutpsd[noe];
 for (int i=0; i<noe; i++) {
    cutpsd[i]=Form("BGO_EnergyG>=%f && BGO_EnergyG<%f",Ebin[i],Ebin[i+1]);
}

 TFile *fout = new TFile("PHe_MC_p_He-1PeV_v6r0p10p12p13_PSD_Histo_Smearing_HighEnergy.root","RECREATE");
 
for(int j =0;j<noe;j++){
//for(int j =16;j<17;j++){

	cout<<"Now in bin:"<<j<<"... "<<endl;

	TH1F* PSDE=new TH1F(Form("PSDE_%d",j),Form("MC_PSD_Charge with BGO energy%.3f-%.3f",Ebin[j],Ebin[j+1]),400,0,40);
	TH1F* PSDE_P=new TH1F(Form("PSDE_P_%d",j),Form("MC_PSD_Charge with BGO energy%.3f-%.3f",Ebin[j],Ebin[j+1]),400,0,40);
	//PSDE_P->SetLineColor(1);
	TH1F* PSDE_He=new TH1F(Form("PSDE_He_%d",j),Form("MC_PSD_Charge with BGO energy%.3f-%.3f",Ebin[j],Ebin[j+1]),400,0,40);
	//PSDE_He->SetLineColor(4);

	PSDE->Sumw2();
	PSDE_P->Sumw2();
	PSDE_He->Sumw2();

	cout<<Form("BGO_EnergyG>=%f && BGO_EnergyG<%f",Ebin[j],Ebin[j+1])<<endl;

	TString nbj= Form(">>+PSDE_%d",j);
	TString nbjP= Form(">>+PSDE_P_%d",j);
	TString nbjHe= Form(">>+PSDE_He_%d",j);

	TString HeliumFillCor ="((PSD_Etrack01+PSD_Etrack23)/2-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+nbj;
	TString HeliumFillCorHe ="((PSD_Etrack01+PSD_Etrack23)/2-"+HeMcMPVf+")*("+HeFSig+"/"+HeMcFSig+")+"+HeMPVf+nbjHe;
	TString ProtonFillCor ="((PSD_Etrack01+PSD_Etrack23)/2-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+nbj;
	TString ProtonFillCorP ="((PSD_Etrack01+PSD_Etrack23)/2-"+PMcMPVf+")*("+PFSig+"/"+PMcFSig+")+"+PMPVf+nbjP;

	cout<<"Begin to Fill...."<<endl;
	
	for (int i=0; i<nsetP; i++) {
	
	sk_p[i]-> Draw(ProtonFillCor,wPN1[i]*ctotCont*cutpsd[j]);
	sk_p[i]-> Draw(ProtonFillCorP,wPN1[i]*ctotCont*cutpsd[j]);
	}

	cout<<"MC_Proton_filling finished...."<<endl;
	
					sk_he[i]-> Draw(HeliumFillCor,wHeN1[i]*ctotCont*cutpsd[j]);
					sk_he[i]-> Draw(HeliumFillCorHe,wHeN1[i]*ctotCont*cutpsd[j]);
					}

	cout<<"MC_Helium_filling finished...."<<endl;
	
	fout->cd();

	PSDE->Write();
	PSDE_P->Write();
	PSDE_He->Write();
	}

 fout->Close();
*/
//~~~~~~~~~~~~~~~~~~~~~~~~~~~

 TH2F *h2NtrigCor6 = new TH2F("h2NtrigCor6", "Ntrig(Eo,Et) - Contained Track FLAT", noe, Ebin, noe, Ebin); 
 TH1F *h1_ForEfficiency = new TH1F("h1_ForEfficiency","MC true energy after selections",noe,Ebin);

h2NtrigCor6->Sumw2();
h1_ForEfficiency->Sumw2();

// =============================


for (int i=0; i<nsetHe; i++) {

//Peso senza AMS 

	sk_he[i]->Draw("((1 - 0.026*(MC_EnergyT/40.)^(-0.32))*(BGO_EnergyG)):(MC_EnergyT)>>+h2NtrigCor6",ctotContMCCor6*wHe[i]*GeoCorr,"colz");
	sk_he[i]->Draw("(MC_EnergyT)>>+h1_ForEfficiency",ctotContMCCor6*wHe[i]*GeoCorr);

     cout<<"He...Cor in:"<<i<<endl;
    }
for (int i=0; i<nsetP; i++) {

//Peso senza AMS 
sk_p[i]->Draw("(BGO_EnergyG):(MC_EnergyT)>>+h2NtrigCor6",ctotContMCCor6*wP[i],"colz");
    sk_p[i]->Draw("(MC_EnergyT)>>+h1_ForEfficiency",ctotContMCCor6*wP[i]);

    cout<<"P...Cor in:"<<i<<endl;
   }


cout<<"h2Ntrig... "<<endl;


 TFile *fout = new TFile("PHe_MC_p_He-1PeV_v6r0p10p12p13.root","RECREATE");

 fout->cd();

 h2NtrigCor6->Write();
 h1_ForEfficiency->Write();

/*
 TFile *fout = new TFile("PHe_MC_p-1PeV_He-1PeV_ChargeSel_AMSWeight.root","RECREATE");
 h2ChargeSel->Write();
 h2ChargeSelAfter->Write();
 */
fout->Close();


cout<<"End script."<<endl;
};



