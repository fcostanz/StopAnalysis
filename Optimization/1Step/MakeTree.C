#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

bool pcp = true;

int MakeTree( int iSample = 0){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/home/fcostanz/Bonsai/Optimization/";

  const int NSamples = 11;  
  TString sample[NSamples];
  sample[0] = "DiLep";
  sample[1] = "OneLep";
  sample[2] = "WJets";
  sample[3] = "Rare";
  sample[4] = "T2tb-mStop175mLSP50";
  sample[5] = "T2tb-mStop200mLSP25";
  sample[6] = "T2tb-mStop375mLSP50";
  sample[7] = "T2tb-mStop250mLSP25";
  sample[8] = "T2tb-mStop325mLSP100";
  sample[9] = "T2tb-mStop450mLSP150";
  sample[10] = "T2tb-mStop550mLSP1";
  
  const int NMtCut = 3;
  const int NNjetCut = 2;
  const int NTopnessCut = 3;
  const int NMt2wCut = 4;
  const int NYCut = 4;
  const int NDphiCut = 3;
  const int NDrlblCut = 3;
  const int NDrlbgCut = 3;
  const int NChiCut = 3;
  const int NMetCut = 3;
  const int NSR = NMtCut * NNjetCut * NTopnessCut * NMt2wCut * NYCut * NDphiCut * NDrlblCut * NDrlbgCut * NChiCut * NMetCut;

  Double_t MtCut[NMtCut]           = {100., 120., 200.};
  Double_t NjetCut[NNjetCut]       = { 3., 4.};
  Double_t TopnessCut[NTopnessCut] = {-15., 5., 8.};
  Double_t Mt2wCut[NMt2wCut]       = {0., 200., 250., 300.};
  Double_t YCut[NYCut]             = {0., 8., 10., 13.};
  Double_t DphiCut[NDphiCut]       = {0., 0.8, 1.};
  Double_t DrlblCut[NDrlblCut]     = {5., 2., 1.5};
  Double_t DrlbgCut[NDrlbgCut]     = {0., 1.5, 2.};
  Double_t ChiCut[NChiCut]         = {999999999., 5., 3.};
  Double_t MetCut[NMetCut]         = {150., 250., 350.};

  Float_t tt[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NChiCut][NMetCut] = {};
  Float_t tb[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NChiCut][NMetCut] = {};
  Float_t bb[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NChiCut][NMetCut] = {};
  Float_t tt2[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NChiCut][NMetCut] = {};
  Float_t tb2[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NChiCut][NMetCut] = {};
  Float_t bb2[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NChiCut][NMetCut] = {};
  
  for ( int imt = 0; imt < NMtCut ; imt++){
    for ( int ijet = 0; ijet < NNjetCut ; ijet++){
      for ( int itopness = 0; itopness < NTopnessCut ; itopness++){
	for ( int imt2w = 0; imt2w < NMt2wCut ; imt2w++){
	  for ( int iy = 0; iy < NYCut ; iy++){
	    for ( int idphi = 0; idphi < NDphiCut ; idphi++){
	      for ( int idrlbl = 0; idrlbl < NDrlblCut ; idrlbl++){
		for ( int idrlbg = 0; idrlbg < NDrlbgCut ; idrlbg++){
		  for ( int ichi = 0; ichi < NChiCut ; ichi++){
		    for ( int imet = 0; imet < NMetCut ; imet++){	     
		      tt[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] = 0.;
		      tb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] = 0.;
		      bb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] = 0.;
		      tt2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] = 0.;
		      tb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] = 0.;
		      bb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] = 0.;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  


  Float_t lumi=19500.;
  Float_t weight = 1.;
  /////////////////////////////////////////////////////
  //  Variable Definition
  ///////////////////////////////////////////////////// 

  Float_t globalWeight = 0.;
  Float_t triggerWeight = 0.;
  Float_t puWeight = 0.;
  Float_t isrWeight = 0.;
  Float_t topPtWeight = 0.;

  Float_t lepFromTop = 0.;
  Float_t charginos = 0.;

  Float_t njets = 0.;
  Float_t jet1 = 0.;
  Float_t jet2 = 0.;
  Float_t jet3 = 0.;
  Float_t jet4 = 0.;

  Float_t nbjets = 0.;
  Float_t bjet1 = 0.;
  Float_t bjetHighestDisc = 0.;
  
  Float_t lPt = 0.;
  Float_t lEta = 0.;
  Float_t lRelIso = 0.;

  Float_t isoTrack = 0.;
  Float_t tauVeto = 0.;
  
  Float_t rawmet = 0.;
  Float_t typeImet = 0.;
  Float_t phiCorrMet = 0.;
  
  Float_t ht = 0.;
  Float_t htRatio = 0.;
  Float_t meff = 0.;
  Float_t y = 0.;
  
  Float_t mt = 0.;
  Float_t mlb1 = 0.;
  Float_t mlb = 0.;
  Float_t m3b = 0.;
  Float_t mt2w = 0.;
  Float_t hadChi2 = 0.;
  Float_t topness = 0.;

  Float_t dphimin = 0.;
  Float_t drlb1 = 0.;
  Float_t drlbmin = 0.;

  Int_t pdgIdLep1 = 0;
  Char_t kinRegionFlag = false;
  Char_t searchRegionFlag = false;

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  
  TString inFileName = mainDir; inFileName += sample[iSample]; inFileName +=".root";
  cout<<inFileName<<endl;

  TFile* inFile = new TFile(inFileName,"READ");
  if (!inFile->IsOpen()){
    std::cout<<"not open"<<std::endl;
  }

  TTree* tree;
  tree= (TTree*)inFile->Get("NoSystematic/bonsai");
  //===========================================
  if(pcp)cout<<"inputs set!"<<endl;
  int N = tree->GetEntries();  cout<<"THERE ARE "<<N<<" EVENTS IN "<<inFileName<<endl;
    
  tree->SetBranchAddress( "GlobalWeight", &globalWeight);
  tree->SetBranchAddress( "TriggerWeight", &triggerWeight);
  tree->SetBranchAddress( "PUWeight", &puWeight);
  tree->SetBranchAddress( "isrWeight", &isrWeight);
  tree->SetBranchAddress( "topPtWeight", &topPtWeight);

  tree->SetBranchAddress( "LepFromTop", &lepFromTop);
  tree->SetBranchAddress( "Charginos", &charginos);

  tree->SetBranchAddress( "njets", &njets);
  tree->SetBranchAddress( "jet1", &jet1);
  tree->SetBranchAddress( "jet2", &jet2);
  tree->SetBranchAddress( "jet3", &jet3);
  tree->SetBranchAddress( "jet4", &jet4);

  tree->SetBranchAddress( "nbjets", &nbjets);
  tree->SetBranchAddress( "bjet1", &bjet1);
  tree->SetBranchAddress( "bjetHighestDisc", &bjetHighestDisc);

  tree->SetBranchAddress( "lPt", &lPt);
  tree->SetBranchAddress( "lEta", &lEta);
  tree->SetBranchAddress( "lRelIso", &lRelIso);

  tree->SetBranchAddress( "phiCorrMet", &phiCorrMet);
  
  tree->SetBranchAddress( "ht", &ht);
  tree->SetBranchAddress( "htRatio", &htRatio);
  tree->SetBranchAddress( "meff", &meff);
  tree->SetBranchAddress( "y", &y);

  tree->SetBranchAddress( "mt", &mt);
  tree->SetBranchAddress( "mlb1", &mlb1);
  tree->SetBranchAddress( "mlb", &mlb);
  tree->SetBranchAddress( "m3b", &m3b);
  tree->SetBranchAddress( "mt2w", &mt2w);
  tree->SetBranchAddress( "hadChi2", &hadChi2);
  tree->SetBranchAddress( "topness", &topness);

  tree->SetBranchAddress( "dphimin", &dphimin);
  tree->SetBranchAddress( "drlb1", &drlb1);
  tree->SetBranchAddress( "drlbmin", &drlbmin);

  tree->SetBranchAddress("pdgIdLep1",&pdgIdLep1);
  tree->SetBranchAddress("kinRegion",&kinRegionFlag);
  tree->SetBranchAddress("searchRegionPost",&searchRegionFlag);


  /////////////////////////////////////////////////////
  //  Output Definition
  ///////////////////////////////////////////////////// 

  TFile* outFile = new TFile( sample[iSample]+"_optimization.root", "RECREATE"); 
  TTree* outTree = new TTree("Optimization", "Optimization subTree");


  /////////////////////////////////////////////////////
  //  Branch Definition
  ///////////////////////////////////////////////////// 

  Float_t MtCut0 = 0.;
  Float_t NjetCut0 = 0.;
  Float_t TopnessCut0 = 0.;
  Float_t Mt2wCut0 = 0.;
  Float_t YCut0 = 0.;
  Float_t DphiCut0 = 0.;
  Float_t DrlblCut0 = 0.;
  Float_t DrlbgCut0 = 0.;
  Float_t ChiCut0 = 0.;
  Float_t MetCut0 = 0.;

  TString SR0("");

  Float_t tt0 = 0.;
  Float_t tb0 = 0.;
  Float_t bb0 = 0.;
  Float_t tt20 = 0.;
  Float_t tb20 = 0.;
  Float_t bb20 = 0.;

  outTree->Branch("mtCut", &MtCut0, "mtCut/F");
  outTree->Branch("njetCut", &NjetCut0, "njetCut/F");
  outTree->Branch("topnessCut", &TopnessCut0, "topnessCut/F");
  outTree->Branch("mt2wCut", &Mt2wCut0, "mt2wCut/F");
  outTree->Branch("yCut", &YCut0, "yCut/F");
  outTree->Branch("dphiCut", &DphiCut0, "dphiCut/F");
  outTree->Branch("drlblCut", &DrlblCut0, "drlblCut/F");
  outTree->Branch("drlbgCut", &DrlbgCut0, "drlbgCut/F");
  outTree->Branch("chiCut", &ChiCut0, "ChiCut/F");
  outTree->Branch("metCut", &MetCut0, "metCut/F");
  outTree->Branch("sr", "TString", &SR0);
  
  outTree->Branch("tt", &tt0, "tt/F");
  outTree->Branch("tb", &tb0, "tb/F");
  outTree->Branch("bb", &bb0, "bb/F");
  outTree->Branch("tt2", &tt20, "tt2/F");
  outTree->Branch("tb2", &tb20, "tb2/F");
  outTree->Branch("bb2", &bb20, "bb2/F");

  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 
  for (int ievt=0;ievt<N;++ievt){
    if (ievt%13454 == 0) {
      cout<<"Event number "<<ievt<<"\r"<<flush;
    }
    
    tree->GetEntry(ievt);

    weight = globalWeight * triggerWeight * puWeight * lumi;
    if (iSample < 2) weight *= topPtWeight;
    if (iSample > 3 ) weight *= isrWeight;
    
    if (lRelIso > 0.1) continue;
    if (!searchRegionFlag) continue;

    for ( int imt = 0; imt < NMtCut ; imt++){
      if (mt < MtCut[imt]) continue;

      for ( int ijet = 0; ijet < NNjetCut ; ijet++){
	if (njets < NjetCut[ijet]) continue;

	for ( int itopness = 0; itopness < NTopnessCut ; itopness++){
	  if (topness < TopnessCut[itopness]) continue;
	  
	  for ( int imt2w = 0; imt2w < NMt2wCut ; imt2w++){
	    if (mt2w < Mt2wCut[imt2w]) continue;

	    for ( int iy = 0; iy < NYCut ; iy++){
	      if (y < YCut[iy]) continue;
	      
	      for ( int idphi = 0; idphi < NDphiCut ; idphi++){
		if (dphimin < DphiCut[idphi]) continue;

		for ( int idrlbl = 0; idrlbl < NDrlblCut ; idrlbl++){
		  if (drlb1 > DrlblCut[idrlbl]) continue;

		  for ( int idrlbg = 0; idrlbg < NDrlbgCut ; idrlbg++){
		    if (drlb1 < DrlbgCut[idrlbg]) continue;
		    
		    for ( int ichi = 0; ichi < NChiCut ; ichi++){
		      if (hadChi2 > ChiCut[ichi]) continue;
		      
		      for ( int imet = 0; imet < NMetCut ; imet++){	     
			if (phiCorrMet < MetCut[imet]) continue;
			
			if (charginos == 0){
			  tt[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] += weight;
			  tt2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] += weight * weight;
			}
			if (charginos == 1){
			  tb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] += weight;
			  tb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] += weight * weight;
			}
			if (charginos == 2){
			  bb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] += weight;
			  bb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet] += weight * weight;
			}
		      }
		    }
		  }
		}
	      }
	    }
	  } 
	}
      }
    }
  }

  for ( int imt = 0; imt < NMtCut ; imt++){
    for ( int ijet = 0; ijet < NNjetCut ; ijet++){
      for ( int itopness = 0; itopness < NTopnessCut ; itopness++){
	for ( int imt2w = 0; imt2w < NMt2wCut ; imt2w++){
	  for ( int iy = 0; iy < NYCut ; iy++){
	    for ( int idphi = 0; idphi < NDphiCut ; idphi++){
	      for ( int idrlbl = 0; idrlbl < NDrlblCut ; idrlbl++){
		for ( int idrlbg = 0; idrlbg < NDrlbgCut ; idrlbg++){
		  for ( int ichi = 0; ichi < NChiCut ; ichi++){
		    for ( int imet = 0; imet < NMetCut ; imet++){	     
		      
		      MtCut0 = MtCut[imt];
		      NjetCut0 = NjetCut[ijet];    
		      TopnessCut0 = TopnessCut[itopness];
		      Mt2wCut0 = Mt2wCut[imt2w];
		      YCut0 = YCut[iy];
		      DphiCut0 =  DphiCut[idphi];
		      DrlblCut0 =  DrlblCut[idrlbl];
		      DrlbgCut0 =  DrlbgCut[idrlbg];
		      ChiCut0 =  ChiCut[ichi];
		      MetCut0 = MetCut[imet];
		      
		      
		      TString tmp = "MT-"; tmp += MtCut[imt]; tmp += "_";
		      tmp += "Njet-"; tmp += NjetCut[ijet]; tmp += "_";
		      tmp += "Topness-"; tmp += TopnessCut[itopness]; tmp += "_";
		      tmp += "Mt2w-"; tmp += Mt2wCut[imt2w]; tmp += "_";
		      tmp += "Y-"; tmp += YCut[iy]; tmp += "_";
		      tmp += "Dphi-"; tmp += DphiCut[idphi]; tmp += "_";
		      tmp += "Drlbl-"; tmp += DrlblCut[idrlbl]; tmp += "_";
		      tmp += "Drlbg-"; tmp += DrlbgCut[idrlbg]; tmp += "_";
		      tmp += "Chi-"; tmp += ChiCut[ichi]; tmp += "_";
		      tmp += "Met-"; tmp += MetCut[imet];
		      
		      SR0 = tmp;
		      
		      tt0 = tt[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet];
		      tb0 = tb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet];
		      bb0 = bb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet];
		      tt20 = tt2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet];
		      tb20 = tb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet];
		      bb20 = bb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi][imet];
		      
		      outTree->Fill();
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  outFile->cd();
  outTree->Write();

  outFile->Close();    
  inFile->Close();
  
  return 0;
}
