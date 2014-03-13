#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >          LorentzED;
using namespace ROOT::Math::VectorUtil;

bool pcp = true;

int MakeTree( int iSample = 0){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/KinVariables/";

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
  
  const int NMtCut = 4;
  const int NNjetCut = 2;
  const int NTopnessCut = 4;
  const int NMt2wCut = 3;
  const int NYCut = 4;
  const int NDphiCut = 3;
  const int NDrlblCut = 3;
  const int NDrlbgCut = 3;
  const int NChiCut = 3;
  const int NMetCut = 4;
  const int NSR = NMtCut * NNjetCut * NTopnessCut * NMt2wCut * NYCut * NDphiCut * NDrlblCut * NDrlbgCut * NMetCut;

  Float_t MtCut[NMtCut]           = {100., 120., 150., 250.};
  Float_t NjetCut[NNjetCut]       = { 3., 4.};
  Float_t TopnessCut[NTopnessCut] = {-15., 3., 5., 8.};
  Float_t Mt2wCut[NMt2wCut]       = {0., 200., 250.};
  Float_t YCut[NYCut]             = {0., 8., 10., 13.};
  Float_t DphiCut[NDphiCut]       = {0., 0.8, 1.};
  Float_t DrlblCut[NDrlblCut]     = {5., 2., 1.5};
  Float_t DrlbgCut[NDrlbgCut]     = {0., 1.5, 2.};
  Float_t ChiCut[NDrlbgCut]       = {9999999., 5., 3.};
  Float_t MetCut[NMetCut]         = {0., 150., 200., 350.};

  Float_t nEvents[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NMetCut] = {};
  Float_t nEvents2[NMtCut][NNjetCut][NTopnessCut][NMt2wCut][NYCut][NDphiCut][NDrlblCut][NDrlbgCut][NMetCut] = {};

  
  for ( int imt = 0; imt < NMtCut ; imt++){
    for ( int ijet = 0; ijet < NNjetCut ; ijet++){
      for ( int itopness = 0; itopness < NTopnessCut ; itopness++){
	for ( int imt2w = 0; imt2w < NMt2wCut ; imt2w++){
	  for ( int iy = 0; iy < NYCut ; iy++){
	    for ( int idphi = 0; idphi < NDphiCut ; idphi++){
	      for ( int idrlbl = 0; idrlbl < NDrlblCut ; idrlbl++){
		for ( int idrlbg = 0; idrlbg < NDrlbgCut ; idrlbg++){
		  //for ( int ichi = 0; ichi < NChiCut ; ichi++){
		    for ( int imet = 0; imet < NMetCut ; imet++){	     
		      nEvents[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][imet] = 0.;
		      nEvents2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][imet] = 0.;
		      
		    }
		    //}
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
  
  Float_t xs = 0.;
  Float_t Ntries = 0.;

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
  Char_t kinRegion = false;

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  
  cout<<"here0b"<<endl;
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
  tree->SetBranchAddress( "NEvents", &Ntries);
  tree->SetBranchAddress( "xs", &xs);

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
  tree->SetBranchAddress("kinRegion",&kinRegion);


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

  Float_t nEvents0 = 0.;
  Float_t nEvents20 = 0.;

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
  
  outTree->Branch("nEvents", &nEvents0, "nEvents/F");
  outTree->Branch("nEvents2", &nEvents20, "nEvents2/F");


  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 

  cout<<"here1"<<endl;

  for (int ievt=0;ievt<N;++ievt){
    if (ievt%13454 == 0) {
      cout<<"Event number "<<ievt<<"\r"<<flush;
    }
    
    tree->GetEntry(ievt);

    weight = globalWeight * lumi;
    if ( iSample == 8)
      weight /= 10.;

    cout<<weight<<" "<<xs<<" "<<Ntries<<endl;

    if (!kinRegion) continue;    

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
		    
		    //for ( int ichi = 0; ichi < NChiCut ; ichi++){
		    //if (hadChi2 > ChiCut[ichi]) continue;

		      for ( int imet = 0; imet < NMetCut ; imet++){	     
			if (phiCorrMet < MetCut[imet]) continue;
		  
			nEvents[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][imet] += weight;
			nEvents2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][imet] += weight * weight;
		      }
		      //}
		  }
		}
	      }
	    }
	  } 
	}
      }
    }
  }

  cout<<nEvents[0][0][0][1][3][2][0][2][2]<<" "<<nEvents[0][1][0][1][3][2][0][2][2]<<endl;

  /////////////////////////////////////////////////////
  //  WriteHisto
  /////////////////////////////////////////////////////
  
  /*  baseDir->cd();

  lpth->Write();
  letah->Write();
  lRelIsoh->Write();
  
  njetsh->Write();
  jet1h->Write();
  jet2h->Write();
  jet3h->Write();
  jet4h->Write();
    
  nbjetsh->Write();
  bjet1h->Write();
  bjetHighDh->Write();
  
  meth->Write();
  
  hth->Write();
  htRatioh->Write();
  meffh->Write();
  yh->Write();
  
  mth->Write();
  mlb1h->Write();
  mlbh->Write();
  m3bh->Write();
  mt2wh->Write();
  hadChi2h->Write();
  topnessh->Write();
  
  dphiminh->Write();
  drlb1h->Write();
  drlbminh->Write();
  */

  Int_t count = 0.;

  for ( int imt = 0; imt < NMtCut ; imt++){
    for ( int ijet = 0; ijet < NNjetCut ; ijet++){
      for ( int itopness = 0; itopness < NTopnessCut ; itopness++){
	for ( int imt2w = 0; imt2w < NMt2wCut ; imt2w++){
	  for ( int iy = 0; iy < NYCut ; iy++){
	    for ( int idphi = 0; idphi < NDphiCut ; idphi++){
	      for ( int idrlbl = 0; idrlbl < NDrlblCut ; idrlbl++){
		for ( int idrlbg = 0; idrlbg < NDrlbgCut ; idrlbg++){
		  //for ( int ichi = 0; ichi < NChiCut ; ichi++){
		    for ( int imet = 0; imet < NMetCut ; imet++){	     
		      
		      MtCut0 = MtCut[imt];
		      NjetCut0 = NjetCut[ijet];    
		      TopnessCut0 = TopnessCut[itopness];
		      Mt2wCut0 = Mt2wCut[imt2w];
		      YCut0 = YCut[iy];
		      DphiCut0 =  DphiCut[idphi];
		      DrlblCut0 =  DrlblCut[idrlbl];
		      DrlbgCut0 =  DrlbgCut[idrlbg];
		      //ChiCut0 =  ChiCut[ichi];
		      MetCut0 = MetCut[imet];
		      
		      
		      TString tmp = "MT-"; tmp += MtCut[imt]; tmp += "_";
		      tmp += "Njet-"; tmp += NjetCut[ijet]; tmp += "_";
		      tmp += "Topness-"; tmp += TopnessCut[itopness]; tmp += "_";
		      tmp += "Mt2w-"; tmp += Mt2wCut[imt2w]; tmp += "_";
		      tmp += "Y-"; tmp += YCut[iy]; tmp += "_";
		      tmp += "Dphi-"; tmp += DphiCut[idphi]; tmp += "_";
		      tmp += "Drlbl-"; tmp += DrlblCut[idrlbl]; tmp += "_";
		      tmp += "Drlbg-"; tmp += DrlbgCut[idrlbg]; tmp += "_";
		      //tmp += "Chi-"; tmp += ChiCut[ichi]; tmp += "_";
		      tmp += "Met-"; tmp += MetCut[imet];
		      
		      SR0 = tmp;
		      
		      nEvents0 = nEvents[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][imet];
		      nEvents20 = nEvents2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][imet];
		      
		      if ( nEvents0 < 3.) count++;
		      
		      outTree->Fill();
		    }
		    //}
		}
	      }
	    }
	  }
	}
      }
    }
  }
  cout<<endl;
  cout<<count<<endl;


  outFile->cd();
  outTree->Write();

  outFile->Close();    
  inFile->Close();
}
