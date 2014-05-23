#include <vector>

#include "TFile.h"
#include "TTree.h"

int SearchRegions(){
  
  TFile* outFile = new TFile( "SearchRegions.root", "RECREATE");
  outFile->cd();
  TTree* outTree = new TTree( "SearchRegions", "SearchRegions");  

  Int_t ID = 0;
  Float_t nJetCut = 0.;
  Float_t mtCut = 0.;
  Float_t dphiCut = 0.;
  Float_t hadChi2Cut = 0.;
  Float_t metCut = 0.;
  Float_t mt2wCut = 0.;
  Float_t centralityCut = 0.;
  Float_t yCut = 0.;
  Float_t mlbCut = 0.;
  Float_t m3Cut = 0.;

  outTree->Branch( "ID", &ID);
  outTree->Branch( "nJetCut", &nJetCut);
  outTree->Branch( "mtCut", &mtCut);
  outTree->Branch( "dphiCut", &dphiCut);
  outTree->Branch( "hadChi2Cut", &hadChi2Cut);
  outTree->Branch( "metCut", &metCut);
  outTree->Branch( "mt2wCut", &mt2wCut);
  outTree->Branch( "centralityCut", &centralityCut);
  outTree->Branch( "yCut", &yCut);
  outTree->Branch( "mlbCut", &mlbCut);
  outTree->Branch( "m3Cut", &m3Cut);

  ID = 62219;
  nJetCut = 4;
  mtCut = 120;
  dphiCut = 0;
  hadChi2Cut = 1e+09;
  centralityCut = -1;
  metCut = 150;
  yCut = 0;
  mlbCut = 100;
  m3Cut = 250;
  mt2wCut = 0;
  outTree->Fill();

  ID = 78744;
  nJetCut = 5;
  mtCut = 120;
  dphiCut = 0;
  hadChi2Cut = 1e+09;
  centralityCut = 0.6;
  metCut = 150;
  yCut = 8;
  mlbCut = 1e+09;
  m3Cut = 250.;  
  mt2wCut = 0;
  outTree->Fill();

  ID = 51849;
  nJetCut = 3;
  mtCut = 120;
  dphiCut = 0.8;
  hadChi2Cut = 1e+09;
  centralityCut = -1;
  metCut = 150;
  yCut = 8;
  mlbCut = 1e+09;
  m3Cut = 250; 
  mt2wCut = 200;
  outTree->Fill();
  
  ID = 4;
  nJetCut = 4;
  mtCut = 120;
  dphiCut = 1.;
  hadChi2Cut = 5.;
  metCut = 150.;  
  mt2wCut = 200.;
  centralityCut = 0.7;
  yCut = 0;
  mlbCut = 1e+09;
  m3Cut = -2.;
  outTree->Fill();

  ID = 98872;
  nJetCut = 3;  
  mtCut = 120;
  dphiCut = 1;
  hadChi2Cut = 1e+09;
  centralityCut = 0.7;
  metCut = 250;
  yCut = 10;
  mlbCut = 1e+09;
  m3Cut = 0.;
  mt2wCut = 200;
  outTree->Fill();

  ID = 98868;
  nJetCut = 3;
  mtCut = 120;
  dphiCut = 1;
  hadChi2Cut = 1e+09;
  centralityCut = 0.7;
  metCut = 300;
  yCut = 10;
  mlbCut = 1e+09;
  m3Cut = 0.;
  mt2wCut = 200;
  outTree->Fill();

  ID = 68727;
  nJetCut = 4;
  mtCut = 120;  
  dphiCut = 1;
  hadChi2Cut = 1e+09;
  centralityCut = 0.7;
  metCut = 250;
  yCut = 10;
  mlbCut = 1e+09;
  m3Cut = 0;
  mt2wCut = 200;
  outTree->Fill();
  
  ID = 52812;
  nJetCut = 3;
  mtCut = 120;
  dphiCut = 0.8;
  centralityCut = -1;
  metCut = 150;
  yCut = 10;
  mlbCut = 1e+09;
  m3Cut = 0;
  mt2wCut = 200;
  outTree->Fill();

  ID = 116392;
  nJetCut = 4;
  mtCut = 120;
  dphiCut = 1;
  hadChi2Cut = 1e+09;
  centralityCut = 0.7;
  metCut = 300;
  yCut = 10;
  mlbCut = 1e+09;
  m3Cut = 0;
  mt2wCut = 200;
  outTree->Fill();
  
  ID = 10;
  nJetCut = 4;
  mtCut = 120;
  dphiCut = -1.;
  hadChi2Cut = 1e+09;
  centralityCut = -1.;
  metCut = 100;
  yCut = 0;
  mlbCut = 1e+09;
  m3Cut = -2.;
  mt2wCut = 0.;
  outTree->Fill();

  ID = 11;
  metCut = 150;
  outTree->Fill();

  ID = 12;
  metCut = 200;
  outTree->Fill();

  ID = 13;
  metCut = 250;
  outTree->Fill();

  ID = 14;
  metCut = 300;
  outTree->Fill();

  outTree->Write();

  return 0;
}
