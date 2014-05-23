#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <vector>

bool pcp = true;

int MakeHistos(Int_t mySR){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/Systematics/";

  std::vector<TString> dirNames;
  dirNames.push_back(TString("NoSystematic"));
  dirNames.push_back(TString("JES_Up"));
  dirNames.push_back(TString("JES_Down"));
  dirNames.push_back(TString("BTagReweight_UpBC"));
  dirNames.push_back(TString("BTagReweight_DownBC"));
  dirNames.push_back(TString("BTagReweight_UpLight"));
  dirNames.push_back(TString("BTagReweight_DownLight"));

  TFile* srFile = new TFile( "../SearchRegions/SearchRegions.root", "READ"); 
  TTree* srTree;  
  srFile->GetObject( "SearchRegions", srTree);
  
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
  Float_t drlblCut = 10.;

  srTree->SetBranchAddress( "ID", &ID);
  srTree->SetBranchAddress( "nJetCut", &nJetCut);
  srTree->SetBranchAddress( "mtCut", &mtCut);
  srTree->SetBranchAddress( "dphiCut", &dphiCut);
  srTree->SetBranchAddress( "hadChi2Cut", &hadChi2Cut);
  srTree->SetBranchAddress( "metCut", &metCut);
  srTree->SetBranchAddress( "mt2wCut", &mt2wCut);
  srTree->SetBranchAddress( "centralityCut", &centralityCut);
  srTree->SetBranchAddress( "yCut", &yCut);
  srTree->SetBranchAddress( "mlbCut", &mlbCut);
  srTree->SetBranchAddress( "m3Cut", &m3Cut);

  Float_t lumi=19500.;
  Float_t weight = 1.;
  /////////////////////////////////////////////////////
  //  Variable Definition
  ///////////////////////////////////////////////////// 

  Float_t xs = 0.;
  Float_t globalWeight = 0.;
  Float_t FE = 0.;
  Float_t triggerWeight = 0.;
  Float_t puWeight = 0.;
  Float_t isrWeight = 0.;
  Float_t topPtWeight = 0.;
  Float_t T2ttL = 0.;
  Float_t T2ttR = 0.;
  
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
  Float_t m3 = 0.;
  Float_t centrality = 0.;
  Float_t mt2w = 0.;
  Float_t hadChi2 = 0.;
  Float_t topness = 0.;

  Float_t dphimin = 0.;
  Float_t drlb1 = 0.;
  Float_t drlbmin = 0.;

  Float_t mStop = 0.;
  Float_t mLSP = 0.;

  Int_t pdgIdLep1 = 0;
  Char_t kinRegionFlag = false;
  Char_t searchRegionPost = false;

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  
  TFile* inFile;
  TTree* inTree;
  TString inFileName = mainDir; inFileName += "/T2tb.root";
  inFile = new TFile(inFileName,"READ");

  TString outFileName = "./"; outFileName += mySR; outFileName += ".root";
  TFile* outFile = new TFile( outFileName, "RECREATE");

  TH2F *tt  = new TH2F(  "tt",  "tt",30,62.5,812.5, 15, -12.5, 362.5);
  TH2F *tb  = new TH2F(  "tb",  "tb",30,62.5,812.5, 15, -12.5, 362.5);
  TH2F *bb  = new TH2F(  "bb",  "bb",30,62.5,812.5, 15, -12.5, 362.5);
  
  TH2F *ttl = new TH2F( "ttl", "ttl",30,62.5,812.5, 15, -12.5, 362.5);
  TH2F *ttr = new TH2F( "ttr", "ttr",30,62.5,812.5, 15, -12.5, 362.5);
  TH2F *sig_tot = new TH2F( "sig_tot", "sig_tot",30,62.5,812.5, 15, -12.5, 362.5);
  

  for(int idir = 0; idir < (int) dirNames.size(); idir++){
    cout<<dirNames.at(idir)<<endl;
    inTree= (TTree*)inFile->Get(dirNames.at(idir) + TString("/bonsai"));

    inTree->SetBranchAddress( "xs", &xs);
    inTree->SetBranchAddress( "GlobalWeight", &globalWeight);
    inTree->SetBranchAddress( "FE", &FE);
    inTree->SetBranchAddress( "TriggerWeight", &triggerWeight);
    inTree->SetBranchAddress( "PUWeight", &puWeight);
    inTree->SetBranchAddress( "isrWeight", &isrWeight);
    inTree->SetBranchAddress( "topPtWeight", &topPtWeight);
    inTree->SetBranchAddress( "T2ttL", &T2ttL);
    inTree->SetBranchAddress( "T2ttR", &T2ttR);
  
    inTree->SetBranchAddress( "LepFromTop", &lepFromTop);
    inTree->SetBranchAddress( "Charginos", &charginos);
    
    inTree->SetBranchAddress( "njets", &njets);
    inTree->SetBranchAddress( "jet1", &jet1);
    inTree->SetBranchAddress( "jet2", &jet2);
    inTree->SetBranchAddress( "jet3", &jet3);
    inTree->SetBranchAddress( "jet4", &jet4);
  
    inTree->SetBranchAddress( "nbjets", &nbjets);
    inTree->SetBranchAddress( "bjet1", &bjet1);
    inTree->SetBranchAddress( "bjetHighestDisc", &bjetHighestDisc);
    
    inTree->SetBranchAddress( "lPt", &lPt);
    inTree->SetBranchAddress( "lEta", &lEta);
    inTree->SetBranchAddress( "lRelIso", &lRelIso);
    
    inTree->SetBranchAddress( "phiCorrMet", &phiCorrMet);
    
    inTree->SetBranchAddress( "ht", &ht);
    inTree->SetBranchAddress( "htRatio", &htRatio);
    inTree->SetBranchAddress( "meff", &meff);
    inTree->SetBranchAddress( "y", &y);
  
    inTree->SetBranchAddress( "mt", &mt);
    inTree->SetBranchAddress( "mlb1", &mlb1);
    inTree->SetBranchAddress( "mlb", &mlb);
    inTree->SetBranchAddress( "m3b", &m3b);
    inTree->SetBranchAddress( "m3", &m3);
    inTree->SetBranchAddress( "centrality", &centrality);
    inTree->SetBranchAddress( "mt2w", &mt2w);
    inTree->SetBranchAddress( "hadChi2", &hadChi2);
    inTree->SetBranchAddress( "topness", &topness);
    
    inTree->SetBranchAddress( "dphimin", &dphimin);
    inTree->SetBranchAddress( "drlb1", &drlb1);
    inTree->SetBranchAddress( "drlbmin", &drlbmin);
  
    inTree->SetBranchAddress( "mStop", &mStop);
    inTree->SetBranchAddress(  "mLSP",  &mLSP);
    
    inTree->SetBranchAddress("pdgIdLep1",&pdgIdLep1);
    inTree->SetBranchAddress("searchRegionPost",&searchRegionPost);    


    /////////////////////////////////////////////////////
    //  Output Definition
    /////////////////////////////////////////////////////
    
    outFile->mkdir(dirNames.at(idir));
    TDirectory* dir = outFile->GetDirectory(dirNames.at(idir));
    dir->cd(); 

    /////////////////////////////////////////////////////
    //  Event Loop
    ///////////////////////////////////////////////////// 

    srTree->GetEntry(mySR);
    mt2wCut -= 0.001;

    cout<<endl;
    cout<<"ID = "<<ID<<";"<<endl;    
    cout<<"nJetCut = "<<nJetCut<<";"<<endl;
    cout<<"mtCut = "<<mtCut<<";"<<endl;
    cout<<"dphiCut = "<<dphiCut<<";"<<endl;
    cout<<"hadChi2Cut = "<<hadChi2Cut<<";"<<endl;
    cout<<"metCut = "<<metCut<<";"<<endl;
    cout<<"mt2wCut = "<<mt2wCut<<";"<<endl;
    cout<<"centralityCut = "<<centralityCut<<";"<<endl;
    cout<<"yCut = "<<yCut<<";"<<endl;
    cout<<"mlbCut = "<<mlbCut<<";"<<endl;
    cout<<"m3Cut = "<<m3Cut<<";"<<endl;
    cout<<"mt2wCut = "<<mt2wCut<<";"<<endl;
    
    srTree->GetEntry(mySR);
    
    tt->Reset();
    ttl->Reset();
    ttr->Reset();
    
    tb->Reset();
    bb->Reset();    
    
    Int_t N = inTree->GetEntries(); cout<<"THERE ARE "<<N<<" EVENTS IN "<<inFileName<<endl;
    /////////////////////////////////////////////////////
    //  Event Loop
    ///////////////////////////////////////////////////// 
    for (int ievt=0;ievt<N;++ievt){
      inTree->GetEntry(ievt);
      
      weight = globalWeight * triggerWeight * puWeight * lumi * isrWeight;
      
      if (lRelIso > 0.1) continue;
      if (!searchRegionPost) continue;
      
      if (mt < mtCut) continue;
      if (njets < nJetCut) continue;
      if (mt2w < mt2wCut) continue;
      if (hadChi2 > hadChi2Cut) continue;
      if (y < yCut) continue;
      if (dphimin < dphiCut) continue; 
      if (phiCorrMet < metCut) continue;
      if (m3 < m3Cut) continue;
      if (centrality < centralityCut) continue;
      if (mlb > mlbCut) continue;
      
      if (charginos == 0){
	tt->Fill( mStop, mLSP, weight);
	ttl->Fill( mStop, mLSP, weight * T2ttL);
	ttr->Fill( mStop, mLSP, weight * T2ttR);
      }	    
      if (charginos == 1)
	tb->Fill( mStop, mLSP, weight);
      if (charginos == 2)
	bb->Fill( mStop, mLSP, weight);

      sig_tot->SetBinContent( sig_tot->FindBin(mStop, mLSP), xs * lumi);
    }
    dir->cd();
    
    tt->Write();
    ttl->Write();
    ttr->Write();
    
    tb->Write();
    bb->Write();
    sig_tot->Write();
  }

  outFile->Close(); 
  inFile->Close();
  srFile->Close();
  
  return 0;
}
