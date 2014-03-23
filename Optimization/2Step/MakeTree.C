#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include <iostream>
#include <iomanip>
#include <vector>

bool pcp = true;

int MakeTree(){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/home/fcostanz/Bonsai/Optimization/";

  const int NSamples = 5;
  
  TString sample[NSamples];
  sample[0] = "T2tb";
  sample[1] = "DiLep";
  sample[2] = "OneLep";
  sample[3] = "WJets";
  sample[4] = "Rare";
  
  std::vector<Float_t> SRs;
  SRs.push_back(105089);
  SRs.push_back(108918);
  SRs.push_back(124470);
  SRs.push_back(1431);
  SRs.push_back(143910);
  SRs.push_back(1513);
  SRs.push_back(159438);
  SRs.push_back(16294);
  SRs.push_back(16297);
  SRs.push_back(163596);
  SRs.push_back(174960);
  SRs.push_back(1755);
  SRs.push_back(1784);
  SRs.push_back(20044);
  SRs.push_back(2027);
  SRs.push_back(21739);
  SRs.push_back(23689);
  SRs.push_back(2377);
  SRs.push_back(25678);
  SRs.push_back(259);
  SRs.push_back(2999);
  SRs.push_back(31357);
  SRs.push_back(31537);
  SRs.push_back(31620);
  SRs.push_back(31812);
  SRs.push_back(325);
  SRs.push_back(3942);
  SRs.push_back(4);
  SRs.push_back(47711);
  SRs.push_back(47736);
  SRs.push_back(47791);
  SRs.push_back(47898);
  SRs.push_back(47952);
  SRs.push_back(48006);
  SRs.push_back(48033);
  SRs.push_back(48114);
  SRs.push_back(48195);
  SRs.push_back(48199);
  SRs.push_back(48222);
  SRs.push_back(48276);
  SRs.push_back(48277);
  SRs.push_back(48280);
  SRs.push_back(48304);
  SRs.push_back(48791);
  SRs.push_back(49167);
  SRs.push_back(49248);
  SRs.push_back(50411);
  SRs.push_back(50492);
  SRs.push_back(51141);
  SRs.push_back(55084);
  SRs.push_back(55087);
  SRs.push_back(55110);
  SRs.push_back(56001);
  SRs.push_back(569);
  SRs.push_back(62211);
  SRs.push_back(62214);
  SRs.push_back(62299);
  SRs.push_back(62376);
  SRs.push_back(62470);
  SRs.push_back(62777);
  SRs.push_back(62874);
  SRs.push_back(63603);
  SRs.push_back(63697);
  SRs.push_back(63750);
  SRs.push_back(63828);
  SRs.push_back(64327);
  SRs.push_back(65043);
  SRs.push_back(66745);
  SRs.push_back(739);
  SRs.push_back(77766);
  SRs.push_back(77781);
  SRs.push_back(77784);
  SRs.push_back(77790);
  SRs.push_back(77815);
  SRs.push_back(77820);
  SRs.push_back(78004);
  SRs.push_back(78006);
  SRs.push_back(78009);
  SRs.push_back(78010);
  SRs.push_back(78180);
  SRs.push_back(78220);
  SRs.push_back(78327);
  SRs.push_back(78933);
  SRs.push_back(81657);
  SRs.push_back(82053);
  SRs.push_back(85572);
  SRs.push_back(893);

  TFile* cutFile = new TFile( "../1Step/"+sample[1]+"_optimization-1Step.root", "READ"); 
  TTree* cutTree;  
  cutFile->GetObject( "Optimization", cutTree);

  Float_t mtCut = 0.;
  Float_t nJetCut = 0.;
  Float_t topnessCut = 0.;
  Float_t mt2wCut = 0.;
  Float_t yCut = 0.;
  Float_t dphiCut = 0.;
  Float_t drlblCut = 0.;
  Float_t drlbgCut = 0.;
  Float_t chi2Cut = 0.;
  Float_t metCut = 0.;

  cutTree->SetBranchAddress( "mtCut", &mtCut);
  cutTree->SetBranchAddress( "nJetCut", &nJetCut);
  cutTree->SetBranchAddress( "topnessCut", &topnessCut);
  cutTree->SetBranchAddress( "mt2wCut", &mt2wCut);
  cutTree->SetBranchAddress( "yCut", &yCut);
  cutTree->SetBranchAddress( "dphiCut", &dphiCut);
  cutTree->SetBranchAddress( "drlblCut", &drlblCut);
  cutTree->SetBranchAddress( "drlbgCut", &drlbgCut);
  cutTree->SetBranchAddress( "chi2Cut", &chi2Cut);
  cutTree->SetBranchAddress( "metCut", &metCut);

  Float_t lumi=19500.;
  Float_t weight = 1.;
  /////////////////////////////////////////////////////
  //  Variable Definition
  ///////////////////////////////////////////////////// 

  Float_t globalWeight = 0.;
  Float_t FE = 0.;
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

  Float_t mStop = 0.;
  Float_t mLSP = 0.;

  Int_t pdgIdLep1 = 0;
  Char_t kinRegionFlag = false;
  Char_t searchRegionFlag = false;

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  
  TFile* inFile[NSamples];
  TTree* inTree[NSamples];
  for (int iSample = 0; iSample < NSamples; iSample++){
    TString inFileName = mainDir; inFileName += "/"; inFileName += sample[iSample]; inFileName +=".root";
    inFile[iSample] = new TFile(inFileName,"READ");

    inTree[iSample] = (TTree*)inFile[iSample]->Get("NoSystematic/bonsai");

    inTree[iSample]->SetBranchAddress( "GlobalWeight", &globalWeight);
    inTree[iSample]->SetBranchAddress( "FE", &FE);
    inTree[iSample]->SetBranchAddress( "TriggerWeight", &triggerWeight);
    inTree[iSample]->SetBranchAddress( "PUWeight", &puWeight);
    inTree[iSample]->SetBranchAddress( "isrWeight", &isrWeight);
    inTree[iSample]->SetBranchAddress( "topPtWeight", &topPtWeight);
    
    inTree[iSample]->SetBranchAddress( "LepFromTop", &lepFromTop);
    inTree[iSample]->SetBranchAddress( "Charginos", &charginos);
    
    inTree[iSample]->SetBranchAddress( "njets", &njets);
    inTree[iSample]->SetBranchAddress( "jet1", &jet1);
    inTree[iSample]->SetBranchAddress( "jet2", &jet2);
    inTree[iSample]->SetBranchAddress( "jet3", &jet3);
    inTree[iSample]->SetBranchAddress( "jet4", &jet4);
    
    inTree[iSample]->SetBranchAddress( "nbjets", &nbjets);
    inTree[iSample]->SetBranchAddress( "bjet1", &bjet1);
    inTree[iSample]->SetBranchAddress( "bjetHighestDisc", &bjetHighestDisc);
    
    inTree[iSample]->SetBranchAddress( "lPt", &lPt);
    inTree[iSample]->SetBranchAddress( "lEta", &lEta);
    inTree[iSample]->SetBranchAddress( "lRelIso", &lRelIso);
    
    inTree[iSample]->SetBranchAddress( "phiCorrMet", &phiCorrMet);
    
    inTree[iSample]->SetBranchAddress( "ht", &ht);
    inTree[iSample]->SetBranchAddress( "htRatio", &htRatio);
    inTree[iSample]->SetBranchAddress( "meff", &meff);
    inTree[iSample]->SetBranchAddress( "y", &y);
    
    inTree[iSample]->SetBranchAddress( "mt", &mt);
    inTree[iSample]->SetBranchAddress( "mlb1", &mlb1);
    inTree[iSample]->SetBranchAddress( "mlb", &mlb);
    inTree[iSample]->SetBranchAddress( "m3b", &m3b);
    inTree[iSample]->SetBranchAddress( "mt2w", &mt2w);
    inTree[iSample]->SetBranchAddress( "hadChi2", &hadChi2);
    inTree[iSample]->SetBranchAddress( "topness", &topness);
    
    inTree[iSample]->SetBranchAddress( "dphimin", &dphimin);
    inTree[iSample]->SetBranchAddress( "drlb1", &drlb1);
    inTree[iSample]->SetBranchAddress( "drlbmin", &drlbmin);

    inTree[iSample]->SetBranchAddress( "mStop", &mStop);
    inTree[iSample]->SetBranchAddress(  "mLSP",  &mLSP);

    inTree[iSample]->SetBranchAddress("pdgIdLep1",&pdgIdLep1);
    inTree[iSample]->SetBranchAddress("kinRegion",&kinRegionFlag);
    inTree[iSample]->SetBranchAddress("searchRegionPost",&searchRegionFlag);    
  }

  /////////////////////////////////////////////////////
  //  Output Definition
  ///////////////////////////////////////////////////// 
  TFile* outFile = new TFile( "Optimization-2Step.root", "RECREATE");
  outFile->cd();
  TTree* outTree = new TTree("Optimization", "Optimization subTree");
  
  Int_t SR = 0;

  TH2F *tt  = new TH2F( "tt", "tt",30,62.5,812.5, 15, -12.5, 362.5);
  TH2F *tb  = new TH2F( "tb", "tb",30,62.5,812.5, 15, -12.5, 362.5);
  TH2F *bb  = new TH2F( "bb", "bb",30,62.5,812.5, 15, -12.5, 362.5);
  
  Float_t bkg[NSamples - 1] = {};
  Float_t bkg2[NSamples - 1] = {};

  outTree->Branch( "SR", &SR, "SR/I");

  outTree->Branch(      "mtCut",      &mtCut,      "mtCut/F");
  outTree->Branch(    "njetCut",    &nJetCut,    "njetCut/F");
  outTree->Branch( "topnessCut", &topnessCut, "topnessCut/F");
  outTree->Branch(    "mt2wCut",    &mt2wCut,    "mt2wCut/F");
  outTree->Branch(       "yCut",       &yCut,       "yCut/F");
  outTree->Branch(    "dphiCut",    &dphiCut,    "dphiCut/F");
  outTree->Branch(   "drlblCut",   &drlblCut,   "drlblCut/F");
  outTree->Branch(   "drlbgCut",   &drlbgCut,   "drlbgCut/F");
  outTree->Branch(    "chi2Cut",    &chi2Cut,    "chi2Cut/F");
  outTree->Branch(     "metCut",     &metCut,     "metCut/F");

  outTree->Branch("tt","TH2F",&tt,32000,0);
  outTree->Branch("tb","TH2F",&tb,32000,0);
  outTree->Branch("bb","TH2F",&bb,32000,0);

  outTree->Branch( "diLep", &bkg[0],  "diLep/F");
  outTree->Branch("oneLep", &bkg[1], "oneLep/F");
  outTree->Branch( "wJets", &bkg[2],  "wJets/F");
  outTree->Branch(  "rare", &bkg[3],   "rare/F");

  outTree->Branch( "diLep2", &bkg2[0],  "diLep2/F");
  outTree->Branch("oneLep2", &bkg2[1], "oneLep2/F");
  outTree->Branch( "wJets2", &bkg2[2],  "wJets2/F");
  outTree->Branch(  "rare2", &bkg2[3],   "rare2/F");

  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 
  for( int isr = 0; isr < (int) SRs.size(); isr++){
    SR = SRs.at(isr);
    cutTree->GetEntry(SR);

    cout<<endl;
    cout<<"SR id "<<SR<<" "<<endl;       
    cout<<"mtCut = "<<mtCut<<endl;
    cout<<"nJetCut = "<<nJetCut<<endl;
    cout<<"topnessCut = "<<topnessCut<<endl;
    cout<<"mt2wCut = "<<mt2wCut<<endl;
    cout<<"yCut = "<<yCut<<endl;
    cout<<"dphiCut = "<<dphiCut<<endl;
    cout<<"drlblCut = "<<drlblCut<<endl;
    cout<<"drlbgCut = "<<drlbgCut<<endl;
    cout<<"chi2Cut = "<<chi2Cut<<endl;
    cout<<"metCut = "<<metCut<<endl;
    
    tt->Reset();
    tb->Reset();
    bb->Reset();    

    for (int iSample = 1; iSample < NSamples; iSample++){
      bkg[iSample-1] = 0.;
      bkg2[iSample-1] = 0.;
    }

    for (int iSample = 0; iSample < NSamples; iSample++){
      Int_t N = inTree[iSample]->GetEntries();

      /////////////////////////////////////////////////////
      //  Event Loop
      ///////////////////////////////////////////////////// 
      for (int ievt=0;ievt<N;++ievt){
	inTree[iSample]->GetEntry(ievt);
 
	weight = globalWeight * triggerWeight * puWeight * lumi;
	if (iSample < 0.001 ) {
	  if (charginos < 0.001 && FE > 0.999) continue;
	  weight *= isrWeight;
	}
      
	if (iSample < 2.001 && iSample > 0.999) weight *= topPtWeight;
      
	if (lRelIso > 0.1) continue;
	if (!searchRegionFlag) continue;
	
	if (mt < mtCut) continue;
	if (njets < nJetCut) continue;
	if (topness < topnessCut) continue;
	if (mt2w < mt2wCut) continue;
	if (y < yCut) continue;
	if (dphimin < dphiCut) continue; 
	if (drlb1 > drlblCut) continue; 
	if (drlb1 < drlbgCut) continue; 
	if (hadChi2 > chi2Cut) continue;
	if (phiCorrMet < metCut) continue;

	if (iSample < 0.001){
	  if (charginos == 0)
	    tt->Fill( mStop, mLSP, weight);
	  if (charginos == 1)
	    tb->Fill( mStop, mLSP, weight);
	  if (charginos == 2)
	    bb->Fill( mStop, mLSP, weight);
	}
	else{
	  bkg[iSample-1] += weight;
	  bkg2[iSample-1] += weight * weight;	  
	}
      }
    }
    outTree->Fill();
  }
  cout<<endl;

  outFile->cd();
  outTree->Write();
  outFile->Close();
  
  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();
  
  cutFile->Close();
  
  return 0;
}
