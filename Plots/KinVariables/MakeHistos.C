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

int MakeHistos( int iSample = 0){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/KinVariables/";

  const int NSamples = 8;  
  TString sample[NSamples];
  sample[0] = "DiLep";
  sample[1] = "OneLep";
  sample[2] = "WJets";
  sample[3] = "Rare";
  sample[4] = "T2tb-mStop175mLSP50";
  sample[5] = "T2tb-mStop200mLSP25";
  sample[6] = "T2tb-mStop325mLSP100";
  sample[7] = "T2tb-mStop550mLSP1";
  

  Float_t lumi=19500.;
  Float_t weight = 1.;

  /////////////////////////////////////////////////////
  //  Variable Definition
  ///////////////////////////////////////////////////// 

  Float_t globalWeight = 0.;
  
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

  TString inFileName = mainDir; inFileName += sample[iSample]; inFileName +=".root";
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

  TFile* outFile = new TFile( sample[iSample]+"_outfile.root", "RECREATE"); 
  TDirectory* baseDir = outFile->mkdir("KinVariabless");
  baseDir->cd();

  /////////////////////////////////////////////////////
  //  Histo Definition
  ///////////////////////////////////////////////////// 

  TH1D* lpth = new TH1D(  "lpt", "lep p_{T} [GeV]",  40, 0., 200.);
  TH1D* letah = new TH1D( "leta", "lep #Eta", 30, -3., 3.);
  TH1D* lRelIsoh = new TH1D( "lRelIso", "lep RelIso", 30, 0., 1.);

  TH1D* njetsh = new TH1D(  "njets", "jets multiplicity",    8, 1.5, 9.5);
  TH1D* jet1h  = new TH1D(  "jet1", "1st jet p_{T} [GeV]",   50, 0., 500.);
  TH1D* jet2h  = new TH1D(  "jet2", "2nd jet p_{T} [GeV]",   50, 0., 500.);
  TH1D* jet3h  = new TH1D(  "jet3", "3rd jet p_{T} [GeV]",   50, 0., 500.);
  TH1D* jet4h  = new TH1D(  "jet4", "4th jet p_{T} [GeV]",   50, 0., 500.);

  TH1D* nbjetsh    = new TH1D( "nbjets", "b jets multiplicity",  6, -0.5, 5.5);  
  TH1D* bjet1h     = new TH1D(     "bjet1", "Leading b jet p_{T} [GeV]", 25, 0., 500.);
  TH1D* bjetHighDh = new TH1D( "bjetHighD", "p_{T} of the highest b disc jet [GeV]",  25, 0., 500.);
  
  TH1D* meth = new TH1D( "MET", "MET [GeV]", 25, 0., 500.);
	
  TH1D* hth      = new TH1D(   "Ht",   "Ht [GeV]",  40, 0., 1000.);
  TH1D* htRatioh = new TH1D(   "HtRatio",   "HtRatio",  20, 0., 1.);
  TH1D* meffh    = new TH1D( "Meff", "Meff [GeV]",  80, 0., 2000.);
  TH1D* yh       = new TH1D(    "Y",    "Y [GeV^{1/2}]",  30, 0.,   30.);

  TH1D* mth  = new TH1D(  "Mt",  "Mt [GeV]", 25, 0., 500.);
  TH1D* mlb1h      = new TH1D(      "mlb1",     "Mlb1 [GeV]", 40, 0., 1000.);
  TH1D* mlbh  = new TH1D(  "mlb", "Mlb [GeV]", 40, 0., 1000.);
  TH1D* m3bh       = new TH1D(       "m3b",      "M3b [GeV]", 40, 0., 1000.);
  TH1D* mt2wh      = new TH1D(      "mt2w",     "MT2W [GeV]", 40, 0.,  500.);
  TH1D* hadChi2h = new TH1D(      "hadChi2",     "hadChi2 [GeV]", 20, 0., 10.);
  TH1D* topnessh = new TH1D(      "topness",     "topness [GeV]", 30, -15., 15.);

  TH1D* dphiminh = new TH1D( "dphimin", "min dPhi (MET, jet1/2)", 17, 0., TMath::Pi() * 17./16.);
  TH1D* drlb1h    = new TH1D(    "drlb1", "dR(lep, bjet1)", 20, 0., 5.);
  TH1D* drlbminh    = new TH1D(    "drlbmin", "min dR(lep, bjet)", 20, 0., 5.);

  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 

  float count=0;

  for (int ievt=0;ievt<N;++ievt){
    if (ievt%13454 == 0) {
      cout<<"Event number "<<ievt<<"\r"<<flush;
    }
    
    tree->GetEntry(ievt);

    weight = globalWeight * lumi;

    if (!kinRegion) continue;    
    //if (abs(pdgIdLep1) != 13) continue;

    //if (charginos != 0) continue;
    //if (lepFromTop != 1) continue;

    if ( mt < 100.) continue;
    if ( njets < 4.) continue;
    if ( mt2w < 200.)  continue;
    if ( y < 13.)  continue;
    if ( dphimin < 1.)  continue;
    if ( drlb1 < 2.) continue;
    if ( phiCorrMet < 200.) continue;
 
    count += weight;

    /////////////////////////////////////////////////////
    //  Histo Filling
    ///////////////////////////////////////////////////// 
    
    lpth->Fill( lPt, weight);
    letah->Fill( lEta, weight);
    lRelIsoh->Fill( lRelIso, weight);
    
    njetsh->Fill( njets, weight);
    jet1h->Fill( jet1, weight);
    jet2h->Fill( jet2, weight);
    jet3h->Fill( jet3, weight);
    jet4h->Fill( jet4, weight);
    
    nbjetsh->Fill( nbjets, weight);
    bjet1h->Fill( bjet1, weight);
    bjetHighDh->Fill( bjetHighestDisc, weight);
    
    meth->Fill( phiCorrMet, weight);
    
    hth->Fill( ht, weight);
    htRatioh->Fill( htRatio, weight);
    meffh->Fill( meff, weight);
    yh->Fill( y, weight);
    
    mth->Fill( mt, weight);
    mlb1h->Fill( mlb1, weight);
    mlbh->Fill( mlb, weight);
    m3bh->Fill( m3b, weight);
    mt2wh->Fill( mt2w, weight);
    hadChi2h->Fill( hadChi2, weight);
    topnessh->Fill( topness, weight);
    
    dphiminh->Fill( dphimin, weight);
    drlb1h->Fill( drlb1, weight);
    drlbminh->Fill( drlbmin, weight);
  }
  cout<<endl;
  cout<<count<<endl;

  /////////////////////////////////////////////////////
  //  WriteHisto
  /////////////////////////////////////////////////////
  
  baseDir->cd();

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

  outFile->Close();    
  inFile->Close();
}
