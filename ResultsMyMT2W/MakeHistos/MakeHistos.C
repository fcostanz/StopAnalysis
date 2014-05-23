#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>

bool pcp = true;
using namespace std;

int MakeHistos(int iSample = 0, int iSR = 6){
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/FullProduction/";

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
  
  srTree->GetEntry(iSR);
  mt2wCut -= 0.001;

  cout<<"ID = "<<ID<<";"<<endl;    
  cout<<"nJetCut = "<<nJetCut<<";"<<endl;
  cout<<"mtCut = "<<mtCut<<";"<<endl;
  cout<<"dphiCut = "<<dphiCut<<";"<<endl;
  cout<<"chi2Cut = "<<dphiCut<<";"<<endl;
  cout<<"metCut = "<<metCut<<";"<<endl;
  cout<<"mt2wCut = "<<mt2wCut<<";"<<endl;
  cout<<"centralityCut = "<<centralityCut<<";"<<endl;
  cout<<"yCut = "<<yCut<<";"<<endl;
  cout<<"mlbCut = "<<mlbCut<<";"<<endl;
  cout<<"m3Cut = "<<m3Cut<<";"<<endl;
  cout<<"mt2wCut = "<<mt2wCut<<";"<<endl;

  const int NSamples = 21;  
  const int NSignals = 4;
  const int NControlRegions = 6;
  const int NLeps = 3;
  const int NDirs = NControlRegions * NLeps;

  TString sample[NSamples];
  sample[0] = "SingleMu";
  sample[1] = "SingleElectron";
  sample[2] = "DiLep";
  sample[3] = "OneLep";
  sample[4] = "WJets";
  sample[5] = "Rare";
  sample[6] = "SemiLep";
  sample[7] = "SingleTop";
  /*sample[7] = "T2tb-mStop175mLSP50";
  sample[8] = "T2tb-mStop200mLSP25";
  sample[9] = "T2tb-mStop325mLSP100";
  sample[10] = "T2tb-mStop550mLSP1";*/
  sample[8] = "DiLepMG";
  sample[8] = "OneLepMG";  
  sample[8] = "DiLepMCatNLO";
  sample[8] = "OneLepMCatNLO";
  sample[9] = "DiLepScaleUp";
  sample[10] = "OneLepScaleUp";
  sample[11] = "DiLepScaleDown";
  sample[12] = "OneLepScaleDown";
  sample[13] = "DiLepMatchUp";
  sample[14] = "OneLepMatchUp";
  sample[15] = "DiLepMatchDown";
  sample[16] = "OneLepMatchDown";
  sample[17] = "DiLepMass166p5";
  sample[18] = "OneLepMass166p5";
  sample[19] = "DiLepMass178p5";
  sample[20] = "OneLepMass178p5";

  double lumi=19500.;
  double weight = 1.;

  bool lepFlag;
  TString lep[NLeps];
  lep[0] = "El";
  lep[1] = "Mu";
  lep[2] = "ElAndMu";

  TString controlRegion[NControlRegions];
  controlRegion[0] = "Preselection";
  controlRegion[1] = "SearchRegionPreIsoTrackVeto";
  controlRegion[2] = "SearchRegionPostIsoTrackVeto";
  controlRegion[3] = "CR1";
  controlRegion[4] = "CR4";
  controlRegion[5] = "CR5";

  bool flag[NDirs];
  TString controlDirName[NDirs];
  TDirectory* controlDir[NDirs];

  for (int ilep = 0; ilep < NLeps; ilep++){
    for ( int iControlRegion = 0; iControlRegion < NControlRegions; iControlRegion++){
      controlDirName[ilep * NControlRegions + iControlRegion ] = lep[ilep];
      controlDirName[ilep * NControlRegions + iControlRegion ] += "-";
      controlDirName[ilep * NControlRegions + iControlRegion ] += controlRegion[iControlRegion];
    }
  }

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  
  std::cout<<"Running over Sample "<<sample[iSample]<<std::endl;

  TString inFileName = mainDir; inFileName += sample[iSample]; inFileName +=".root";
  TFile* inFile = new TFile(inFileName,"READ");
  if (!inFile->IsOpen()){
    std::cout<<"not open"<<std::endl;
  }

  /////////////////////////////////////////////////////
  //  Tree Definition
  /////////////////////////////////////////////////////  
  TTree* tree;
  tree= (TTree*)inFile->Get("NoSystematic/bonsai");
  int N = tree->GetEntries();  cout<<"THERE ARE "<<N<<" EVENTS IN "<<inFileName<<endl;
  
  Float_t globalWeight = 0.;
  Float_t triggerWeight = 0.;
  Float_t npv = 0.;
  Float_t ngoodpv = 0.;
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
  Float_t bdiscH = 0.;

  Float_t lPt = 0.;
  Float_t lEta = 0.;
  Float_t lRelIso = 0.;

  Float_t isoTrack = 0.;
  Float_t tauVeto = 0.;
  
  Float_t rawmet = 0.;
  Float_t typeImet = 0.;
  Float_t phiCorrMet = 0.;
  
  Float_t ht = 0.;
  Float_t ht3 = 0.;
  Float_t ht4 = 0.;
  Float_t ht5 = 0.;
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

  Int_t pdgIdLep1 = 0;
  Int_t pdgIdLep2 = 0;

  Char_t kinRegion = false;
  Char_t searchRegionPre = false;
  Char_t searchRegionPost = false;
  Char_t CR1 = false;
  Char_t CR4 = false;
  Char_t CR5 = false;

  tree->SetBranchAddress( "GlobalWeight", &globalWeight);
  tree->SetBranchAddress( "TriggerWeight", &triggerWeight);
  tree->SetBranchAddress( "NPV", &npv);
  tree->SetBranchAddress( "NgoodPV", &ngoodpv);
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
  tree->SetBranchAddress( "discH", &bdiscH);

  tree->SetBranchAddress( "lPt", &lPt);
  tree->SetBranchAddress( "lEta", &lEta);
  tree->SetBranchAddress( "lRelIso", &lRelIso);

  tree->SetBranchAddress( "phiCorrMet", &phiCorrMet);
  
  tree->SetBranchAddress( "ht", &ht);
  tree->SetBranchAddress( "ht3", &ht3);
  tree->SetBranchAddress( "ht4", &ht4);
  tree->SetBranchAddress( "ht5", &ht5);
  tree->SetBranchAddress( "htRatio", &htRatio);
  tree->SetBranchAddress( "meff", &meff);
  tree->SetBranchAddress( "y", &y);

  tree->SetBranchAddress( "mt", &mt);
  tree->SetBranchAddress( "mlb1", &mlb1);
  tree->SetBranchAddress( "mlb", &mlb);
  tree->SetBranchAddress( "m3b", &m3b);
  tree->SetBranchAddress( "m3", &m3);
  tree->SetBranchAddress( "centrality", &centrality);
  tree->SetBranchAddress( "mt2w", &mt2w);
  tree->SetBranchAddress( "hadChi2", &hadChi2);
  tree->SetBranchAddress( "topness", &topness);

  tree->SetBranchAddress( "dphimin", &dphimin);
  tree->SetBranchAddress( "drlb1", &drlb1);
  tree->SetBranchAddress( "drlbmin", &drlbmin);

  tree->SetBranchAddress("pdgIdLep1",&pdgIdLep1);
  tree->SetBranchAddress("pdgIdLep2",&pdgIdLep2);

  tree->SetBranchAddress("kinRegion",&kinRegion);
  tree->SetBranchAddress("searchRegionPre",&searchRegionPre);
  tree->SetBranchAddress("searchRegionPost",&searchRegionPost);
  tree->SetBranchAddress("CR1",&CR1);
  tree->SetBranchAddress("CR4",&CR4);
  tree->SetBranchAddress("CR5",&CR5);

  /////////////////////////////////////////////////////
  //  Output Definition
  /////////////////////////////////////////////////////

  //Branching Ratio
 
  TString outFileName = "./"; outFileName += iSR; outFileName += "/"; outFileName += sample[iSample]; outFileName += ".root";
  TFile* outFile = new TFile( outFileName, "RECREATE");
  outFile->cd();

  TH1D* wTopPtReweighth[NDirs];
  TH1D* woTopPtReweighth[NDirs];

  TH1D* npvh[NDirs];
  TH1D* ngoodpvh[NDirs];

  TH1D* lpth[NDirs];
  TH1D* letah[NDirs];
  TH1D* lrelisoh[NDirs];

  TH1D* njetsh[NDirs];
  TH1D* jet1h[NDirs];
  TH1D* jet2h[NDirs];
  TH1D* jet3h[NDirs];
  TH1D* jet4h[NDirs];

  TH1D* nbjetsh[NDirs];
  TH1D* bjet1h[NDirs];
  TH1D* bjetHighDh[NDirs];
  TH1D* bdiscHh[NDirs];
  
  TH1D* hth[NDirs];
  TH1D* ht3h[NDirs];
  TH1D* ht4h[NDirs];
  TH1D* ht5h[NDirs];
  TH1D* htratioh[NDirs];
  TH1D* meth[NDirs];
  TH1D* meffh[NDirs];
  TH1D* yh[NDirs];

  TH1D* mth[NDirs];
  TH1D* mt1h[NDirs];
  TH1D* mt2h[NDirs];
  TH1D* mt3h[NDirs];
  TH1D* mt4h[NDirs];
  TH1D* mt5h[NDirs];
  TH1D* mlb1h[NDirs];
  TH1D* mlbh[NDirs];
  TH1D* m3bh[NDirs];
  TH1D* m3h[NDirs];
  TH1D* centralityh[NDirs];
  TH1D* mt2wh[NDirs];
  TH1D* hadchi2h[NDirs];
  TH1D* topnessh[NDirs];

  TH1D* dphiminh[NDirs];
  TH1D* drlb1h[NDirs];
  TH1D* drlbminh[NDirs];

  TString dirName = "";
  //if (outFile->GetDirectory(dirName)) outFile->Delete(dirName + ";*");
  outFile->mkdir( dirName);
  TDirectory* SRDir = outFile->GetDirectory( dirName);
  SRDir->cd();
  
  for ( int iDir = 0; iDir < NDirs; iDir++){    
    SRDir->mkdir( controlDirName[iDir]);
    controlDir[iDir] = SRDir->GetDirectory( controlDirName[iDir]);
    controlDir[iDir]->cd();

    wTopPtReweighth[iDir] = new TH1D( "wTopPtReweight", "wTopPtReweight", 1, 0., 1.);
    woTopPtReweighth[iDir] = new TH1D( "woTopPtReweight", "woTopPtReweight", 1, 0., 1.);

    npvh[iDir] = new TH1D( "npv", "NPV", 51, -0.5, 50.5);
    ngoodpvh[iDir] = new TH1D( "ngoodpv", "NgoodPV", 51, -0.5, 50.5);

    lpth[iDir]     = new TH1D( "lpt", "lep p_{T} [GeV]",  12, 25., 500.);
    letah[iDir]    = new TH1D( "leta", "lep #Eta", 30, -3., 3.);
    lrelisoh[iDir] = new TH1D( "lRelIso", "lep RelIso", 30, 0., 1.);
    
    njetsh[iDir] = new TH1D(  "njets", "jets multiplicity",    10, -0.5, 9.5);
    jet1h[iDir]  = new TH1D(  "jet1", "1st jet p_{T} [GeV]",   25, 0., 500.);
    jet2h[iDir]  = new TH1D(  "jet2", "2nd jet p_{T} [GeV]",   25, 0., 500.);
    jet3h[iDir]  = new TH1D(  "jet3", "3rd jet p_{T} [GeV]",   25, 0., 500.);
    jet4h[iDir]  = new TH1D(  "jet4", "4th jet p_{T} [GeV]",   25, 0., 500.);
    
    nbjetsh[iDir]    = new TH1D( "nbjets", "b jets multiplicity",  6, -0.5, 5.5);  
    bjet1h[iDir]     = new TH1D(     "bjet1", "Leading b jet p_{T} [GeV]", 25, 0., 500.);
    bjetHighDh[iDir] = new TH1D( "bjetHighD", "p_{T} of the highest b disc jet [GeV]",  12, 0., 500.);
    bdiscHh[iDir]    = new TH1D( "bdisc", "bdisc",  20, 0., 1.);  

    meth[iDir] = new TH1D( "MET", "MET [GeV]", 12, 0., 400.);
    
    hth[iDir]      = new TH1D(   "Ht",   "Ht [GeV]",  20, 0., 1000.);
    ht3h[iDir]     = new TH1D(   "Ht3",   "Ht3 [GeV]",  20, 0., 1000.);
    ht4h[iDir]     = new TH1D(   "Ht4",   "Ht4 [GeV]",  20, 0., 1000.);
    ht5h[iDir]     = new TH1D(   "Ht5",   "Ht5 [GeV]",  20, 0., 1000.);
    htratioh[iDir] = new TH1D(   "HtRatio",   "HtRatio",  20, 0., 1.);
    meffh[iDir]    = new TH1D( "Meff", "Meff [GeV]",  40, 0., 1000.);
    yh[iDir]       = new TH1D(    "Y",    "Y [GeV^{1/2}]",  15, 0.,   30.);
    
    mth[iDir]      = new TH1D(  "Mt",  "Mt [GeV]", 30, 0., 300.);
    mt1h[iDir]     = new TH1D( "Mt1Jet",  "Mt1 [GeV]", 30, 0., 300.);
    mt2h[iDir]     = new TH1D( "Mt2Jet",  "Mt2 [GeV]", 30, 0., 300.);
    mt3h[iDir]     = new TH1D( "Mt3Jet",  "Mt3 [GeV]", 30, 0., 300.);
    mt4h[iDir]     = new TH1D( "Mt4Jet",  "Mt4 [GeV]", 30, 0., 300.);
    mt5h[iDir]     = new TH1D( "Mt5+Jet",  "Mt5 [GeV]", 30, 0., 300.);
    mlb1h[iDir]    = new TH1D(      "mlb1",     "Mlb1 [GeV]", 10, 0., 500.);
    mlbh[iDir]     = new TH1D(  "mlb", "Mlb [GeV]", 10, 0., 500.);
    m3bh[iDir]     = new TH1D(       "m3b",      "M3b [GeV]", 15, 0., 500.);
    m3h[iDir]      = new TH1D(        "m3",       "M3 [GeV]", 15, 0., 500.);
    centralityh[iDir] = new TH1D( "centrality", "centrality", 10, 0., 1.);
    mt2wh[iDir]    = new TH1D(      "mt2w",     "MT2W [GeV]", 15, 0.,  500.);
    hadchi2h[iDir] = new TH1D(      "hadChi2",     "hadChi2 [GeV]", 20, 0., 10.);
    topnessh[iDir] = new TH1D(      "topness",     "topness [GeV]", 30, -15., 15.);
    
    dphiminh[iDir] = new TH1D( "dphimin", "min dPhi (MET, jet1/2)", 17, 0., TMath::Pi() * 17./16.);
    drlb1h[iDir]   = new TH1D(    "drlb1", "dR(lep, bjet1)", 20, 0., 5.);
    drlbminh[iDir] = new TH1D(    "drlbmin", "min dR(lep, bjet)", 20, 0., 5.);
  }
  outFile->cd();

  for (int ievt=0;ievt<N;++ievt){    
    tree->GetEntry(ievt);
    
    if (ievt%134563 == 0) cout<<"Event number "<<ievt<<"\r"<<flush;    

    if (iSample == 0 || iSample == 1) weight = 1.;           
    else weight = globalWeight * triggerWeight * puWeight * topPtWeight * lumi;
    //else weight = globalWeight * triggerWeight * puWeight * lumi;

    if (lRelIso > 0.1) continue;
    
    for (int ilep = 0; ilep < NLeps; ilep++){
      lepFlag = true;
      if (ilep == 0) lepFlag = (abs(pdgIdLep1) == 11);
      if (ilep == 1) lepFlag = (abs(pdgIdLep1) == 13);

      flag[ilep * NControlRegions + 0] = (searchRegionPost || CR1) && lepFlag;
      flag[ilep * NControlRegions + 1] = searchRegionPre && lepFlag;
      flag[ilep * NControlRegions + 2] = searchRegionPost && lepFlag;
      flag[ilep * NControlRegions + 3] = CR1 && lepFlag;
      flag[ilep * NControlRegions + 4] = CR4 && lepFlag;
      flag[ilep * NControlRegions + 5] = CR5 && lepFlag;
    }
    
    for ( int iDir = 0; iDir < NDirs; iDir++){
      if (!flag[iDir])
	continue;
      
      /////////////////////////////////////////////////////
      //  Histo Filling
      /////////////////////////////////////////////////////
      

      bool allCuts    = (njets - nJetCut) > -0.0001 && mt2w > mt2wCut && y > yCut && dphimin > dphiCut && drlb1 < drlblCut && phiCorrMet > metCut && m3 > m3Cut;
      bool allButNJet =                                mt2w > mt2wCut && y > yCut && dphimin > dphiCut && drlb1 < drlblCut && phiCorrMet > metCut && m3 > m3Cut;
      bool allButMt2w = (njets - nJetCut) > -0.0001 &&                   y > yCut && dphimin > dphiCut && drlb1 < drlblCut && phiCorrMet > metCut && m3 > m3Cut;
      bool allButY    = (njets - nJetCut) > -0.0001 && mt2w > mt2wCut &&             dphimin > dphiCut && drlb1 < drlblCut && phiCorrMet > metCut && m3 > m3Cut;
      bool allButDphi = (njets - nJetCut) > -0.0001 && mt2w > mt2wCut && y > yCut &&                      drlb1 < drlblCut && phiCorrMet > metCut && m3 > m3Cut;
      bool allButDrlb = (njets - nJetCut) > -0.0001 && mt2w > mt2wCut && y > yCut && dphimin > dphiCut &&                     phiCorrMet > metCut && m3 > m3Cut;
      bool allButMet  = (njets - nJetCut) > -0.0001 && mt2w > mt2wCut && y > yCut && dphimin > dphiCut && drlb1 < drlblCut &&                        m3 > m3Cut;
      bool allButM3   = (njets - nJetCut) > -0.0001 && mt2w > mt2wCut && y > yCut && dphimin > dphiCut && drlb1 < drlblCut && phiCorrMet > metCut              ;

      bool allButCentrality = allCuts && mlb1 < mlbCut                               && hadChi2 < hadChi2Cut;
      bool allButMlb        = allCuts &&                  centrality > centralityCut && hadChi2 < hadChi2Cut; 
      bool allButHadChi2    = allCuts && mlb1 < mlbCut && centrality > centralityCut                        ;

      allCuts    &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut;
      allButNJet &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut;
      allButMt2w &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut;
      allButY    &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut; 
      allButDphi &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut; 
      allButDrlb &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut; 
      allButMet  &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut;
      allButM3   &= centrality > centralityCut && mlb1 < mlbCut && hadChi2 < hadChi2Cut; 
    
      //if (allButMt2w) mt2wh[iDir]->Fill( mt2w, weight);

      if (allButNJet) {
	njetsh[iDir]->Fill( njets, weight);
	
	if      (fabs(njets - 1) < .0001) mt1h[iDir]->Fill( mt, weight);
	else if (fabs(njets - 2) < .0001) mt2h[iDir]->Fill( mt, weight);
	else if (fabs(njets - 3) < .0001) mt3h[iDir]->Fill( mt, weight);	
	else if (fabs(njets - 4) < .0001) mt4h[iDir]->Fill( mt, weight);
	else mt5h[iDir]->Fill( mt, weight);
      }
      /*if (allButY) yh[iDir]->Fill( y, weight);
      if (allButDphi) dphiminh[iDir]->Fill( dphimin, weight);
      if (allButDrlb) drlb1h[iDir]->Fill( drlb1, weight);
      if (allButMet) meth[iDir]->Fill( phiCorrMet, weight);
      if (allButM3) m3h[iDir]->Fill( m3, weight);
      if (allButCentrality) centralityh[iDir]->Fill( centrality, weight);
      if (allButMlb) mlb1h[iDir]->Fill( mlb1, weight);
      if (allButHadChi2) hadchi2h[iDir]->Fill( hadChi2, weight);*/

      if (allCuts) {
	wTopPtReweighth[iDir]->Fill(0., weight);
	woTopPtReweighth[iDir]->Fill(0., weight / topPtWeight);

	npvh[iDir]->Fill( npv, weight);
	ngoodpvh[iDir]->Fill( ngoodpv, weight);

	lpth[iDir]->Fill( lPt, weight);
	letah[iDir]->Fill( lEta, weight);
	lrelisoh[iDir]->Fill( lRelIso, weight);
	
	jet1h[iDir]->Fill( jet1, weight);
	jet2h[iDir]->Fill( jet2, weight);
	jet3h[iDir]->Fill( jet3, weight);
	jet4h[iDir]->Fill( jet4, weight);
	
	nbjetsh[iDir]->Fill( nbjets, weight);
	bjet1h[iDir]->Fill( bjet1, weight);
	bjetHighDh[iDir]->Fill( bjetHighestDisc, weight);
	bdiscHh[iDir]->Fill( bdiscH, weight);
	
	hth[iDir]->Fill( ht, weight);
	ht3h[iDir]->Fill( ht3, weight);
	ht4h[iDir]->Fill( ht4, weight);
	ht5h[iDir]->Fill( ht5, weight);
	htratioh[iDir]->Fill( htRatio, weight);
	meffh[iDir]->Fill( meff, weight);
      
	mth[iDir]->Fill( mt, weight);
	mlb1h[iDir]->Fill( mlb1, weight);
	m3bh[iDir]->Fill( m3b, weight);

	topnessh[iDir]->Fill( topness, weight);
      
	drlbminh[iDir]->Fill( drlbmin, weight);

	mt2wh[iDir]->Fill( mt2w, weight);

	yh[iDir]->Fill( y, weight);
	dphiminh[iDir]->Fill( dphimin, weight);
	drlb1h[iDir]->Fill( drlb1, weight);
	meth[iDir]->Fill( phiCorrMet, weight);
	m3h[iDir]->Fill( m3, weight);
	centralityh[iDir]->Fill( centrality, weight);
	mlb1h[iDir]->Fill( mlb1, weight);
	hadchi2h[iDir]->Fill( hadChi2, weight);	
      }
    }
  }
  
  outFile->Write();
  outFile->Close();
  
  inFile->Close();
  
  return 0;
}
