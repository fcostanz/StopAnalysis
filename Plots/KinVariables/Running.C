#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<double> >          LorentzED;
using namespace ROOT::Math::VectorUtil;

void IntPlot (TH1D* plot, TH1D* intPlot, float weight=1.);

float gettrigweight(int id1, float pt, float eta);

bool pcp = true;

int Running(){
  
  gSystem->Load("libRA4bDict.so");

  Jet::SetWP("8TeV"); 

  TH1::SetDefaultSumw2(true);
  
  if(pcp)cout<<"going to set inputs"<<endl;

  const int NSamples = 6;
  const int NMETcuts = 6;
  const int NSigRegions = 4;
  
  TString sample[NSamples];
  sample[0] = "Data";
  sample[1] = "WJets";
  sample[2] = "TTJets";
  sample[3] = "QCD";
  sample[4] = "SingleTop";
  sample[5] = "Rare";
  
  double lumi=19500;
  double weight = 1.;

  bool METFlag[NMETcuts];
  double METcut[NMETcuts] = { 0., 50., 60., 80., 100., 150.};

  bool sigRegionFlag[NSigRegions];
  TString sigRegion[NSigRegions];
  sigRegion[0] = "Mu";
  sigRegion[1] = "Ele1.5";
  sigRegion[2] = "Ele2.5";
  sigRegion[3] = "MuOrEle2.5";

  TFile* outFile[NSamples];
  TDirectory* baseDir[NSamples];
  TDirectory* metDir[NSamples][NMETcuts][NSigRegions];

  TH1D* lpth[NSamples][NMETcuts][NSigRegions];
  TH1D* letah[NSamples][NMETcuts][NSigRegions];

  TH1D* njetsh[NSamples][NMETcuts][NSigRegions];
  TH1D* nbjetsh[NSamples][NMETcuts][NSigRegions];

  TH1D* jet1h[NSamples][NMETcuts][NSigRegions];
  TH1D* bjet1h[NSamples][NMETcuts][NSigRegions];
  TH1D* bjetHighDh[NSamples][NMETcuts][NSigRegions];
  
  TH1D* meth[NSamples][NMETcuts][NSigRegions];
  TH1D* mth[NSamples][NMETcuts][NSigRegions];

  TH1D* hth[NSamples][NMETcuts][NSigRegions];
  TH1D* meffh[NSamples][NMETcuts][NSigRegions];
  TH1D* yh[NSamples][NMETcuts][NSigRegions];

  TH1D* mlb1h[NSamples][NMETcuts][NSigRegions];
  TH1D* mlb_hemih[NSamples][NMETcuts][NSigRegions];
  TH1D* m3bh[NSamples][NMETcuts][NSigRegions];
  TH1D* mt2wh[NSamples][NMETcuts][NSigRegions];

  TH1D* dphiminh[NSamples][NMETcuts][NSigRegions];
  TH1D* drlbh[NSamples][NMETcuts][NSigRegions];
  
  for ( int iSample = 0; iSample < NSamples; iSample++){
    outFile[iSample] = new TFile( sample[iSample]+"/NoSub/"+sample[iSample]+"_NoSub_STOP_NoTail.root", "RECREATE");
    outFile[iSample]->cd();
    baseDir[iSample] = outFile[iSample]->mkdir("ControlPlots");
    baseDir[iSample]->cd();
    for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
      for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
	TString dirName = "MET"; dirName += METcut[iMETcut]; dirName += "_"; dirName += sigRegion[iSigRegion];
	metDir[iSample][iMETcut][iSigRegion] = baseDir[iSample]->mkdir( dirName);
	metDir[iSample][iMETcut][iSigRegion]->cd();


	lpth[iSample][iMETcut][iSigRegion] = new TH1D(  "lpt", "Lepton's pt",  50, 0., 500.);
	letah[iSample][iMETcut][iSigRegion] = new TH1D( "leta", "Lepton's eta", 30, -3., 3.);
	
	njetsh[iSample][iMETcut][iSigRegion]  = new TH1D(  "njets", "jets multiplicity",    14, 0., 14.);
	nbjetsh[iSample][iMETcut][iSigRegion] = new TH1D( "nbjets", "b jets multiplicity",  10, 0., 10.);
	
	jet1h[iSample][iMETcut][iSigRegion]      = new TH1D(      "jet1", "Leading jet pt",   50, 0., 500.);
	bjet1h[iSample][iMETcut][iSigRegion]     = new TH1D(     "bjet1", "Leading b jet pt", 50, 0., 500.);
	bjetHighDh[iSample][iMETcut][iSigRegion] = new TH1D( "bjetHighD", "pt of the highest b disc jet",  50, 0., 500.);
	 
	meth[iSample][iMETcut][iSigRegion] = new TH1D( "MET", "MET", 50, 0., 500.);
	mth[iSample][iMETcut][iSigRegion]  = new TH1D(  "Mt",  "Mt", 50, 0., 500.);
	
	hth[iSample][iMETcut][iSigRegion]   = new TH1D(   "Ht",   "Ht",  100, 0., 1000.);
	meffh[iSample][iMETcut][iSigRegion] = new TH1D( "Meff", "Meff", 2000, 0., 2000.);
	yh[iSample][iMETcut][iSigRegion]    = new TH1D(    "Y",    "Y",   40, 0.,   20.);
	
	mlb1h[iSample][iMETcut][iSigRegion]      = new TH1D(      "mlb1",     "Mlb1", 100, 0., 1000.);
	mlb_hemih[iSample][iMETcut][iSigRegion]  = new TH1D(  "mlb_hemi", "Mlb_hemi", 100, 0., 1000.);
	m3bh[iSample][iMETcut][iSigRegion]       = new TH1D(       "m3b",      "M3b", 100, 0., 1000.);
	mt2wh[iSample][iMETcut][iSigRegion]      = new TH1D(      "mt2w",     "MT2W", 100, 0., 1000.);
	
	dphiminh[iSample][iMETcut][iSigRegion] = new TH1D( "dphimin", "min dPhi (MET, jet1/2)", 31, 0., TMath::Pi() * 31./30.);
	drlbh[iSample][iMETcut][iSigRegion]    = new TH1D(    "drlb", "dR lep bjet", 20, 0., 10.);
	 
      }
    }
  }

  for ( int iSample = 0; iSample < NSamples; iSample++) {

    std::cout<<"Running over Sample "<<sample[iSample]<<std::endl;
    TString inFileName = "tree/"; inFileName += sample[iSample]; inFileName +=".root";
    
    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile->IsOpen()){
      std::cout<<"not open"<<std::endl;
    }
    TTree* tree;
    tree= (TTree*)inFile->Get("/skimmingTree/prova2/skimmingTree");
    //===========================================
    if(pcp)cout<<"inputs set!"<<endl;
    int N = tree->GetEntries();  cout<<"THERE ARE "<<N<<" EVENTS IN "<<inFileName<<endl;
    
    EventInfo*              info;
    std::vector< Particle>* tracks;
    std::vector< Muon>*     muons;
    std::vector< Electron>* electrons;
    std::vector< Tau>*      taus;
    std::vector< Jet>*      jets;
    std::vector< GenJet>*   genJets;
    LorentzV*               unclusteredEnergy;
    LorentzV*               rawMET;
    LorentzV*               typeIMET;
    LorentzV*               typeIPhiCorrMET;
    LorentzV*               mvaMET;
    LorentzV*               caloMET;
    
    tree->SetBranchAddress(            "Info.",            &info);
    tree->SetBranchAddress(           "Tracks",          &tracks);
    tree->SetBranchAddress(            "Muons",           &muons);
    tree->SetBranchAddress(        "Electrons",       &electrons);
    tree->SetBranchAddress(             "Taus",            &taus);
    tree->SetBranchAddress(             "Jets",            &jets);
    tree->SetBranchAddress(          "GenJets",         &genJets);
    tree->SetBranchAddress(          "RawMET.",          &rawMET);
    tree->SetBranchAddress(        "TypeIMET.",        &typeIMET);
    tree->SetBranchAddress( "TypeIPhiCorrMET.", &typeIPhiCorrMET);
    tree->SetBranchAddress(          "MvaMET.",          &mvaMET);
    tree->SetBranchAddress(         "CaloMET.",         &caloMET); 
    
    int count = 0;

    Event event;
    for (int ievt=0;ievt<N;++ievt){
      
      if (ievt%13454 == 0) {
	cout<<"Event number "<<ievt<<"\r"<<flush;
      }
      
      tree->GetEntry(ievt);

      event.SetInfo( *info);
      event.SetTracks( *tracks);
      event.SetMuons( *muons);
      event.SetElectrons( *electrons);
      event.SetTaus( *taus);
      event.SetJets( *jets);
      event.SetGenJets( *genJets);
      event.SetCaloMET( *caloMET);
      event.SetTypeIPhiCorrMET( *typeIPhiCorrMET);

      /*
      cout<<"HT = "<<event.HT();
      cout<<"; HTratio = "<<event.HTratio();
      cout<<"; Meff = "<<event.Meff();
      cout<<"; Y = "<<event.Y();
      cout<<"; MT = "<<event.MT();
      cout<<"; nbjets = "<<event.nBJets();     
      cout<<"; Mlb1 = "<<event.Mlb1();
      cout<<"; Mlb_hemi = "<<event.Mlb_hemi();
      cout<<"; M3b = "<<event.M3b();
      cout<<"; MT2W = "<<event.MT2W();
      cout<<"; dp = "<<event.DeltaPhiMinj12m();
      cout<<"; dr = "<<event.DeltaRlb1();
      cout<<endl;
      */

      if (iSample == 0){
	weight = 1.;
	if ( event.Info()->badLaserFilter) continue;
	if ( event.Info()->badXTalLaserCorrectionFilter) continue;
	if ( event.Info()->trackingPOGFilter) continue;
      }
      else
	weight = event.Info()->EventWeight * lumi;

      if ( fabs(DeltaPhi( *event.CaloMET(), *event.TypeIPhiCorrMET()) ) > 1.5) continue;
     
      if ( event.nJets() < 4) continue;
      if ( event.nBJets() == 0) continue;

      for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++)
	METFlag[iMETcut] = event.TypeIPhiCorrMET()->Pt() > METcut[iMETcut];
      
      sigRegionFlag[0] = event.Muons("Selected")->size() == 1 && event.Electrons("Selected")->size() == 0 && fabs( event.FirstLepton()->Eta()) < 2.1;
      sigRegionFlag[1] = event.Muons("Selected")->size() == 0 && event.Electrons("Selected")->size() == 1 && fabs( event.FirstLepton()->Eta()) < 1.5;
      sigRegionFlag[2] = event.Muons("Selected")->size() == 0 && event.Electrons("Selected")->size() == 1 && fabs( event.FirstLepton()->Eta()) < 2.5;
      sigRegionFlag[3] = event.Muons("Selected")->size() + event.Electrons("Selected")->size() == 1 ;

      if ( iSample !=0 && event.Electrons("Selected")->size() == 1 && fabs( event.FirstLepton()->Eta()) >= 2.1 && fabs( event.FirstLepton()->Eta()) < 2.5)
	weight *= gettrigweight( 11, event.FirstLepton()->Pt(), event.FirstLepton()->Eta());

      
      for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
	for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
	  if ( sigRegionFlag[iSigRegion] && METFlag[iMETcut]) {

	    lpth[iSample][iMETcut][iSigRegion]->Fill( event.FirstLepton()->Pt(), weight);
	    letah[iSample][iMETcut][iSigRegion]->Fill( event.FirstLepton()->Eta(), weight);     
	    
	    njetsh[iSample][iMETcut][iSigRegion]->Fill( event.nJets(), weight);  
	    nbjetsh[iSample][iMETcut][iSigRegion]->Fill( event.nBJets(), weight);  

	    jet1h[iSample][iMETcut][iSigRegion]->Fill( event.Jets("Selected")->at(0)->Pt(), weight); 
	    bjet1h[iSample][iMETcut][iSigRegion]->Fill( event.BJets()->at(0)->Pt(), weight);  
	    bjetHighDh[iSample][iMETcut][iSigRegion]->Fill( event.BJetsBDiscOrdered()->at(0)->Pt(), weight);	    
	                                               
	    meth[iSample][iMETcut][iSigRegion]->Fill( event.TypeIPhiCorrMET()->Pt(), weight);     
	    mth[iSample][iMETcut][iSigRegion]->Fill( event.MT(), weight);      
	                                               
	    hth[iSample][iMETcut][iSigRegion]->Fill( event.HT(), weight);      
	    meffh[iSample][iMETcut][iSigRegion]->Fill( event.Meff(), weight);    
	    yh[iSample][iMETcut][iSigRegion]->Fill( event.Y(), weight);       
	                                               
	    mlb1h[iSample][iMETcut][iSigRegion]->Fill( event.Mlb1(), weight);    
	    mlb_hemih[iSample][iMETcut][iSigRegion]->Fill( event.Mlb_hemi(), weight);
	    m3bh[iSample][iMETcut][iSigRegion]->Fill( event.M3b(), weight);     
	    mt2wh[iSample][iMETcut][iSigRegion]->Fill( event.MT2W(), weight);    
	                                               
	    dphiminh[iSample][iMETcut][iSigRegion]->Fill( event.DeltaPhiMinj12m(), weight); 
	    drlbh[iSample][iMETcut][iSigRegion]->Fill( event.DeltaRlb1(), weight);    

	  }
	}
      }
    }
    inFile->Close();
    delete inFile;
  }

  for ( int iSample = 0; iSample < NSamples; iSample++){
    baseDir[iSample]->Write();
  }

  for ( int iSample = 0; iSample < NSamples; iSample++){
    for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
      for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){

	delete lpth[iSample][iMETcut][iSigRegion];
	delete letah[iSample][iMETcut][iSigRegion];
	
	delete njetsh[iSample][iMETcut][iSigRegion];
	delete nbjetsh[iSample][iMETcut][iSigRegion];
	                                           
	delete jet1h[iSample][iMETcut][iSigRegion];     
	delete bjet1h[iSample][iMETcut][iSigRegion];    
	delete bjetHighDh[iSample][iMETcut][iSigRegion];
	
	delete meth[iSample][iMETcut][iSigRegion];
	delete mth[iSample][iMETcut][iSigRegion];
	
	delete hth[iSample][iMETcut][iSigRegion];
	delete meffh[iSample][iMETcut][iSigRegion];
	delete yh[iSample][iMETcut][iSigRegion];
	                                           
	delete mlb1h[iSample][iMETcut][iSigRegion];     
	delete mlb_hemih[iSample][iMETcut][iSigRegion];
	delete m3bh[iSample][iMETcut][iSigRegion];
	delete mt2wh[iSample][iMETcut][iSigRegion];     
	
	delete dphiminh[iSample][iMETcut][iSigRegion];
	delete drlbh[iSample][iMETcut][iSigRegion];    
      }

      //delete metDir[iSample][iMETcut][iSigRegion];
    }

    delete baseDir[iSample];
    
    outFile[iSample]->Close();
    delete outFile[iSample];
  } 
}
  
  //cout<<count<<endl;

  

  /*const int NSamples = 5;
  const int NMETcuts = 3;
  const int NSigRegions = 7;
  

  TString sample[NSamples];
  sample[0] = "Data";
  sample[1] = "WJets";
  sample[2] = "TTJets";
  
  sample[3] = "SingleTop";
  sample[4] = "Rare";
  

  double lumi=19500;
  double weight = 1.;

  double triggerWeight = -1.;
  double puWeight = -1.;
  double globalWeight = -1.;
  double nbjets = -1.;
  double isMu = -1.;
  double lPt = -1.;
  double lEta = -1.;
  double met = -1.;
  double mt = -1.;
  double metDPhi = -1.;	

  bool METFlag[NMETcuts];
  double METcut[NMETcuts] = { 100., 150., 200.};

  bool sigRegionFlag[NSigRegions];
  TString sigRegion[NSigRegions];
  sigRegion[0] = "AllMuons";
  sigRegion[1] = "MuEta0.8";
  sigRegion[2] = "MuEta1.5";
  sigRegion[3] = "MuEndCap";
  sigRegion[4] = "EleEta2.1";
  sigRegion[5] = "EleBarrel";
  sigRegion[6] = "EleEndCap";

  TFile* outFile[NSamples];
  TDirectory* baseDir[NSamples];
  TDirectory* metDir[NSamples][NMETcuts][NSigRegions];
  TH1D* mtHisto[NSamples][NMETcuts][NSigRegions];
  TH1D* metHisto[NSamples][NMETcuts][NSigRegions];
  TH1D* lptHisto[NSamples][NMETcuts][NSigRegions];
  TH1D* letaHisto[NSamples][NMETcuts][NSigRegions];
  TH1D* njetsHisto[NSamples][NMETcuts][NSigRegions];
  TH1D* nbjetsHisto[NSamples][NMETcuts][NSigRegions];

  for ( int iSample = 0; iSample < NSamples; iSample++){
    outFile[iSample] = new TFile( sample[iSample]+"/NoSub/"+sample[iSample]+"_NoSub_STOP_NoTail.root", "RECREATE");
    outFile[iSample]->cd();
    baseDir[iSample] = outFile[iSample]->mkdir("ControlPlots");
    baseDir[iSample]->cd();
    for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
      for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
	TString dirName = "MET"; dirName += METcut[iMETcut]; dirName += "_"; dirName += sigRegion[iSigRegion];
	metDir[iSample][iMETcut][iSigRegion] = baseDir[iSample]->mkdir( dirName);
	metDir[iSample][iMETcut][iSigRegion]->cd();

	//TString name = sample[iSample]; name+="MET"; name+=METcut[iMETcut]; name+=sigRegion[iSigRegion];
	mtHisto[iSample][iMETcut][iSigRegion]      = new TH1D(      "Mt", "Mt", 50, 0., 500.);
	metHisto[iSample][iMETcut][iSigRegion]     = new TH1D(     "MET", "MET", 50, 0., 500.);
	lptHisto[iSample][iMETcut][iSigRegion]     = new TH1D(     "lpt", "Pt of the lepton", 25, 0., 200.);
	letaHisto[iSample][iMETcut][iSigRegion]    = new TH1D(    "leta", "Eta of the lepton", 25, -3., 3.);
	nbjetsHisto[iSample][iMETcut][iSigRegion]  = new TH1D(  "nbjets", "Number of b-tagged jets", 8, 0., 8.);
      }
    }
  }

  for ( int iSample = 0; iSample < NSamples; iSample++) {
    
    std::cout<<"Running over Sample "<<sample[iSample]<<std::endl;
    TString inFileName = "tree/"; inFileName += sample[iSample]; inFileName +=".root";
    
    TFile* inFile = new TFile(inFileName,"READ");
    if (!inFile->IsOpen()){
      std::cout<<"not open"<<std::endl;
    }
    TTree* tree;
    tree= (TTree*)inFile->Get("bonsai");
    
    tree->SetBranchAddress( "trigWeight", &triggerWeight);
    tree->SetBranchAddress( "puWeight", &puWeight);
    tree->SetBranchAddress( "globalWeight", &globalWeight);
    tree->SetBranchAddress( "bjets", &nbjets);
    tree->SetBranchAddress( "isMu", &isMu);
    tree->SetBranchAddress( "lPt", &lPt);
    tree->SetBranchAddress( "lEta", &lEta);
    tree->SetBranchAddress( "met", &met);
    tree->SetBranchAddress( "mt", &mt);
    tree->SetBranchAddress( "metDPhi", &metDPhi);

    int N=tree->GetEntries();
    std::cout<<"looping over "<<N<<" entries"<<std::endl;

    for (int ievt=0;ievt<N;++ievt){
      
      if (ievt%13454 == 0) {
	cout<<"Event number "<<ievt<<"\r"<<flush;
      }
      tree->GetEntry(ievt);
      
      if (iSample == 0)
	weight = 1.;
      else
	weight = triggerWeight * puWeight * globalWeight * lumi;
      
      if ( nbjets > 0.) continue;
      if ( fabs(metDPhi) > 1.5) continue;
      if ( lPt < 30.) continue;
     
      for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++)
	METFlag[iMETcut] = met > METcut[iMETcut];
      
      sigRegionFlag[0] = isMu && fabs(lEta) < 2.1;
      sigRegionFlag[1] = isMu && fabs(lEta) < 0.8;
      sigRegionFlag[2] = isMu && fabs(lEta) < 1.5;
      sigRegionFlag[3] = isMu && fabs(lEta) > 1.5 && fabs(lEta) < 2.1;
      sigRegionFlag[4] = !isMu && fabs(lEta) < 2.1;
      sigRegionFlag[5] = !isMu && fabs(lEta) < 1.4442;
      sigRegionFlag[6] = !isMu && fabs(lEta) < 2.1 && fabs(lEta) > 1.566;  
      
      for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
	for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
	  if ( sigRegionFlag[iSigRegion] && METFlag[iMETcut]) {
	    mtHisto[iSample][iMETcut][iSigRegion]->Fill( mt, weight);
	    metHisto[iSample][iMETcut][iSigRegion]->Fill( met, weight);
	    lptHisto[iSample][iMETcut][iSigRegion]->Fill( lPt, weight);
	    letaHisto[iSample][iMETcut][iSigRegion]->Fill( lEta, weight);
	    nbjetsHisto[iSample][iMETcut][iSigRegion]->Fill( nbjets, weight);
	  }
	}
      }
    }
    inFile->Close();
    delete inFile;
  }
  cout<<"Event number "<<N<<endl;
  
  double peakWeight[NMETcuts][NSigRegions] = {};
  
  double dataPeak[NMETcuts][NSigRegions] = {};
  double wJetsPeak[NMETcuts][NSigRegions] = {};
  double otherMCPeak[NMETcuts][NSigRegions] = {};
  
  int xbin1 = mtHisto[0][0][0]->FindBin( 50.0001);
  int xbin2 = mtHisto[0][0][0]->FindBin( 79.9999);

  for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
    for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
      dataPeak[iMETcut][iSigRegion] = mtHisto[0][iMETcut][iSigRegion]->Integral( xbin1, xbin2);
      wJetsPeak[iMETcut][iSigRegion] = mtHisto[1][iMETcut][iSigRegion]->Integral( xbin1, xbin2);
      otherMCPeak[iMETcut][iSigRegion] = 0.;
      for ( int iSample = 2; iSample < NSamples; iSample++){
	otherMCPeak[iMETcut][iSigRegion] += mtHisto[iSample][iMETcut][iSigRegion]->Integral( xbin1, xbin2);
      }
      
      cout<<dataPeak[iMETcut][iSigRegion]<<" "<<wJetsPeak[iMETcut][iSigRegion]<<" "<<otherMCPeak[iMETcut][iSigRegion]<<endl;

      peakWeight[iMETcut][iSigRegion] = ( dataPeak[iMETcut][iSigRegion] - otherMCPeak[iMETcut][iSigRegion])/ wJetsPeak[iMETcut][iSigRegion];
      
      if (iSigRegion < 4)
	peakWeight[iMETcut][iSigRegion] = peakWeight[iMETcut][0];
      //peakWeight[iMETcut][iSigRegion] = peakWeight[iMETcut][iSigRegion];
      else 	
	peakWeight[iMETcut][iSigRegion] = peakWeight[iMETcut][4];
      //peakWeight[iMETcut][iSigRegion] = peakWeight[iMETcut][iSigRegion];

      mtHisto[1][iMETcut][iSigRegion]->Scale(peakWeight[iMETcut][iSigRegion]);
      metHisto[1][NMETcuts][NSigRegions]->Scale(peakWeight[iMETcut][iSigRegion]);
      lptHisto[1][NMETcuts][NSigRegions]->Scale(peakWeight[iMETcut][iSigRegion]);
      letaHisto[1][NMETcuts][NSigRegions]->Scale(peakWeight[iMETcut][iSigRegion]);
      nbjetsHisto[1][NMETcuts][NSigRegions] ->Scale(peakWeight[iMETcut][iSigRegion]);
    }
  }  

  double SF[NMETcuts][NSigRegions] = {};
  double SFErr[NMETcuts][NSigRegions] = {};

  double dataTail[NMETcuts][NSigRegions] = {};
  double dataTailErr[NMETcuts][NSigRegions] = {};
  double MCTail[NMETcuts][NSigRegions] = {};
  double MCTailErr[NMETcuts][NSigRegions] = {};

  xbin1 = mtHisto[0][0][0]->FindBin( 120.0001);
  xbin2 = mtHisto[0][0][0]->GetNbinsX()+1;

  for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
    for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
      dataTail[iMETcut][iSigRegion] = mtHisto[0][iMETcut][iSigRegion]->IntegralAndError( xbin1, xbin2, dataTailErr[iMETcut][iSigRegion]);
      
      for ( int iSample = 1; iSample < NSamples; iSample++){
	double err = 0.;
	MCTail[iMETcut][iSigRegion] += mtHisto[iSample][iMETcut][iSigRegion]->IntegralAndError( xbin1, xbin2, err);
	MCTailErr[iMETcut][iSigRegion] += err*err;
      }
      MCTailErr[iMETcut][iSigRegion] = sqrt(MCTailErr[iMETcut][iSigRegion]);

      SF[iMETcut][iSigRegion] = dataTail[iMETcut][iSigRegion] / MCTail[iMETcut][iSigRegion];

      SFErr[iMETcut][iSigRegion] = SF[iMETcut][iSigRegion] * sqrt( dataTailErr[iMETcut][iSigRegion] * dataTailErr[iMETcut][iSigRegion] /  dataTail[iMETcut][iSigRegion] / dataTail[iMETcut][iSigRegion] + MCTailErr[iMETcut][iSigRegion] * MCTailErr[iMETcut][iSigRegion] / MCTail[iMETcut][iSigRegion] / MCTail[iMETcut][iSigRegion]);

    }
  }


  for ( int iSample = 0; iSample < NSamples; iSample++){
    for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
      for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
        mtHisto[iSample][iMETcut][iSigRegion]->Rebin(2);
      }
    }
  }

  std::cout<<"\t";
  for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++)
    cout<<sigRegion[iSigRegion]<<"\t";
  std::cout<<std::endl;
  for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
    std::cout<<"MET>"<<METcut[iMETcut]<<"\t";
    for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++)
      std::cout<<peakWeight[iMETcut][iSigRegion]<<"\t";
    std::cout<<std::endl;
  }

  std::cout<<"\t";
  for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++)
    cout<<sigRegion[iSigRegion]<<"\t";
  std::cout<<std::endl;
  for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
    std::cout<<"MET>"<<METcut[iMETcut]<<"\t";
    for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++)
      std::cout<<std::setprecision(3)<<SF[iMETcut][iSigRegion]<<"+-"<<SFErr[iMETcut][iSigRegion]<<"\t";
    std::cout<<std::endl;
  }

   
  for ( int iSample = 0; iSample < NSamples; iSample++){
    baseDir[iSample]->Write();
  }

  for ( int iSample = 0; iSample < NSamples; iSample++){
    for ( int iMETcut = 0; iMETcut < NMETcuts; iMETcut++){
      for ( int iSigRegion = 0; iSigRegion < NSigRegions; iSigRegion++){
	delete mtHisto[iSample][iMETcut][iSigRegion];
	delete metHisto[iSample][iMETcut][iSigRegion];
	delete lptHisto[iSample][iMETcut][iSigRegion];
	delete letaHisto[iSample][iMETcut][iSigRegion];
	delete nbjetsHisto[iSample][iMETcut][iSigRegion];
      }
      delete metDir[iSample][iMETcut][iSigRegion];
    }
    delete baseDir[iSample];

    outFile[iSample]->Close();
    delete outFile[iSample];
  } 

  
}

*/
void IntPlot (TH1D* plot, TH1D* intPlot, float weight){

  TH1::SetDefaultSumw2(true);
  
  int Ntries = 0;
  float totIntegral = 0.;
  float integral = 0.;
  float p=0.;
  float err = 0.;

  Ntries = plot->GetEntries();
  totIntegral = plot->Integral(0, plot->GetNbinsX()+1);

  cout<<Ntries<<" "<<totIntegral<<endl;

  for (int ibin=0; ibin<plot->GetNbinsX()+2; ibin++) {
    integral = plot->Integral(ibin, plot->GetNbinsX()+1);
 
    intPlot->SetBinContent(ibin, integral);

    p=integral/totIntegral;
    err = sqrt (Ntries * p * (1 - p) ) * weight;
    //err = sqrt(integral * weight);

    intPlot->SetBinError(ibin, err);

  } 
     
}
float gettrigweight(int id1, float pt, float eta){

  //electron efficiencies
  if ( abs(id1)==11 ) {
    if ( fabs(eta)<1.5) {
      if ( pt>=20 && pt<22 ) return 0.00;
      if ( pt>=22 && pt<24 ) return 0.00;
      if ( pt>=24 && pt<26 ) return 0.00;
      if ( pt>=26 && pt<28 ) return 0.08;
      if ( pt>=28 && pt<30 ) return 0.61;
      if ( pt>=30 && pt<32 ) return 0.86;
      if ( pt>=32 && pt<34 ) return 0.88;
      if ( pt>=34 && pt<36 ) return 0.90;
      if ( pt>=36 && pt<38 ) return 0.91;
      if ( pt>=38 && pt<40 ) return 0.92;
      if ( pt>=40 && pt<50 ) return 0.94;
      if ( pt>=50 && pt<60 ) return 0.95;
      if ( pt>=60 && pt<80 ) return 0.96;
      if ( pt>=80 && pt<100 ) return 0.96;
      if ( pt>=100 && pt<150 ) return 0.96;
      if ( pt>=150 && pt<200 ) return 0.97;
      if ( pt>=200 ) return 0.97;
    } else if ( fabs(eta)>=1.5 && fabs(eta)<2.5) {
      if ( pt>=20 && pt<22 ) return 0.00;
      if ( pt>=22 && pt<24 ) return 0.00;
      if ( pt>=24 && pt<26 ) return 0.02;
      if ( pt>=26 && pt<28 ) return 0.18;
      if ( pt>=28 && pt<30 ) return 0.50;
      if ( pt>=30 && pt<32 ) return 0.63;
      if ( pt>=32 && pt<34 ) return 0.68;
      if ( pt>=34 && pt<36 ) return 0.70;
      if ( pt>=36 && pt<38 ) return 0.72;
      if ( pt>=38 && pt<40 ) return 0.74;
      if ( pt>=40 && pt<50 ) return 0.76;
      if ( pt>=50 && pt<60 ) return 0.77;
      if ( pt>=60 && pt<80 ) return 0.78;
      if ( pt>=80 && pt<100 ) return 0.80;
      if ( pt>=100 && pt<150 ) return 0.79;
      if ( pt>=150 && pt<200 ) return 0.76;
      if ( pt>=200 ) return 0.81;
    }
  } else if ( abs(id1)==13 ) {//muon efficiencies

    if ( fabs(eta)<0.8 ) {
      if (pt>=20 && pt<22)  return  0.00;     
      if (pt>=22 && pt<24)  return  0.03;      
      if (pt>=24 && pt<26)  return  0.87;
      if (pt>=26 && pt<28)  return  0.90;
      if (pt>=28 && pt<30)  return  0.91;
      if (pt>=30 && pt<32)  return  0.91;
      if (pt>=32 && pt<34)  return  0.92;
      if (pt>=34 && pt<36)  return  0.93;
      if (pt>=36 && pt<38)  return  0.93;
      if (pt>=38 && pt<40)  return  0.93;
      if (pt>=40 && pt<50)  return  0.94;
      if (pt>=50 && pt<60)  return  0.95;
      if (pt>=60 && pt<80)  return  0.95;
      if (pt>=80 && pt<100) return 0.94;
      if (pt>=100 && pt<150) return 0.94;
      if (pt>=150 && pt<200) return 0.93;
      if (pt>=200) return 0.92;
    } else if ( fabs(eta)>=0.8 && fabs(eta)<1.5 ) {
      if (pt>=20 && pt<22)  return  0.00;
      if (pt>=22 && pt<24)  return  0.05;
      if (pt>=24 && pt<26)  return  0.78;
      if (pt>=26 && pt<28)  return  0.81;
      if (pt>=28 && pt<30)  return  0.81;
      if (pt>=30 && pt<32)  return  0.81;
      if (pt>=32 && pt<34)  return  0.82;
      if (pt>=34 && pt<36)  return  0.82;
      if (pt>=36 && pt<38)  return  0.83;
      if (pt>=38 && pt<40)  return  0.83;
      if (pt>=40 && pt<50)  return  0.84;
      if (pt>=50 && pt<60)  return  0.84;
      if (pt>=60 && pt<80)  return  0.84;
      if (pt>=80 && pt<100) return 0.84;
      if (pt>=100 && pt<150) return 0.84;
      if (pt>=150 && pt<200) return 0.84;
      if (pt>=200) return 0.82;
    } else if ( fabs(eta)>=1.5 && fabs(eta)<2.1 ) {
      if (pt>=20 && pt<22)  return  0.00;
      if (pt>=22 && pt<24)  return  0.11;
      if (pt>=24 && pt<26)  return  0.76;
      if (pt>=26 && pt<28)  return  0.78;
      if (pt>=28 && pt<30)  return  0.79;
      if (pt>=30 && pt<32)  return  0.80;
      if (pt>=32 && pt<34)  return  0.80;
      if (pt>=34 && pt<36)  return  0.81;
      if (pt>=36 && pt<38)  return  0.81;
      if (pt>=38 && pt<40)  return  0.82;
      if (pt>=40 && pt<50)  return  0.82;
      if (pt>=50 && pt<60)  return  0.83;
      if (pt>=60 && pt<80)  return  0.83;
      if (pt>=80 && pt<100) return 0.83;
      if (pt>=100 && pt<150) return 0.83;
      if (pt>=150 && pt<200) return 0.82;
      if (pt>=200) return 0.82;
    }
  }//end check for muons

  return 1.;

}


