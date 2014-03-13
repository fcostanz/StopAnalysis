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

int T2tbCheck(int iSample){
  
  TH1::SetDefaultSumw2(true);
  
  if(pcp)cout<<"going to set inputs"<<endl;

  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/KinVariables/";

  const int NSignals = 7;
  
  TString sample[NSignals];
  sample[0]="T2tb-mStop175mLSP50";
  sample[1]="T2tb-mStop200mLSP25";
  sample[2]="T2tb-mStop375mLSP50";
  sample[3]="T2tb-mStop250mLSP25";
  sample[4]="T2tb-mStop325mLSP100";
  sample[5]="T2tb-mStop450mLSP150";
  sample[6]="T2tb-mStop550mLSP1";

  double lumi=19500.;

  TFile* outFile = new TFile( sample[iSample]+"_outfile.root", "RECREATE"); 

  outFile->cd();
  TH1D* njets_NoChargino = new TH1D( "njets_NoChargino", "jets multiplicity",   7, 0.5, 7.5);
  TH1D* lpt_NoChargino   = new TH1D(   "lpt_NoChargino",   "lep p_{T} [GeV]",  20, 0., 100.);
  
  TH1D* njets_1Chargino_lepFromTop = new TH1D( "njets_1Chargino_lepFromTop", "jets multiplicity",   7, 0.5, 7.5);
  TH1D* lpt_1Chargino_lepFromTop   = new TH1D(   "lpt_1Chargino_lepFromTop",   "lep p_{T} [GeV]",  20, 0., 100.);
  
  TH1D* njets_1Chargino_lepFromChi = new TH1D( "njets_1Chargino_lepFromChi", "jets multiplicity",   7, 0.5, 7.5);
  TH1D* lpt_1Chargino_lepFromChi   = new TH1D(   "lpt_1Chargino_lepFromChi",   "lep p_{T} [GeV]",  20, 0., 100.);

  TH1D* nLepFromTop_1Chargino = new TH1D(  "nLepFromTop_1Chargino", "lep From Top", 4, -0.5, 3.5);

  
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
    
  Float_t weight = 0.;
  Float_t globalWeight = 0.;
  
  Float_t lepFromTop = 0.;
  Float_t charginos = 0.;
  
  Float_t lPt = 0.;
  Float_t njets = 0.;
  Float_t mt = 0.;

  Char_t kinRegion = false;

  tree->SetBranchAddress( "GlobalWeight", &globalWeight);

  tree->SetBranchAddress( "LepFromTop", &lepFromTop);
  tree->SetBranchAddress( "Charginos", &charginos);

  tree->SetBranchAddress( "lPt", &lPt);
  tree->SetBranchAddress( "njets", &njets);
  tree->SetBranchAddress( "mt", &mt);

  tree->SetBranchAddress("kinRegion",&kinRegion);

  for (int ievt=0;ievt<N;++ievt){
    if (ievt%13454 == 0) {
      cout<<"Event number "<<ievt<<"\r"<<flush;
    }
    
    tree->GetEntry(ievt);

    weight = globalWeight * lumi;

    if (!kinRegion)
      continue;

    if ( mt < 120.)
      continue;

    if (charginos == 0){
      njets_NoChargino->Fill(njets, weight);
      lpt_NoChargino->Fill(lPt, weight);
    }

    else if (charginos == 1){
      nLepFromTop_1Chargino->Fill( lepFromTop, weight);
      
      if ( lepFromTop == 0){
	njets_1Chargino_lepFromChi->Fill(njets, weight);
	lpt_1Chargino_lepFromChi->Fill(lPt, weight);
      }
      else if ( lepFromTop == 1){
	njets_1Chargino_lepFromTop->Fill(njets, weight);
	lpt_1Chargino_lepFromTop->Fill(lPt, weight);
      }
      else{
	cout<<"Houston, we've got a problem!"<<endl;
      }
    }
  }

  outFile->cd();
  nLepFromTop_1Chargino->Write();
        
  njets_NoChargino->Write();
  lpt_NoChargino->Write();
  
  njets_1Chargino_lepFromChi->Write();
  lpt_1Chargino_lepFromChi->Write();

  njets_1Chargino_lepFromTop->Write();
  lpt_1Chargino_lepFromTop->Write();

  outFile->Write();
  outFile->Close();
    
  inFile->Close();
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


