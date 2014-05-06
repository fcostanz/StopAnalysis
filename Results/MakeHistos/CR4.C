#include "TROOT.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TCollection.h"
#include "TKey.h"

//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

bool pcp = true;

int CR4(){
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;

  TFile* srFile = new TFile( "../SearchRegions/SearchRegions.root", "READ"); 
  TTree* srTree;  
  srFile->GetObject( "SearchRegions", srTree);
  
  Int_t ID = 0;
  Float_t nJetCut = 0.;
  Float_t mtCut = 0.;
  Float_t dphiCut = 0.;
  Float_t centralityCut = 0.;
  Float_t metCut = 0.;
  Float_t yCut = 0.;
  Float_t mlbCut = 0.;
  Float_t m3Cut = 0.;
  Float_t mt2wCut = 0.;
  Float_t drlblCut = 10.;

  srTree->SetBranchAddress( "ID", &ID);
  srTree->SetBranchAddress( "nJetCut", &nJetCut);
  srTree->SetBranchAddress( "mtCut", &mtCut);
  srTree->SetBranchAddress( "dphiCut", &dphiCut);
  srTree->SetBranchAddress( "centralityCut", &centralityCut);
  srTree->SetBranchAddress( "metCut", &metCut);
  srTree->SetBranchAddress( "yCut", &yCut);
  srTree->SetBranchAddress( "mlbCut", &mlbCut);
  srTree->SetBranchAddress( "m3Cut", &m3Cut);
  srTree->SetBranchAddress( "mt2wCut", &mt2wCut);

  const int NSamples = 5;  
  const int NLep = 3;

  TString sample[NSamples];
  sample[0] = "Data";
  sample[1] = "DiLep";
  sample[2] = "OneLep";
  sample[3] = "WJets";
  sample[4] = "Rare";

  bool lepFlag[NLep];
  TString lep[NLep];
  lep[0] = "El";
  lep[1] = "Mu";
  lep[2] = "ElAndMu";

  Int_t n1 = 0;
  Int_t n2 = 0;
  Int_t n3 = 0;
  Int_t n4 = 0;
    	
  Double_t integral = 0.;
  Double_t error = 0.;

  TString outFileName = "./CR4.root";
  TFile* outFile = new TFile( outFileName, "RECREATE");
  outFile->cd();

  Float_t SFLep[NLep] = {};

  Float_t SFLepErr[NLep] = {};

  Float_t N[3] = {};  Float_t NErr[3] = {};
  Float_t M[3] = {};  Float_t MErr[3] = {};
  Float_t SFNJet[3] = {}; Float_t SFNJetErr[3] = {};
  
  Float_t K[2] = {};  Float_t KErr[2] = {};

  Float_t alpha = 0.;  Float_t alphaErr = 0.;

  Float_t tailData = 0.; Float_t tailDataErr = 0.;
  Float_t tailMC = 0.; Float_t tailMCErr = 0.;
  Float_t tailDiLep = 0.; Float_t tailDiLepErr = 0.;

  Float_t SF = 0.; Float_t SFErr = 0.;
  
  TTree* tree[3];  
  for(int ilep = 0; ilep < NLep; ilep++){
    tree[ilep] = new TTree( lep[ilep], lep[ilep]);

    tree[ilep]->Branch( "N2", &N[0]);
    tree[ilep]->Branch( "N2Err", &NErr[0]);
    tree[ilep]->Branch( "N3", &N[1]);
    tree[ilep]->Branch( "N3Err", &NErr[1]);
    tree[ilep]->Branch( "N4", &N[2]);
    tree[ilep]->Branch( "N4Err", &NErr[2]);

    tree[ilep]->Branch( "M2", &M[0]);
    tree[ilep]->Branch( "M2Err", &MErr[0]);
    tree[ilep]->Branch( "M3", &M[1]);
    tree[ilep]->Branch( "M3Err", &MErr[1]);
    tree[ilep]->Branch( "M4", &M[2]);
    tree[ilep]->Branch( "M4Err", &MErr[2]);

    tree[ilep]->Branch( "SF2", &SFNJet[0]);
    tree[ilep]->Branch( "SF2Err", &SFNJetErr[0]);
    tree[ilep]->Branch( "SF3", &SFNJet[1]);
    tree[ilep]->Branch( "SF3Err", &SFNJetErr[1]);
    tree[ilep]->Branch( "SF4", &SFNJet[2]);
    tree[ilep]->Branch( "SF4Err", &SFNJetErr[2]);
  
    tree[ilep]->Branch( "K3", &K[1]);
    tree[ilep]->Branch( "K3Err", &KErr[1]);
    tree[ilep]->Branch( "K4", &K[2]);
    tree[ilep]->Branch( "K4Err", &KErr[2]);

    tree[ilep]->Branch( "alpha", &alpha);
    tree[ilep]->Branch( "alphaErr", &alphaErr);

    tree[ilep]->Branch( "tailData", &tailData); 
    tree[ilep]->Branch( "tailDataErr", &tailDataErr);
    tree[ilep]->Branch( "tailMC", &tailMC);
    tree[ilep]->Branch( "tailMCErr", &tailMCErr);
    tree[ilep]->Branch( "tailDiLep", &tailDiLep);
    tree[ilep]->Branch( "tailDiLepErr", &tailDiLepErr);
    
    tree[ilep]->Branch( "SF", &SF);
    tree[ilep]->Branch( "SFErr", &SFErr);
  }

  TString histoName = "";

  TH1D* nh[NSamples][3];
  TH1D* Nh[3];
  TH1D* Mh[3];
  TH1D* SFNJeth[3];  
  for ( int i = 0; i < 3; i++){
    histoName = "N"; histoName += i + 1;          
    Nh[i] = new TH1D( histoName, histoName, 1, 0., 1.);
    
    histoName = "M"; histoName += i + 1;          
    Mh[i] = new TH1D( histoName, histoName, 1, 0., 1.);

    histoName = "SFNJet"; histoName += i + 1;          
    SFNJeth[i] = new TH1D( histoName, histoName, 1, 0., 1.);

    histoName = "n"; histoName += i + 1;
    for ( int iSample = 0; iSample < NSamples; iSample++)
      nh[iSample][i] = new TH1D( histoName + sample[iSample], histoName + sample[iSample], 1, 0., 1.);
  }
  
  TH1D* Kh[2];
  for ( int i = 0; i < 3; i++){
    histoName = "K"; histoName += i + 2; 
    Kh[i] = new TH1D( histoName, histoName, 1, 0., 1.);
  }

  TH1D* alphah = new TH1D( "alpha", "alpha", 1, 0., 1.);

  TH1D* tailh[NSamples];

  TH1D* SFh = new TH1D( "SF", "SF", 1, 0., 1.);

  for ( int iSample = 0; iSample < NSamples; iSample++)
    tailh[iSample] = new TH1D( TString("tail") + sample[iSample], TString("tail") + sample[iSample], 1, 0., 1.);

  int NSR = srTree->GetEntries();

  TFile* inFile[NSamples];
  TDirectory* inBaseDir[NSamples];
  TH1D* mtDatah;
  TH1D* mtMCh;  
  TH1D* mth[NSamples]; 
  TH1D* mt1h[NSamples];
  TH1D* mt2h[NSamples];
  TH1D* mt3h[NSamples];
  TH1D* mt4h[NSamples];
  TH1D* njetsh[NSamples];

  for (int iSample = 0; iSample < NSamples; iSample++){
    TString inFileName = "./MakeHistos/"+sample[iSample]+".root";
    inFile[iSample]= new TFile( inFileName,"READ");
  }

  TH1D* numh = new TH1D( "num", "num", 1, 0., 1.);      
  TH1D* denh = new TH1D( "den", "den", 1, 0., 1.);
  TH1D* tmph = new TH1D( "tmp", "tmp", 1, 0., 1.);

  for (int ilep = 0; ilep < NLep; ilep++){
    for ( int iSR = 0; iSR < NSR; iSR++){      
      srTree->GetEntry(iSR);
      
      TString baseDirName = ""; baseDirName += iSR;
      for (int iSample = 0; iSample < NSamples; iSample++)
	inBaseDir[iSample] =  inFile[iSample]->GetDirectory(baseDirName);
      
      ////////////////////
      // Cleaning
      ////////////////////
      
      numh->Clear(); numh->Reset();
      denh->Clear(); denh->Reset();
      tmph->Clear(); tmph->Reset();

      for ( int i = 0; i < 3; i++){
	Nh[i]->Clear();   Nh[i]->Reset();
	Mh[i]->Clear();   Mh[i]->Reset();
	SFNJeth[i]->Clear(); SFNJeth[i]->Reset();
	
	for ( int iSample = 0; iSample < NSamples; iSample++){
	  nh[iSample][i]->Clear();
	  nh[iSample][i]->Reset();
	}
      }

      for ( int i = 0; i < 2; i++){
	Kh[i]->Clear();
	Kh[i]->Reset();
      }
	
      alphah->Clear();  alphah->Reset();
      
      SFh->Clear(); SFh->Reset();
      for ( int iSample = 0; iSample < NSamples; iSample++){
	tailh[iSample]->Clear();
	tailh[iSample]->Reset();	
      }

      ////////////////////
      // Reading input
      ////////////////////
  
      for (int iSample = 0; iSample < NSamples; iSample++){
	inBaseDir[iSample]->GetObject( lep[ilep] + TString("-CR4/Mt1Jet"), mt1h[iSample]);
	inBaseDir[iSample]->GetObject( lep[ilep] + TString("-CR4/Mt2Jet"), mt2h[iSample]);	
	inBaseDir[iSample]->GetObject( lep[ilep] + TString("-CR4/Mt3Jet"), mt3h[iSample]);	
	inBaseDir[iSample]->GetObject( lep[ilep] + TString("-CR4/Mt4Jet"), mt4h[iSample]);	

	inBaseDir[iSample]->GetObject( lep[ilep] + TString("-CR4/njets"), njetsh[iSample]);

	nh[iSample][0]->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(1)));
	nh[iSample][0]->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(1)));
	tmph->Clear(); tmph->Reset();
	tmph->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(2)));
	tmph->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(2)));
	nh[iSample][0]->Add(tmph);

	nh[iSample][1]->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(3)));
	nh[iSample][1]->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(3)));

	n1 = njetsh[iSample]->FindBin(4);
	n2 = njetsh[iSample]->GetNbinsX() + 1;
	integral = njetsh[iSample]->IntegralAndError( n1, n2, error);
	nh[iSample][2]->SetBinContent(1, integral);
	nh[iSample][2]->SetBinError(1, error);
      }

      ///////////////////////////////////////////
      // N, M and K
      ///////////////////////////////////////////

      for ( int i = 0; i < 3; i++){
	Nh[i]->Add(nh[0][i]);
	Nh[i]->Add(nh[2][i], -1.);
	Nh[i]->Add(nh[3][i], -1.);
	Nh[i]->Add(nh[4][i], -1.);

	Mh[i]->Add(nh[1][i]);

	SFNJeth[i]->Add(Nh[i]);
	SFNJeth[i]->Divide(Mh[i]);

	N[i] = Nh[i]->GetBinContent(1);
	NErr[i] = Nh[i]->GetBinError(1);
	
	M[i] = Mh[i]->GetBinContent(1);
	MErr[i] = Mh[i]->GetBinError(1);
	
	SFNJet[i] = SFNJeth[i]->GetBinContent(1);
	SFNJetErr[i] = SFNJeth[i]->GetBinError(1);
      }
      
      for ( int i = 0; i < 2; i++){
	Kh[i]->Add(SFNJeth[i+1]);
	Kh[i]->Divide(SFNJeth[0]);
	
	K[i] = Kh[i]->GetBinContent(1);
	KErr[i] = Kh[i]->GetBinError(1);     
      }

      ///////////////////////////////////////////////
      // Mt shape
      ///////////////////////////////////////////////
      /*
      for (int iSample = 1; iSample < NSamples; iSample++){
	mt->Add(mt1h[iSample]);
	
	
	mtMCh = new TH1D( mt);
	
	mt->Clear(); mt->Reset();
	mt->Add(mt2h[iSample]);
	mt->Multiply(mt2h[iSample]);	


      */


      ///////////////////////////////////////////////
      // Fill outTree
      ///////////////////////////////////////////////

      tree[ilep]->Fill();
      
      //////////////////////////////////
      // Pritn some output
      //////////////////////////////////

      cout<<ID<<"\t";
      cout<<iSR<<" "<<ilep<<"\t";
      for ( int i = 0; i < 3; i++)
	cout<<SFNJet[i]<<"+-"<<SFNJetErr[i]<<"\t";
      cout<<K[0]<<"+-"<<KErr[0]<<"\t";
      cout<<K[1]<<"+-"<<KErr[1]<<endl;        
    }      
  }

  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();

  outFile->cd();
  for (int ilep = 0; ilep < NLep; ilep++)
    tree[ilep]->Write();
  
  return 0;
}
