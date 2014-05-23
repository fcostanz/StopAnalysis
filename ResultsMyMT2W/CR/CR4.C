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
  const int NJet = 4;

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

  Float_t N[NJet] = {};  Float_t NErr[NJet] = {};
  Float_t M[NJet] = {};  Float_t MErr[NJet] = {};
  Float_t SFNJet[NJet] = {}; Float_t SFNJetErr[NJet] = {};
  
  Float_t K[NJet - 1] = {};  Float_t KErr[NJet - 1] = {};
  Float_t KComb = 0.; Float_t KCombErr = 0.;

  Float_t alpha = 0.;  Float_t alphaErr = 0.;

  Float_t tailData = 0.; Float_t tailDataErr = 0.;
  Float_t tailMC = 0.; Float_t tailMCErr = 0.;

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
    tree[ilep]->Branch( "N5", &N[3]);
    tree[ilep]->Branch( "N5Err", &NErr[3]);    

    tree[ilep]->Branch( "M2", &M[0]);
    tree[ilep]->Branch( "M2Err", &MErr[0]);
    tree[ilep]->Branch( "M3", &M[1]);
    tree[ilep]->Branch( "M3Err", &MErr[1]);
    tree[ilep]->Branch( "M4", &M[2]);
    tree[ilep]->Branch( "M4Err", &MErr[2]);
    tree[ilep]->Branch( "M5", &M[3]);
    tree[ilep]->Branch( "M5Err", &MErr[3]);  

    tree[ilep]->Branch( "SF2", &SFNJet[0]);
    tree[ilep]->Branch( "SF2Err", &SFNJetErr[0]);
    tree[ilep]->Branch( "SF3", &SFNJet[1]);
    tree[ilep]->Branch( "SF3Err", &SFNJetErr[1]);
    tree[ilep]->Branch( "SF4", &SFNJet[2]);
    tree[ilep]->Branch( "SF4Err", &SFNJetErr[2]);
    tree[ilep]->Branch( "SF5", &SFNJet[3]);
    tree[ilep]->Branch( "SF5Err", &SFNJetErr[3]);
    
    tree[ilep]->Branch( "K3", &K[0]);
    tree[ilep]->Branch( "K3Err", &KErr[0]);
    tree[ilep]->Branch( "K4", &K[1]);
    tree[ilep]->Branch( "K4Err", &KErr[1]);
    tree[ilep]->Branch( "K5", &K[2]);
    tree[ilep]->Branch( "K5Err", &KErr[2]);
    tree[ilep]->Branch( "K", &KComb);
    tree[ilep]->Branch( "KErr", &KCombErr);

    tree[ilep]->Branch( "alpha", &alpha);
    tree[ilep]->Branch( "alphaErr", &alphaErr);

    tree[ilep]->Branch( "tailData", &tailData); 
    tree[ilep]->Branch( "tailDataErr", &tailDataErr);
    tree[ilep]->Branch( "tailMC", &tailMC);
    tree[ilep]->Branch( "tailMCErr", &tailMCErr);
    
    tree[ilep]->Branch( "SF", &SF);
    tree[ilep]->Branch( "SFErr", &SFErr);
  }

  TString histoName = "";

  TH1D* nh[NSamples][NJet];
  TH1D* Nh[NJet];
  TH1D* Mh[NJet];
  TH1D* SFNJeth[NJet];  
  for ( int i = 0; i < NJet; i++){
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
  
  TH1D* NCombh = new TH1D( "NComb", "NComb", 1, 0., 1.);
  TH1D* MCombh = new TH1D( "MComb", "MComb", 1, 0., 1.);

  TH1D* Kh[NJet - 1];
  for ( int i = 0; i < NJet - 1; i++){
    histoName = "K"; histoName += i + 3; 
    Kh[i] = new TH1D( histoName, histoName, 1, 0., 1.);
  }
  TH1D* KCombh = new TH1D( "KComb", "KComb", 1, 0., 1.);  

  TH1D* alphah = new TH1D( "alpha", "alpha", 1, 0., 1.);

  TH1D* tailh[NSamples];

  TH1D* SFh = new TH1D( "SF", "SF", 1, 0., 1.);

  for ( int iSample = 0; iSample < NSamples; iSample++)
    tailh[iSample] = new TH1D( TString("tail") + sample[iSample], TString("tail") + sample[iSample], 1, 0., 1.);

  int NSR = srTree->GetEntries();

  TFile* inFile[NSamples];
  TDirectory* inBaseDir[NSamples];
  TDirectory* metDir[NSamples];
  TH1D* mtDatah;
  TH1D* mtMCh;
  TH1D* mtnh[NSamples][NJet + 1]; 
  TH1D* njetsh[NSamples];

  TH1D* wToph = new TH1D( "wToph", "wToph", 1, 0., 1.);
  TH1D* woToph = new TH1D( "woToph", "woToph", 1, 0., 1.);

  for (int iSample = 0; iSample < NSamples; iSample++){
    TString inFileName = "../MakeHistos/"+sample[iSample]+".root";
    inFile[iSample]= new TFile( inFileName,"READ");
  }

  TH1D* numh = new TH1D( "num", "num", 1, 0., 1.);      
  TH1D* denh = new TH1D( "den", "den", 1, 0., 1.);
  TH1D* tmph = new TH1D( "tmp", "tmp", 1, 0., 1.);

  for (int ilep = 0; ilep < NLep; ilep++){
    for ( int iSR = 0; iSR < NSR; iSR++){      
      srTree->GetEntry(iSR);
      
      if (metCut > 260) metCut = 250.;
      

      TString baseDirName = ""; baseDirName += iSR;
      for (int iSample = 0; iSample < NSamples; iSample++)
	inBaseDir[iSample] =  inFile[iSample]->GetDirectory(baseDirName);
      
      TString metDirName = "";       
      Float_t metTmp = metCut;
      for ( int jSR = NSR-5; jSR < NSR; jSR++){
	srTree->GetEntry(jSR);
	if ( fabs(metCut - metTmp) < 0.001){
	  metDirName += jSR;
	  break;
	}
      }

      srTree->GetEntry(iSR);
      for (int iSample = 0; iSample < NSamples; iSample++)
	metDir[iSample] =  inFile[iSample]->GetDirectory(metDirName);

      ////////////////////
      // Cleaning
      ////////////////////
      
      numh->Clear(); numh->Reset();
      denh->Clear(); denh->Reset();
      tmph = new TH1D( "tmp", "tmp", 1, 0., 1.);
      tmph->Clear(); tmph->Reset();

      for ( int i = 0; i < NJet; i++){
	Nh[i]->Clear();   Nh[i]->Reset();
	Mh[i]->Clear();   Mh[i]->Reset();
	SFNJeth[i]->Clear(); SFNJeth[i]->Reset();
	
	for ( int iSample = 0; iSample < NSamples; iSample++){
	  nh[iSample][i]->Clear();
	  nh[iSample][i]->Reset();
	}
      }
      NCombh->Clear();   NCombh->Reset();
      MCombh->Clear();   MCombh->Reset();
      
      for ( int i = 0; i < NJet - 1; i++){
	Kh[i]->Clear();
	Kh[i]->Reset();
      }
      KCombh->Clear();  KCombh->Reset();

      alphah->Clear();  alphah->Reset();
      
      SFh->Clear(); SFh->Reset();
      for ( int iSample = 0; iSample < NSamples; iSample++){
	tailh[iSample]->Clear();
	tailh[iSample]->Reset();	
      }

      ////////////////////
      // Reading input
      ////////////////////

      histoName = lep[ilep]; histoName += "-CR4/wTopPtReweight";	
      inBaseDir[1]->GetObject( histoName, wToph);
      inBaseDir[2]->GetObject( histoName, tmph);
      wToph->Add( tmph);     
      
      histoName = lep[ilep]; histoName += "-CR4/woTopPtReweight";	
      inBaseDir[1]->GetObject( histoName, woToph);
      inBaseDir[2]->GetObject( histoName, tmph);
      woToph->Add( tmph);

      for (int iSample = 0; iSample < NSamples; iSample++){
	for ( int i = 0; i < NJet; i++){
	  histoName = lep[ilep]; histoName += "-CR4/Mt"; histoName += (i + 1); histoName += "Jet";
	  inBaseDir[iSample]->GetObject( histoName, mtnh[iSample][i]);
	  
	  if (iSample == 1 || iSample == 2) {
	    if (wToph->GetBinContent(1) > 0.)
	      mtnh[iSample][i]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
	  }
	}
	histoName = lep[ilep]; histoName += "-CR4/Mt5+Jet";
	inBaseDir[iSample]->GetObject( histoName, mtnh[iSample][4]);
	if (iSample == 1 || iSample == 2) {
	  if (wToph->GetBinContent(1) > 0.)
	    mtnh[iSample][4]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
	}

	metDir[iSample]->GetObject( lep[ilep] + TString("-CR4/njets"), njetsh[iSample]);

	nh[iSample][0]->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(1)));
	nh[iSample][0]->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(1)));
	tmph->Clear(); tmph->Reset();
	tmph->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(2)));
	tmph->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(2)));
	nh[iSample][0]->Add(tmph);

	nh[iSample][1]->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(3)));
	nh[iSample][1]->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(3)));

	nh[iSample][2]->SetBinContent( 1, njetsh[iSample]->GetBinContent( njetsh[iSample]->FindBin(4)));
	nh[iSample][2]->SetBinError( 1, njetsh[iSample]->GetBinError( njetsh[iSample]->FindBin(4)));

	n1 = njetsh[iSample]->FindBin(5);
	n2 = njetsh[iSample]->GetNbinsX() + 1;
	integral = njetsh[iSample]->IntegralAndError( n1, n2, error);
	nh[iSample][3]->SetBinContent(1, integral);
	nh[iSample][3]->SetBinError(1, error);
      }	

      ///////////////////////////////////////////
      // N, M and K
      ///////////////////////////////////////////
      for ( int i = 0; i < NJet; i++){
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

      for ( int i = 0; i < NJet-1; i++){
	Kh[i]->Add(SFNJeth[i+1]);
	Kh[i]->Divide(SFNJeth[0]);

	K[i] = Kh[i]->GetBinContent(1);
	KErr[i] = Kh[i]->GetBinError(1);     
      }

      Int_t ijet = 0;

      if ( fabs(nJetCut - 3) < 0.001 ) ijet = 1;
      if ( fabs(nJetCut - 4) < 0.001 ) ijet = 2;
      if ( fabs(nJetCut - 5) < 0.001 ) ijet = 3;
      
      for ( int i = ijet; i < NJet; i++){
	NCombh->Add(nh[0][i]);
	NCombh->Add(nh[2][i], -1.);
	NCombh->Add(nh[3][i], -1.);
	NCombh->Add(nh[4][i], -1.);
	
	MCombh->Add(nh[1][i]);
      }
      
      KCombh->Add(NCombh);
      KCombh->Divide(MCombh);
      KCombh->Divide(SFNJeth[0]);

      KComb = KCombh->GetBinContent(1);
      KCombErr = KCombh->GetBinError(1);

      ///////////////////////////////////////////////
      // Mt shape
      ///////////////////////////////////////////////

      mtDatah = new TH1D(*mtnh[0][1]);
      for ( int i = 2; i < NJet + 1; i++)
	mtDatah->Add(mtnh[0][i]);

      mtMCh = new TH1D(*mtnh[0][0]);
      mtMCh->Clear(); mtMCh->Reset(); 

      n1 = 0;      
      n2 = mtDatah->GetNbinsX() + 1;

      for (int iSample = 1; iSample < NSamples; iSample++){
	for ( int i = 1; i < NJet + 1; i++){
	  TH1D* tmpMth = new TH1D(*mtnh[iSample][i]);
	  if( (i - 2.) >= -0.001){
	    tmpMth->Scale(K[i-2]);
	    Int_t j = i - 2;
	    for (int ibin = 0; ibin < tmph->GetNbinsX()+2; ibin++){
	      Float_t Mt = mtnh[iSample][i]->GetBinContent( ibin);
	      Float_t errMt = mtnh[iSample][i]->GetBinError( ibin);	      
	      
	      error = sqrt( Mt * Mt * KErr[j] * KErr[j] + K[j] * K[j] * errMt * errMt);	      
	      tmpMth->SetBinError( ibin, error);
	    }
	  }

	  mtMCh->Add(tmpMth);
 
	  delete tmpMth;
	}
      }
	
      n1 = 0;      
      n2 = mtDatah->GetNbinsX() + 1;

      integral = mtDatah->IntegralAndError( n1, n2, error);
      alphah->SetBinContent( 1, integral);
      alphah->SetBinError( 1, error);

      integral = mtMCh->IntegralAndError( n1, n2, error);
      tmph->SetBinContent( 1, integral);
      tmph->SetBinError( 1, error);

      alphah->Divide(tmph);
      
      alpha = alphah->GetBinContent(1);
      alphaErr = alphah->GetBinError(1);

      n3 = mtDatah->FindBin(mtCut);      
      n4 = mtDatah->GetNbinsX() + 1;

      integral = mtDatah->IntegralAndError( n3, n4, error);
      SFh->SetBinContent( 1, integral);
      SFh->SetBinError( 1, error);

      tailData = integral;
      tailDataErr = error;

      integral = mtMCh->IntegralAndError( n3, n4, error);
      tmph->SetBinContent( 1, integral);
      tmph->SetBinError( 1, error);

      tmph->Multiply(alphah);

      tailMC = tmph->GetBinContent(1);
      tailMCErr = tmph->GetBinError(1);

      SFh->Divide(tmph);
      
      SF = SFh->GetBinContent(1);
      SFErr = SFh->GetBinError(1);

      delete mtDatah;
      delete mtMCh;

      ///////////////////////////////////////////////
      // Fill outTree
      ///////////////////////////////////////////////

      tree[ilep]->Fill();
      
      //////////////////////////////////
      // Pritn some output
      //////////////////////////////////

      cout<<ID<<"\t";
      cout<<iSR<<" "<<ilep<<"\t";
      //for ( int i = 0; i < NJet; i++)
      //cout<<SFNJet[i]<<"+-"<<SFNJetErr[i]<<"\t";
      cout<<K[0]<<"+-"<<KErr[0]<<"\t";
      cout<<K[1]<<"+-"<<KErr[1]<<"\t";  
      cout<<K[2]<<"+-"<<KErr[2]<<"\t";
      cout<<KComb<<"+-"<<KCombErr<<"\t";
      cout<<tailData<<"+-"<<tailDataErr<<"\t";
      cout<<tailMC<<"+-"<<tailMCErr<<"\t";
      cout<<alpha<<"+-"<<alphaErr<<"\t";
      cout<<SF<<"+-"<<SFErr<<endl;    
    }      
  }
  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();
  
  outFile->cd();
  for (int ilep = 0; ilep < NLep; ilep++)
    tree[ilep]->Write();

  return 0;
}
