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

int CR5(){
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

  TString outFileName = "./CR5.root";
  TFile* outFile = new TFile( outFileName, "RECREATE");
  outFile->cd();

  Float_t SFLep[NLep] = {};

  Float_t SFLepErr[NLep] = {};

  Float_t peakSFPre = 0.; Float_t peakSFPreErr = 0.;
  Float_t peakSFPost = 0.; Float_t peakSFPostErr = 0.;
  
  Float_t tailData = 0.; Float_t tailDataErr = 0.;
  Float_t tailMC = 0.; Float_t tailMCErr = 0.;
  Float_t tailDiLep = 0.; Float_t tailDiLepErr = 0.;

  Float_t SF = 0.; Float_t SFErr = 0.;
  
  TTree* tree[3];  
  for(int ilep = 0; ilep < NLep; ilep++){
    tree[ilep] = new TTree( lep[ilep], lep[ilep]);

    tree[ilep]->Branch( "peakSFPre", &peakSFPre);
    tree[ilep]->Branch( "peakSFPreErr", &peakSFPreErr);
    tree[ilep]->Branch( "peakSFPost", &peakSFPost);
    tree[ilep]->Branch( "peakSFPostErr", &peakSFPostErr);

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

  TH1D* peakPreh[NSamples];
  TH1D* peakSFPreh[NSamples];

  TH1D* peakPosth[NSamples];
  TH1D* peakSFPosth[NSamples];

  TH1D* tailh[NSamples];

  TH1D* SFh = new TH1D( "SF", "SF", 1, 0., 1.);

  for ( int iSample = 0; iSample < NSamples; iSample++){
    peakPreh[iSample] = new TH1D( TString("peakPre") + sample[iSample], TString("peakPre") + sample[iSample], 1, 0., 1.);
    peakSFPreh[iSample] = new TH1D( TString("peakSFPre") + sample[iSample], TString("peakSFPre") + sample[iSample], 1, 0., 1.);

    peakPosth[iSample] = new TH1D( TString("peakPost") + sample[iSample], TString("peakPost") + sample[iSample], 1, 0., 1.);
    peakSFPosth[iSample] = new TH1D( TString("peakSFPost") + sample[iSample], TString("peakSFPost") + sample[iSample], 1, 0., 1.);

    tailh[iSample] = new TH1D( TString("tail") + sample[iSample], TString("tail") + sample[iSample], 1, 0., 1.);
  }

  int N = srTree->GetEntries();

  TFile* inFile[NSamples];
  TDirectory* inBaseDir[NSamples];
  TH1D* mth[NSamples];
  for (int iSample = 0; iSample < NSamples; iSample++){
    TString inFileName = "./MakeHistos/"+sample[iSample]+".root";
    inFile[iSample]= new TFile( inFileName,"READ");
  }

  TH1D* numh = new TH1D( "num", "num", 1, 0., 1.);      
  TH1D* denh = new TH1D( "den", "den", 1, 0., 1.);
  TH1D* tmph = new TH1D( "tmp", "tmp", 1, 0., 1.);

  for (int ilep = 0; ilep < NLep; ilep++){
    for ( int iSR = 0; iSR < N; iSR++){      
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

      SFh->Clear(); SFh->Reset();
      for ( int iSample = 0; iSample < NSamples; iSample++){
	peakPreh[iSample]->Clear();   peakPreh[iSample]->Reset();
	peakSFPreh[iSample]->Clear(); peakSFPreh[iSample]->Reset();
	
	peakPosth[iSample]->Clear();   peakPosth[iSample]->Reset();
	peakSFPosth[iSample]->Clear(); peakSFPosth[iSample]->Reset();
	
	tailh[iSample]->Clear(); tailh[iSample]->Reset();
      }

      ////////////////////
      // Reading input
      ////////////////////

      for (int iSample = 0; iSample < NSamples; iSample++){
	TString histoName = lep[ilep]; histoName += "-SearchRegionPreIsoTrackVeto/Mt";
	inBaseDir[iSample]->GetObject(histoName, mth[iSample]);

	n1 = mth[0]->FindBin(50.);
	n2 = mth[0]->FindBin(79.99); 
	n3 = mth[0]->FindBin(mtCut);
	n4 = mth[0]->GetNbinsX() + 1;

	integral = mth[iSample]->IntegralAndError( n1, n2, error);
	peakPreh[iSample]->SetBinContent(1, integral);
	peakPreh[iSample]->SetBinError(1, error);
      }

      for (int iSample = 0; iSample < NSamples; iSample++){
	TString histoName = lep[ilep]; histoName += "-CR5/Mt";
	inBaseDir[iSample]->GetObject(histoName, mth[iSample]);

	n1 = mth[0]->FindBin(50.);
	n2 = mth[0]->FindBin(79.99); 
	n3 = mth[0]->FindBin(mtCut);
	n4 = mth[0]->GetNbinsX() + 1;

	integral = mth[iSample]->IntegralAndError( n1, n2, error);
	peakPosth[iSample]->SetBinContent(1, integral);
	peakPosth[iSample]->SetBinError(1, error);
	
	integral = mth[iSample]->IntegralAndError( n3, n4, error);
	tailh[iSample]->SetBinContent(1, integral);
	tailh[iSample]->SetBinError(1, error);
      }

      ///////////////////////////////////////////
      // Mt peak weight
      ///////////////////////////////////////////

      // peak SF Pre

      peakSFPreh[1]->Add(peakPreh[0]);
      peakSFPreh[1]->Add(peakPreh[3], -1.); 
      peakSFPreh[1]->Add(peakPreh[4], -1.);

      denh->Clear(); denh->Reset();
      denh->Add(peakPreh[1]);
      denh->Add(peakPreh[2]);

      peakSFPreh[1]->Divide(denh);

      peakSFPre = peakSFPreh[1]->GetBinContent(1);
      peakSFPreErr = peakSFPreh[1]->GetBinError(1);

      peakSFPreh[2]->SetBinContent( 1, peakSFPre);
      peakSFPreh[2]->SetBinError( 1, peakSFPreErr);

      // peak SF Post

      peakSFPosth[2]->Add(peakPosth[0]);
      tmph->Clear(); tmph->Reset();
      peakSFPosth[2]->Add(peakPosth[3], -1.);
      peakSFPosth[2]->Add(peakPosth[4], -1.);

      
      tmph->Add(peakPosth[1]);
      tmph->Multiply(peakSFPreh[1]);
      tmph->Add(peakPosth[2]);
      peakSFPosth[2]->Divide(tmph);
      
      peakSFPost = peakSFPosth[2]->GetBinContent(1);
      peakSFPostErr = peakSFPosth[2]->GetBinError(1);

      ///////////////////////////////////////////
      // Tail and SF
      ///////////////////////////////////////////

      tailData = tailh[0]->GetBinContent(1);
      tailDataErr = tailh[0]->GetBinError(1);
      
      tmph->Clear(); tmph->Reset();
      denh->Clear(); denh->Reset();
      denh->Add(tailh[1]);
      denh->Multiply(peakSFPreh[1]);
      tmph->Add(tailh[2]);
      tmph->Multiply(peakSFPosth[2]);
      denh->Add(tmph);
      denh->Add(tailh[3]);
      denh->Add(tailh[4]);
      tailMC = denh->GetBinContent(1);
      tailMCErr = denh->GetBinError(1);
      
      tailDiLep = tailh[2]->GetBinContent(1);
      tailDiLepErr = tailh[2]->GetBinError(1);

      SFh->Add(tailh[0]);
      SFh->Divide(denh);
      SF = SFh->GetBinContent(1);
      SFErr = SFh->GetBinError(1);

      ///////////////////////////////////////////////
      // Fill outTree
      ///////////////////////////////////////////////

      tree[ilep]->Fill();
      
      //////////////////////////////////
      // Pritn some output
      //////////////////////////////////

      cout<<ID<<"\t";
      cout<<iSR<<" "<<ilep<<"\t";
      cout<<peakSFPre<<"+-"<<peakSFPreErr<<"\t";
      cout<<peakSFPost<<"+-"<<peakSFPostErr<<"\t";
      cout<<SF<<"+-"<<SFErr<<endl;      
    }      
  }

  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();

  for (int ilep = 0; ilep < NLep; ilep++){
    outFile->cd();
    tree[ilep]->Write();
  }

  return 0;
}
