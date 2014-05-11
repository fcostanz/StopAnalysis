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

int weightedAvarage( int n, TH1D** inh, TH1D* res);

int CR1(){
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

  TString outFileName = "./CR1.root";
  TFile* outFile = new TFile( outFileName, "RECREATE");
  outFile->cd();

  Float_t SFWLep[NLep] = {};
  Float_t SFWLepErr[NLep] = {};

  Float_t peakSF = 0.; Float_t peakSFErr = 0.;
  
  Float_t tailData = 0.; Float_t tailDataErr = 0.;
  Float_t tailDataW = 0.; Float_t tailDataWErr = 0.;
  Float_t tailMC = 0.; Float_t tailMCErr = 0.;
  Float_t tailW = 0.; Float_t tailWErr = 0.;

  Float_t SFW = 0.; Float_t SFWErr = 0.;
  Float_t SFAll = 0.; Float_t SFAllErr = 0.;
  Float_t SFComb = 0.; Float_t SFCombErr = 0.;

  Float_t RW = 0.; Float_t RWErr = 0.;
  Float_t RWCorr = 0.; Float_t RWCorrErr = 0.;

  Float_t RTop = 0.; Float_t RTopErr = 0.;
  Float_t RTopPes = 0.; Float_t RTopPesErr = 0.;
  Float_t RTopOpt = 0.; Float_t RTopOptErr = 0.;
  Float_t RTopComb = 0.; Float_t RTopCombErr = 0.;
  
  TTree* tree[3];  
  for(int ilep = 0; ilep < NLep; ilep++){
    tree[ilep] = new TTree( lep[ilep], lep[ilep]);

    tree[ilep]->Branch( "peakSF", &peakSF);
    tree[ilep]->Branch( "peakSFErr", &peakSFErr);

    tree[ilep]->Branch( "tailData", &tailData); 
    tree[ilep]->Branch( "tailDataErr", &tailDataErr);
    tree[ilep]->Branch( "tailDataW", &tailDataW);
    tree[ilep]->Branch( "tailDataWErr", &tailDataWErr);
    tree[ilep]->Branch( "tailMC", &tailMC);
    tree[ilep]->Branch( "tailMCErr", &tailMCErr);
    tree[ilep]->Branch( "tailW", &tailW);
    tree[ilep]->Branch( "tailWErr", &tailWErr);
    
    tree[ilep]->Branch( "SFW", &SFW);
    tree[ilep]->Branch( "SFWErr", &SFWErr);
    tree[ilep]->Branch( "SFAll", &SFAll);
    tree[ilep]->Branch( "SFAllErr", &SFAllErr);
    tree[ilep]->Branch( "SFComb", &SFComb);
    tree[ilep]->Branch( "SFCombErr", &SFCombErr);

    tree[ilep]->Branch( "RW", &RW);
    tree[ilep]->Branch( "RWErr", &RWErr);
    tree[ilep]->Branch( "RWCorr", &RWCorr);
    tree[ilep]->Branch( "RWCorrErr", &RWCorrErr);

    tree[ilep]->Branch( "RTop", &RTop);
    tree[ilep]->Branch( "RTopErr", &RTopErr);  
    tree[ilep]->Branch( "RTopPes", &RTopPes);
    tree[ilep]->Branch( "RTopPesErr", &RTopPesErr);
    tree[ilep]->Branch( "RTopOpt", &RTopOpt);
    tree[ilep]->Branch( "RTopOptErr", &RTopOptErr); 
    tree[ilep]->Branch( "RTopComb", &RTopComb);
    tree[ilep]->Branch( "RTopCombErr", &RTopCombErr);
  }
 
  TString histoName = "";
  TH1D* peakh[NSamples];
  TH1D* tailh[NSamples];
  TH1D* peakSFh[NSamples];
  for ( int iSample = 0; iSample < NSamples; iSample++){
    peakh[iSample] = new TH1D( TString("peak") + sample[iSample], TString("peak") + sample[iSample], 1, 0., 1.);
    tailh[iSample] = new TH1D( TString("tail") + sample[iSample], TString("tail") + sample[iSample], 1, 0., 1.);
    peakSFh[iSample] = new TH1D( TString("peakSF") + sample[iSample], TString("peakSF") + sample[iSample], 1, 0., 1.);
  }

  TH1D* RWh = new TH1D( "RW", "RW", 1, 0., 1.);
  TH1D* RToph = new TH1D( "RTop", "RTop", 1, 0., 1.);
  
  TH1D* SFh[3];
  for ( int isf = 0; isf < 3; isf++){
    histoName = "SF"; histoName += isf;
    SFh[isf] = new TH1D( histoName, histoName, 1, 0., 1.);
  }
  
  TH1D* SFCombh[NLep];
  for ( int ilep = 0; ilep < NLep; ilep++){
    histoName = "SFComb"; histoName += ilep;
    SFCombh[ilep] = new TH1D( histoName, histoName, 1, 0., 1.);
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
    SFCombh[ilep]->Clear();
    
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

      for ( int iSample = 0; iSample < NSamples; iSample++){
	peakh[iSample]->Clear();   peakh[iSample]->Reset();
	tailh[iSample]->Clear();   tailh[iSample]->Reset();
	peakSFh[iSample]->Clear(); peakSFh[iSample]->Reset();
      }

      RWh->Clear(); RWh->Reset();
      RToph->Clear(); RToph->Reset();
      
      for ( int isf = 0; isf < 3; isf++){
	SFh[isf]->Clear(); 
	SFh[isf]->Reset();
      }

      ////////////////////
      // Reading input
      ////////////////////

      for (int iSample = 0; iSample < NSamples; iSample++){
	TString histoName = lep[ilep]; histoName += "-CR1/Mt";
	inBaseDir[iSample]->GetObject(histoName, mth[iSample]);

	n1 = mth[0]->FindBin(50.);
	n2 = mth[0]->FindBin(79.99); 
	n3 = mth[0]->FindBin(mtCut);
	n4 = mth[0]->GetNbinsX() + 1;


	integral = mth[iSample]->IntegralAndError( n1, n2, error);
	peakh[iSample]->SetBinContent(1, integral);
	peakh[iSample]->SetBinError(1, error);
	
	integral = mth[iSample]->IntegralAndError( n3, n4, error);
	tailh[iSample]->SetBinContent(1, integral);
	tailh[iSample]->SetBinError(1, error);
      }

      ///////////////////////////////////////////
      // Mt peak weight
      ///////////////////////////////////////////

      peakSFh[3]->Add(peakh[0]);
      //peakSFh[3]->Add(peakh[1], -1.); 
      //peakSFh[3]->Add(peakh[2], -1.);
      peakSFh[3]->Add(peakh[4], -1.); 


      tmph->Clear(); tmph->Reset();
      tmph->Add(peakh[1]);
      tmph->Add(peakh[2]);
      tmph->Add(peakh[3]);
      //peakSFh[3]->Divide(peakh[3]);
      peakSFh[3]->Divide(tmph);

      peakSF = peakSFh[3]->GetBinContent(1);
      peakSFErr = peakSFh[3]->GetBinError(1);

      ///////////////////////////////////////////
      // integral in the tail
      ///////////////////////////////////////////

      tailData = tailh[0]->GetBinContent(1);
      tailDataErr = tailh[0]->GetBinError(1);
      
      tmph->Clear(); tmph->Reset();
      tmph->Add(tailh[0]);
      tmph->Add(tailh[1], -1.);
      tmph->Add(tailh[2], -1.);
      tmph->Add(tailh[4], -1.);
      tailDataW = tmph->GetBinContent(1);
      tailDataWErr = tmph->GetBinError(1);
      
      tmph->Clear(); tmph->Reset();
      tmph->Add(tailh[3]);
      tmph->Multiply(peakSFh[3]);
      tmph->Add(tailh[1]);
      tmph->Add(tailh[2]);
      tmph->Add(tailh[4]);
      tailMC = tmph->GetBinContent(1);
      tailMCErr = tmph->GetBinError(1);
      
      tailW = tailh[3]->GetBinContent(1);
      tailWErr = tailh[3]->GetBinError(1);

      ///////////////////////////////////////////
      // peak-to-tail ratio
      ///////////////////////////////////////////

      RWh->Add(tailh[3]);
      denh->Clear(); denh->Reset();
      denh->Add(peakh[3]);
      denh->Multiply(peakSFh[3]);
      RWh->Divide(denh);
      
      RToph->Add(tailh[2]);
      RToph->Divide(peakh[2]);

      ///////////////////////////////////////////
      // SF W
      ///////////////////////////////////////////
      
      // SFW      
      SFh[0]->Clear(); SFh[0]->Reset();
      denh->Clear(); denh->Reset();
      
      SFh[0]->Add(tailh[0]);
      SFh[0]->Add(tailh[1], -1.);
      SFh[0]->Add(tailh[2], -1.);
      SFh[0]->Add(tailh[4], -1.);

      denh->Add(tailh[3]);
      denh->Multiply(peakSFh[3]);
      
      SFh[0]->Divide(denh);

      // SFAll
      SFh[1]->Clear(); SFh[1]->Reset();
      denh->Clear(); denh->Reset();

      SFh[1]->Add(tailh[0]);

      denh->Add(tailh[3]);
      denh->Multiply(peakSFh[3]);
      denh->Add(tailh[1]);
      denh->Add(tailh[2]);
      denh->Add(tailh[4]);
      
      SFh[1]->Divide(denh);

      // SF Avarage;      
      SFh[2]->SetBinContent(1, (SFh[0]->GetBinContent(1) + SFh[1]->GetBinContent(1)) / 2. );
      SFh[2]->SetBinError(1, (SFh[0]->GetBinError(1) + SFh[1]->GetBinError(1)) / 2. + fabs(SFh[0]->GetBinContent(1) - SFh[1]->GetBinContent(1)) / 2.);

      SFW = SFh[0]->GetBinContent(1);
      SFWErr = SFh[0]->GetBinError(1);
      SFAll = SFh[1]->GetBinContent(1);
      SFAllErr = SFh[1]->GetBinError(1);
      SFComb = SFh[2]->GetBinContent(1);
      SFCombErr = SFh[2]->GetBinError(1);

      SFWLep[ilep] += SFComb / SFCombErr / SFCombErr;
      SFWLepErr[ilep] += 1. /SFCombErr / SFCombErr;

      ///////////////////////////////////////////
      // peak-to-tail ratio
      ///////////////////////////////////////////

      RWh->Add(tailh[3]);
      denh->Clear(); denh->Reset();
      denh->Add(peakh[3]);
      denh->Multiply(peakSFh[3]);
      RWh->Divide(denh);
      
      RToph->Add(tailh[2]);
      RToph->Divide(peakh[2]);

      RW = RWh->GetBinContent(1);
      RWErr = RWh->GetBinError(1);
      RWCorr = 1.3 * RW;
      RWCorrErr = 1.3 * RWErr;
      
      RTop = RToph->GetBinContent(1);
      RTopErr = RToph->GetBinError(1);
      RTopPes = RWh->GetBinContent(1) * 1.3;
      RTopPesErr = RWh->GetBinError(1) * 1.3;
      RTopOpt = RToph->GetBinContent(1) * 1.3;
      RTopOptErr = RToph->GetBinError(1) * 1.3;
      RTopComb = (RTopOpt + RTopPes) /2.;
      RTopCombErr = fabs(RTopOpt - RTopPes) / 2.;
      RTopCombErr = sqrt(RTopCombErr * RTopCombErr + (RTopOptErr - RTopPesErr) * (RTopOptErr - RTopPesErr) / 4.);

      
      ///////////////////////////////////////////////
      // Fill outTree
      ///////////////////////////////////////////////

      tree[ilep]->Fill();
      
      //////////////////////////////////
      // Pritn some output
      //////////////////////////////////

      cout<<ID<<"\t";
      cout<<iSR<<" "<<ilep<<"\t";
      cout<<peakSF<<"+-"<<peakSFErr<<"\t";
      cout<<RW<<"+-"<<RWErr<<"\t";
      cout<<SFComb<<"+-"<<SFCombErr<<"\t";
      cout<<RTop<<"+-"<<RTopErr<<"\t";
      cout<<RTopComb<<"+-"<<RTopCombErr<<endl;    
    }      
  }

  for (int ilep = 0; ilep < NLep; ilep++){
    SFWLep[ilep] /= SFWLepErr[ilep];
    SFWLepErr[ilep] = sqrt(1. / SFWLepErr[ilep]);
    cout<<SFWLep[ilep]<<"+-"<<SFWLepErr[ilep]<<endl;
  }

  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();

  for (int ilep = 0; ilep < NLep; ilep++){
    outFile->cd();
    tree[ilep]->Write();
  }

  return 0;
}


int weightedAvarage( int n, TH1D** inh, TH1D* res){
  TH1::SetDefaultSumw2(true);  

  Float_t num = 0.;
  Float_t den = 0.;

  for ( int i = 0; i < n; i++){
    Float_t tmp = 0.;
    tmp = inh[i]->GetBinContent(1);
    tmp /= inh[i]->GetBinError(1) * inh[i]->GetBinError(1);

    num += tmp;
    den += 1. / (inh[i]->GetBinError(1) * inh[i]->GetBinError(1));
  }

  res->SetBinContent( 1, num / den);
  res->SetBinError( 1, 1. / sqrt(den)); 

  return 0;
}
