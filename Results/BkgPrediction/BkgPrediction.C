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

int BkgPrediction(){
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

  TString outFileName = "./BkgPrediction.root";
  TFile* outFile = new TFile( outFileName, "RECREATE");
  outFile->cd();

  Float_t Pre[NSamples + 1] = {};     Float_t PreErr[NSamples + 1] = {};
  Float_t Post[NSamples + 1] = {};    Float_t PostErr[NSamples + 1] = {};

  Float_t SFPre = 0.;  Float_t SFPreErr = 0.;
  Float_t SFPost = 0.; Float_t SFPostErr = 0.;

  Float_t Tail[NSamples + 1] = {};      Float_t TailErr[NSamples + 1] = {};
  Float_t TailCorr[NSamples + 1] = {};  Float_t TailCorrErr[NSamples + 1] = {};
  
  Float_t NJetsUnc = 0.;   Float_t NJetsUncPercent = 0.;   // rerun
  Float_t CR4CR5Unc = 0.;  Float_t CR4CR5UncPercent = 0.;  // just there
  Float_t DiLepStat = 0.;  Float_t DiLepStatPercent = 0.;  // just there

  Float_t WXsUnc = 0.;     Float_t WXsUncPercent = 0.;     // rerun
  Float_t WSFUnc = 0.;     Float_t WSFUncPercent = 0.;     // just there
  Float_t WStat = 0.;      Float_t WStatPercent = 0.;      // just there

  Float_t RTopUnc = 0.;    Float_t RTopUncPercent = 0.;    // just there
  Float_t OneLepStat = 0.; Float_t OneLepStatPercent = 0.; // just there

  Float_t RareXsUnc = 0.;  Float_t RareXsUncPercent = 0.;  // rerun

  Float_t TotUnc = 0.;     Float_t TotUncPercent = 0.; // just there

  TTree* tree[3];  
  for(int ilep = 0; ilep < NLep; ilep++){
    tree[ilep] = new TTree( lep[ilep], lep[ilep]);

    for (int iSample = 0; iSample < NSamples; iSample++){
      tree[ilep]->Branch( sample[iSample] + "Pre", &Pre[iSample]);
      tree[ilep]->Branch( sample[iSample] + "PreErr", &PreErr[iSample]);

      tree[ilep]->Branch( sample[iSample] + "Post", &Post[iSample]);
      tree[ilep]->Branch( sample[iSample] + "PostErr", &PostErr[iSample]);

      tree[ilep]->Branch( TString("Tail") + sample[iSample]  , &Tail[iSample]);
      tree[ilep]->Branch( TString("Tail") + sample[iSample] + TString("Err") , &TailErr[iSample]);

      tree[ilep]->Branch( TString("Tail") + sample[iSample] + TString("Corr") , &Tail[iSample]);
      tree[ilep]->Branch( TString("Tail") + sample[iSample] + TString("Err") , &TailErr[iSample]);   
    }   
    tree[ilep]->Branch( "AllMCPre", &Pre[NSamples + 1]);
    tree[ilep]->Branch( "AllMCErr", &PreErr[NSamples + 1]);
    
    tree[ilep]->Branch( "AllMCPost", &Post[NSamples + 1]);
    tree[ilep]->Branch( "AllMCPostErr", &PostErr[NSamples + 1]);

    tree[ilep]->Branch( "TailAllMC", &Tail[NSamples + 1]);
    tree[ilep]->Branch( "TailAllMC", &TailErr[NSamples + 1]);
    
    tree[ilep]->Branch( "TailAllMCCorr", &Tail[NSamples + 1]);
    tree[ilep]->Branch( "TailAllMCCorrErr", &TailErr[NSamples + 1]);

    tree[ilep]->Branch( "SFPre", &SFPre);
    tree[ilep]->Branch( "SFPreErr", &SFPreErr);
    tree[ilep]->Branch( "SFPost", &SFPost);
    tree[ilep]->Branch( "SFPostErr", &SFPostErr);
  }

  TString histoName = "";

  TH1D* peakPreh[NSamples + 1];
  TH1D* peakPosth[NSamples + 1];
  TH1D* SFPreh = new TH1D( "SFPre", "SFPost", 1, 0., 1.);
  TH1D* SFPosth = new TH1D( "SFPost", "SFPost", 1, 0., 1.);

  TH1D* tailh[NSamples + 1];
  TH1D* tailCorrh[NSamples + 1];

  TH1D* RWh = new TH1D( "RW", "RW", 1, 0., 1.);
  TH1D* RToph = new TH1D( "RTop", "RTop", 1, 0., 1.);

  TH1D* SFh = new TH1D( "SF", "SF", 1, 0., 1.);

  for ( int iSample = 0; iSample < NSamples; iSample++){
    peakPreh[iSample] = new TH1D( TString("peakPre") + sample[iSample], TString("peakPre") + sample[iSample], 1, 0., 1.);
    peakPosth[iSample] = new TH1D( TString("peakPost") + sample[iSample], TString("peakPost") + sample[iSample], 1, 0., 1.);

    tailh[iSample] = new TH1D( TString("tail") + sample[iSample], TString("tail") + sample[iSample], 1, 0., 1.);
    tailCorrh[iSample] = new TH1D( TString("tailCorr") + sample[iSample], TString("tailCorr") + sample[iSample], 1, 0., 1.);  
  }
  peakPreh[NSamples] = new TH1D( "peakPreAllMC", "peakPreAllMC", 1, 0., 1.);
  peakPosth[NSamples] = new TH1D("peakPostAllMC", "peakPostAllMC", 1, 0., 1.);
  
  tailh[NSamples] = new TH1D( "tailAllMC", "tailAllMC", 1, 0., 1.);
  tailCorrh[NSamples] = new TH1D( "tailCorrAllMC", "tailCorrAllMC", 1, 0., 1.);  

  TFile* cr1File = new TFile( "../MakeHistos/CR1.root", "READ");
  TFile* cr4File = new TFile( "../MakeHistos/CR4.root", "READ");
  TFile* cr5File = new TFile( "../MakeHistos/CR5.root", "READ");

  TTree* cr1Tree;
  TTree* cr4Tree;
  TTree* cr5Tree;

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
    cr1File->GetObject( lep[ilep], cr1Tree);
    cr4File->GetObject( lep[ilep], cr4Tree);
    cr5File->GetObject( lep[ilep], cr5Tree);

    if ( cr1Tree->GetEntries() != N ||
	 cr4Tree->GetEntries() != N ||
	 cr5Tree->GetEntries() != N ){
      cout<<"Error: not consistent number of entries"<<endl;
    }
    
    Float_t RW = 0.;   Float_t RWErr = 0.;
    Float_t RTop = 0.; Float_t RTopErr = 0.;

    cr1Tree->SetBranchAddress( "RWCorr", &RW);
    cr1Tree->SetBranchAddress( "RWCorrErr", &RWErr);
    
    cr1Tree->SetBranchAddress( "RTopComb", &RTop);
    cr1Tree->SetBranchAddress( "RTopCombErr", &RTopErr);

    Float_t K = 0.;   Float_t KErr = 0.;
    cr4Tree->SetBranchAddress( "K", &K);
    cr4Tree->SetBranchAddress( "KErr", &KErr);

    for ( int iSR = 0; iSR < N; iSR++){    
      cr1Tree->GetEntry(iSR);
      cr4Tree->GetEntry(iSR);
      cr5Tree->GetEntry(iSR);
  
      srTree->GetEntry(iSR);
      
      RWh->SetBinContent( 1, RW);
      RWh->SetBinError( 1, RWErr);

      RToph->SetBinContent( 1, RTop);
      RToph->SetBinError( 1, RTopErr);
      

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
      for ( int iSample = 0; iSample < NSamples + 1; iSample++){
	peakPreh[iSample]->Clear();   peakPreh[iSample]->Reset();
	peakPosth[iSample]->Clear();   peakPosth[iSample]->Reset();

	tailh[iSample]->Clear(); tailh[iSample]->Reset();
	tailCorrh[iSample]->Clear(); tailCorrh[iSample]->Reset();
      }

      SFPreh->Clear(); SFPreh->Reset();
      SFPosth->Clear(); SFPosth->Reset();

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
	TString histoName = lep[ilep]; histoName += "-SearchRegionPostIsoTrackVeto/Mt";
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

      SFPreh->Add(peakPreh[0]);
      SFPreh->Add(peakPreh[4], -1.);

      denh->Clear(); denh->Reset();      
      denh->Add(peakPreh[1]);
      denh->Add(peakPreh[2]);
      denh->Add(peakPreh[3]);

      SFPreh->Divide(denh);
      
      // peak SF Post

      SFPosth->Add(peakPosth[0]);
      tmph->Clear(); tmph->Reset();
      tmph->Add(peakPosth[1]);
      tmph->Multiply(SFPreh);
      SFPosth->Add(tmph, -1.); 
      SFPosth->Add(peakPosth[4], -1.);

      denh->Clear(); denh->Reset();
      denh->Add(peakPosth[2]);
      denh->Add(peakPosth[3]);
      SFPosth->Divide(denh);
      
      ///////////////////////////////////////////
      // Tail and SF
      ///////////////////////////////////////////

      tailCorrh[0]->Add(tailh[0]);

      tailCorrh[1]->Add(tailh[1]);
      tailCorrh[1]->Multiply(SFPreh);
      
      tailCorrh[2]->Add(peakPosth[2]);      
      tailCorrh[2]->Scale(RToph->GetBinContent(1.));
      tailCorrh[2]->Multiply(SFPosth);
      
      tailCorrh[3]->Add(peakPosth[3]);
      tailCorrh[3]->Multiply(RWh);
      tailCorrh[3]->Multiply(SFPosth);

      tailCorrh[4]->Add(tailh[4]);


     for (int iSample = 1; iSample < NSamples; iSample++){
	peakPreh[NSamples]->Add(peakPreh[iSample]);
	peakPosth[NSamples]->Add(peakPosth[iSample]);
	
	tailh[NSamples]->Add(tailh[iSample]);
	tailCorrh[NSamples]->Add(tailCorrh[iSample]);
      }

     CR4CR5Unc = tailCorrh[1]->GetBinContent(1) * 0.30;
     CR4CR5UncPercent = CR4CR5Unc / tailCorrh[NSamples]->GetBinContent(1) * 100;
     DiLepStat = tailCorrh[1]->GetBinError(1);
     DiLepStatPercent = DiLepStat / tailCorrh[NSamples]->GetBinContent(1) * 100;
     

     WSFUnc = (tailCorrh[2]->GetBinContent(1) + tailCorrh[3]->GetBinContent(1)) * 0.30;
     WSFUncPercent = WSFUnc / tailCorrh[NSamples]->GetBinContent(1) * 100;
     WStat = tailCorrh[3]->GetBinError(1);
     WStatPercent = WStat / tailCorrh[NSamples]->GetBinContent(1) * 100;

     RTopUnc = tailCorrh[2]->GetBinContent(1) * RToph->GetBinError(1.) / RToph->GetBinContent(1.);
     RTopUncPercent = RTopUnc / tailCorrh[NSamples]->GetBinContent(1) * 100;
     OneLepStat = tailCorrh[2]->GetBinError(1);
     OneLepStatPercent = OneLepStat / tailCorrh[NSamples]->GetBinContent(1) * 100;
     
     TotUnc = sqrt( CR4CR5Unc * CR4CR5Unc + DiLepStat * DiLepStat + 
		    WSFUnc * WSFUnc + WStat * WStat +
		    RTopUnc * RTopUnc + OneLepStat * OneLepStat);
     
     TotUncPercent = TotUnc / tailCorrh[NSamples]->GetBinContent(1) * 100;

     
      ///////////////////////////////////////////////
      // Fill outTree
      ///////////////////////////////////////////////

      for (int iSample = 0; iSample < NSamples + 1; iSample++){
	Pre[iSample] = peakPreh[iSample]->GetBinContent(1);
	PreErr[iSample] = peakPreh[iSample]->GetBinError(1);

	Post[iSample] = peakPosth[iSample]->GetBinContent(1);
	PostErr[iSample] = peakPosth[iSample]->GetBinError(1);

	Tail[iSample] = tailh[iSample]->GetBinContent(1);
	TailErr[iSample] = tailh[iSample]->GetBinError(1);

	TailCorr[iSample] = tailCorrh[iSample]->GetBinContent(1);
	TailCorrErr[iSample] = tailCorrh[iSample]->GetBinError(1);
      }

      SFPre = SFPreh->GetBinContent(1);
      SFPreErr = SFPreh->GetBinError(1);
      
      SFPost = SFPosth->GetBinContent(1);
      SFPostErr = SFPosth->GetBinError(1);
      
      tree[ilep]->Fill();

      //////////////////////////////////
      // Pritn some output
      //////////////////////////////////

      cout<<ID<<"\t";
      cout<<iSR<<" "<<ilep<<"\t";
      cout<<SFPre<<"+-"<<SFPreErr<<"\t";
      cout<<SFPost<<"+-"<<SFPostErr<<"\t";
      cout<<Tail[0]<<"\t";
      cout<<Tail[NSamples]<<"+-"<<TailErr[NSamples]<<"\t";
      cout<<TailCorr[NSamples]<<"+-"<<TotUnc<<endl;
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



class BkgEstimation{

public:

  BkgEstimation(){};





  TH1D* peakPreData;
  TH1D* peakPostData;
  TH1D* peakPostData;

  std::vector<TH1D*> peakPreMC;
  std::vector<TH1D*> peakPostMC;
  std::vector<TH1D*> peakPostMC;

  


  TH1D* SFPreh;
  TH1D* SFPosth;

  TH1D* RWh;
  TH1D* RToph;

  TH1D* tmph;

  int Scale( Double_t SF, Int_t iBkg);
};

int BkgEstimation::GetBkgAndError(){

  TH1D* SFPreh;
  getSFPre( peakPreh, SFPreh);

  TH1D* SFPosth;
  getSFPost( peakPosth, SFPreh, SFPosth);

  getBkg( TH1D** peakPosth, TH1D** tailPosth, TH1D* SFPreh, TH1D* SFPreh, Double_t& bgk);
}
  
  
int BkgEstimation::GetSFPre( TH1D** peakPreh, TH1D* SFPreh){
  TH1D* numh = new TH1D( "numh", "numh", 1., 0., 1.);
  TH1D* denh = new TH1D( "denh", "denh", 1., 0., 1.);
  TH1D* tmph = new TH1D( "tmph", "tmph", 1., 0., 1.);

  SFPreh = new TH1D( "SFPre", "SFPre", 1., 0., 1.);
  
  SFPreh->Add(peakPreh[0]);
  SFPreh->Add(peakPreh[4], -1.);
  
  denh->Add(peakPreh[1]);
  denh->Add(peakPreh[2]);
  denh->Add(peakPreh[3]);
  
  SFPreh->Divide(denh);

  delete numh;
  delete denh;
  delete tmph;

  return 0;
}
  
int BkgEstimation::GetSFPre( TH1D** peakPosth, TH1D* SFPreh, TH1D* SFPosth){
  TH1D* numh = new TH1D( "numh", "numh", 1., 0., 1.);
  TH1D* denh = new TH1D( "denh", "denh", 1., 0., 1.);
  TH1D* tmph = new TH1D( "tmph", "tmph", 1., 0., 1.);

  SFPosth = new TH1D( "SFPost", "SFPost", 1., 0., 1.);

  SFPosth->Add(peakPosth[0]);

  tmph->Add(peakPosth[1]);
  tmph->Multiply(SFPreh);
  SFPosth->Add(tmph, -1.); 
  SFPosth->Add(peakPosth[4], -1.);
  
  denh->Add(peakPosth[2]);
  denh->Add(peakPosth[3]);
  SFPosth->Divide(denh);
  
  delete numh;
  delete denh;
  delete tmph;

  return 0;
}

int  getBkg( TH1D** peakPosth, TH1D** tailPosth, TH1D* SFPreh, TH1D* SFPosth, Double_t& bgk){
  TH1D* tailCorrh[5];
  for ( int i = 0; i < 5; i++)
    tailCorrh[i] = new TH1D( TString("") + i,  TString("") + i, 1, 0., 1.);

  tailCorrh[1]->Add(tailPosth[1]);
  tailCorrh[1]->Multiply(SFPreh);
  
  tailCorrh[2]->Add(peakPosth[2]);      
  tailCorrh[2]->Scale(RToph->GetBinContent(1.));
  tailCorrh[2]->Multiply(SFPosth);
  
  tailCorrh[3]->Add(peakPosth[3]);
  tailCorrh[3]->Multiply(RWh);
  tailCorrh[3]->Multiply(SFPosth);
  
  tailCorrh[4]->Add(tailh[4]);

  delete numh;
  delete denh;
  delete tmph;

  return 0;
}
