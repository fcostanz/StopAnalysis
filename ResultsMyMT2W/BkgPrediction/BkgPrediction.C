#include "TROOT.h"
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TCollection.h"
#include "TKey.h"

#include <algorithm>    // std::swap

//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

bool pcp = true;

class BkgEstimation{

public:
  BkgEstimation(){
    this->Init();
  };

  ~BkgEstimation(){
    for (int ibkg = 0; ibkg < 4; ibkg++){
      delete peakPreBkgCorrh[ibkg];
      delete peakPostBkgCorrh[ibkg];
      delete srBkgCorrh[ibkg];
    }
    delete peakPreAllBkgh;
    delete peakPreAllBkgCorrh;
    delete peakPostAllBkgh;
    delete peakPostAllBkgCorrh;
    delete srAllBkgh;
    delete srAllBkgCorrh;
    
    delete SFPreh;
    delete SFPosth;

    delete RWh;
    delete RToph;

    delete tmph;
  };

  BkgEstimation(const BkgEstimation& copy){
    this->Init();

    peakPreDatah = new TH1D( *copy.peakPreDatah);
    peakPostDatah = new TH1D( *copy.peakPostDatah);
    srDatah = new TH1D( *copy.srDatah);

    for (int ibkg = 0; ibkg < 4; ibkg++){
      peakPreBkgh[ibkg] = new TH1D( *copy.peakPreBkgh[ibkg]);
      peakPostBkgh[ibkg] = new TH1D( *copy.peakPostBkgh[ibkg]);
      srBkgh[ibkg] = new TH1D( *copy.srBkgh[ibkg]);
    }

    RWh = new TH1D( *copy.RWh);
    RToph = new TH1D( *copy.RToph);

    CR4CR5 = copy.CR4CR5;
    K = copy.K;    
  };

  BkgEstimation& operator= (const BkgEstimation& copy){
    peakPreDatah = new TH1D( *copy.peakPreDatah);
    peakPostDatah = new TH1D( *copy.peakPostDatah);
    srDatah = new TH1D( *copy.srDatah);

    for (int ibkg = 0; ibkg < 4; ibkg++){
      peakPreBkgh[ibkg] = new TH1D( *copy.peakPreBkgh[ibkg]);
      peakPostBkgh[ibkg] = new TH1D( *copy.peakPostBkgh[ibkg]);
      srBkgh[ibkg] = new TH1D( *copy.srBkgh[ibkg]);
    }

    RWh = new TH1D( *copy.RWh);
    RToph = new TH1D( *copy.RToph);
    
    CR4CR5 = copy.CR4CR5;
    K = copy.K;    

    return *this;
  };

  Int_t Init(){
    TString histoname = "";
    for (int ibkg = 0; ibkg < 4; ibkg++){
      histoname = "peakPreBkg"; histoname += ibkg; histoname += "Corr";
      peakPreBkgCorrh[ibkg] = new TH1D( histoname, histoname, 1, 0., 1.);

      histoname = "peakPostBkg"; histoname += ibkg; histoname += "Corr";
      peakPostBkgCorrh[ibkg] = new TH1D( histoname, histoname, 1, 0., 1.);

      histoname = "srBkg"; histoname += ibkg; histoname += "Corr";
      srBkgCorrh[ibkg] = new TH1D( histoname, histoname, 1, 0., 1.);
    }

    peakPreAllBkgh = new TH1D( "peakPreAllBkg", "peakPreAllBkg", 1, 0., 1.);
    peakPreAllBkgCorrh = new TH1D( "peakPreAllBkgCorr", "peakPreAllBkgCorr", 1, 0., 1.);
    
    peakPostAllBkgh = new TH1D( "peakPostAllBkg", "peakPostAllBkg", 1, 0., 1.);
    peakPostAllBkgCorrh = new TH1D( "peakPostAllBkgCorr", "peakPostAllBkgCorr", 1, 0., 1.);
    
    srAllBkgh = new TH1D( "srAllBkg", "srAllBkg", 1, 0., 1.);
    srAllBkgCorrh = new TH1D( "srAllBkgCorr", "srAllBkgCorr", 1, 0., 1.);

    SFPreh = new TH1D( "SFPre", "SFPre", 1, 0., 1.);      
    SFPosth = new TH1D( "SFPost", "SFPost", 1, 0., 1.);
    
    tmph = new TH1D( "tmp", "tmp", 1, 0., 1.);
   
    CR4CR5 = 1.;
    K = 1.;
   
    return 0;
  };


  Int_t SFPre(){
    tmph->Clear(); tmph->Reset();
    peakPreAllBkgh->Clear(); peakPreAllBkgh->Reset();
    peakPreAllBkgCorrh->Clear(); peakPreAllBkgCorrh->Reset();
    SFPreh->Clear(); SFPreh->Reset();

    for (int ibkg = 0; ibkg < 4; ibkg++){
      peakPreBkgCorrh[ibkg]->Clear();
      peakPreBkgCorrh[ibkg]->Reset();
      peakPreAllBkgh->Add(peakPreBkgh[ibkg]);
    }

    SFPreh->Add(peakPreDatah);
    SFPreh->Add(peakPreBkgh[2], -1.);
    SFPreh->Add(peakPreBkgh[3], -1.);
    
    tmph->Add(peakPreBkgh[0]);
    tmph->Add(peakPreBkgh[1]);
    //tmph->Add(peakPreBkgh[2]);
  
    SFPreh->Divide(tmph);

    //SFPreh->SetBinContent(1,1.);

    for (int ibkg = 0; ibkg < 3; ibkg++){
      peakPreBkgCorrh[ibkg]->Add(peakPreBkgh[ibkg]);
      peakPreBkgCorrh[ibkg]->Multiply(SFPreh);
    }
    peakPreBkgCorrh[3]->Add(peakPreBkgh[3]);
    for (int ibkg = 0; ibkg < 4; ibkg++)
      peakPreAllBkgCorrh->Add(peakPreBkgCorrh[ibkg]);
    
    return 0;
  };

  Int_t SFPost(){
    this->SFPre();

    tmph->Clear(); tmph->Reset();
    peakPostAllBkgh->Clear(); peakPostAllBkgh->Reset();
    peakPostAllBkgCorrh->Clear(); peakPostAllBkgCorrh->Reset();
    SFPosth->Clear(); SFPosth->Reset();
    
    for (int ibkg = 0; ibkg < 4; ibkg++){
      peakPostBkgCorrh[ibkg]->Clear();
      peakPostBkgCorrh[ibkg]->Reset();
      peakPostAllBkgh->Add(peakPostBkgh[ibkg]);
    }

    SFPosth->Add(peakPostDatah);
    SFPosth->Add(peakPostBkgh[2], -1.);
    SFPosth->Add(peakPostBkgh[3], -1.);
    
    tmph->Clear(); tmph->Reset();
    tmph->Add(peakPostBkgh[0]);
    tmph->Multiply(SFPreh);
    tmph->Add(peakPostBkgh[1]);
    //tmph->Add(peakPostBkgh[2]);
    SFPosth->Divide(tmph);

    peakPostBkgCorrh[0]->Add(peakPostBkgh[0]);
    peakPostBkgCorrh[0]->Multiply(SFPreh);
    for (int ibkg = 1; ibkg < 3; ibkg++){
      peakPostBkgCorrh[ibkg]->Add(peakPostBkgh[ibkg]);
      peakPostBkgCorrh[ibkg]->Multiply(SFPosth);
    }
    peakPostBkgCorrh[3]->Add(peakPostBkgh[3]);
    for (int ibkg = 0; ibkg < 4; ibkg++)
      peakPostAllBkgCorrh->Add(peakPostBkgCorrh[ibkg]);

    return 0;
  };

  int BkgAndError( Double_t& bkg, Double_t& error){
    this->SFPost();

    tmph->Clear(); tmph->Reset();

    for (int ibkg = 0; ibkg < 4; ibkg++){
      srBkgCorrh[ibkg]->Clear();
      srBkgCorrh[ibkg]->Reset();
    }
    srAllBkgh->Clear(); srAllBkgh->Reset();
    srAllBkgCorrh->Clear(); srAllBkgCorrh->Reset();

    for (int ibkg = 0; ibkg < 4; ibkg++)
      srAllBkgh->Add(srBkgh[ibkg]);

    srBkgCorrh[0]->Add(srBkgh[0]);
    srBkgCorrh[0]->Multiply(SFPreh);

    srBkgCorrh[1]->Add(peakPostBkgh[1]);      
    srBkgCorrh[1]->Scale(RToph->GetBinContent(1));
    srBkgCorrh[1]->Scale(SFPosth->GetBinContent(1));
    
    srBkgCorrh[2]->Add(peakPostBkgh[2]);      
    srBkgCorrh[2]->Scale(RWh->GetBinContent(1));
    srBkgCorrh[2]->Scale(SFPosth->GetBinContent(1));

    srBkgCorrh[3]->Add(srBkgh[3]);
    
    for (int ibkg = 0; ibkg < 4; ibkg++)
      srAllBkgCorrh->Add(srBkgCorrh[ibkg]);

    bkg = srAllBkgCorrh->GetBinContent(1.);
    error = srAllBkgCorrh->GetBinError(1.);

    return 0;
  }
  
  int Scale( Double_t SF, Int_t ibkg, Double_t& bkg, Double_t& error){
    if ( ibkg > 3) return 1;

    peakPreBkgh[ibkg]->Scale(SF);
    peakPostBkgh[ibkg]->Scale(SF);
    srBkgh[ibkg]->Scale(SF);

    this->BkgAndError( bkg, error);

    return 0;
  };

  //////////////////////////////////
  // Stat Uncertainties
  //////////////////////////////////
  
  // Mt Peak  

  int MtPeakStat( Float_t& unc){    
    Double_t diLep = srBkgh[0]->GetBinContent(1);
    Double_t oneLep = srBkgh[1]->GetBinContent(1);
    Double_t w = srBkgh[2]->GetBinContent(1);
    
    Double_t sfPreErr = SFPreh->GetBinError(1);
    Double_t sfPostErr = SFPosth->GetBinError(1);

    unc = sqrt( diLep * diLep * sfPreErr * sfPreErr +
		oneLep * oneLep * sfPostErr * sfPostErr +
		w * w * sfPostErr * sfPostErr );
    
    return 0;
  };

  int MtPeakStat( Float_t& unc, Float_t& percent){
    this->MtPeakStat( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };

  // DiLep

  int DiLepStat( Float_t& unc){    
    unc = srBkgCorrh[0]->GetBinError(1);    

    return 0;
  };

  int DiLepStat( Float_t& unc, Float_t& percent){
    this->DiLepStat( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };
  
  int OneLepStat( Float_t& unc){    
    unc = srBkgCorrh[1]->GetBinError(1);    

    return 0;
  };

  // OneLep

  int OneLepStat( Float_t& unc, Float_t& percent){
    this->OneLepStat( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };

  int WStat( Float_t& unc){    
    unc = srBkgCorrh[2]->GetBinError(1);    

    return 0;
  };

  // W

  int WStat( Float_t& unc, Float_t& percent){
    this->WStat( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };

  //////////////////////////////////
  // Sys Uncertainties
  //////////////////////////////////

  // DiLep

  int CR4CR5Unc( Float_t& unc){
    unc = srBkgCorrh[0]->GetBinContent(1) * CR4CR5;
    
    return 0;
  };

  int CR4CR5Unc( Float_t& unc, Float_t& percent){
    this->CR4CR5Unc( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };
  
  int NJetsModellingUnc( const Float_t& K_, Float_t& unc){
    K = K_;

    Double_t bCen;
    Double_t bCenErr;
    this->BkgAndError( bCen, bCenErr);    
    
    Double_t b = 0.;
    Double_t e = 0.;
    BkgEstimation* NJetsModelling = new BkgEstimation( *this);
    NJetsModelling->Scale( K, 1, b, e);

    unc = fabs( bCen - b)/ 2.;

    return 0;
  };

  int NJetsModellingUnc( const Float_t& K_, Float_t& unc, Float_t& percent){
    this->NJetsModellingUnc( K_, unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };

  // OneLep

  int RTopUnc( Float_t& unc){    
    unc = peakPostBkgh[1]->GetBinContent(1) * SFPosth->GetBinContent(1) * RToph->GetBinError(1.);

    return 0;
  };

  int RTopUnc( Float_t& unc, Float_t& percent){
    this->RTopUnc( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };

  // W

  int WSFUnc( Float_t& unc){    
    unc = peakPostBkgh[2]->GetBinContent(1) * SFPosth->GetBinContent(1) * RWh->GetBinError(1.) * 0.3;
    unc += peakPostBkgh[1]->GetBinContent(1) * SFPosth->GetBinContent(1) * RToph->GetBinError(1.) * 0.3;
    
    return 0;
  };

  int WSFUnc( Float_t& unc, Float_t& percent){
    this->WSFUnc( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };
  
  int WXSUnc( Float_t& unc){    
    Double_t bCen;
    Double_t bCenErr;
    this->BkgAndError( bCen, bCenErr);    
    
    Double_t bUp = 0.;
    Double_t eUp = 0.;
    BkgEstimation* WXSUp = new BkgEstimation( *this);
    WXSUp->Scale( 1.5, 2, bUp, eUp);

    Double_t bDown = 0.;
    Double_t eDown = 0.;
    BkgEstimation* WXSDown = new BkgEstimation( *this);
    WXSDown->Scale( 0.5, 2, bDown, eDown);

    Double_t b1 = bDown;
    Double_t b2 = bCen;
    Double_t b3 = bUp;

    if (b1>b2) swap(b1,b2);
    if (b2>b3) swap(b2,b3);
    if (b1>b2) swap(b1,b2);
    
    unc = fabs(b3 - b1)/2.;

    return 0;
  };

  int WXSUnc( Float_t& unc, Float_t& percent){
    this->WXSUnc( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };  

  // Rare
      
  int RareXSUnc( Float_t& unc){    
    Double_t bCen;
    Double_t bCenErr;
    this->BkgAndError( bCen, bCenErr);    
    
    Double_t bUp = 0.;
    Double_t eUp = 0.;
    BkgEstimation* RareXSUp = new BkgEstimation( *this);
    RareXSUp->Scale( 1.5, 3, bUp, eUp);

    Double_t bDown = 0.;
    Double_t eDown = 0.;
    BkgEstimation* RareXSDown = new BkgEstimation( *this);
    RareXSDown->Scale( 0.5, 3, bDown, eDown);

    Double_t b1 = bDown;
    Double_t b2 = bCen;
    Double_t b3 = bUp;

    if (b1>b2) swap(b1,b2);
    if (b2>b3) swap(b2,b3);
    if (b1>b2) swap(b1,b2);
    
    unc = fabs(b3 - b1)/2.;

    return 0;
  };

  int RareXSUnc( Float_t& unc, Float_t& percent){
    this->RareXSUnc( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };


  //////////////////////////////////
  // Tot Uncertainty
  //////////////////////////////////

  int TotUnc( Float_t& unc){
    Float_t mtPeakStat;
    Float_t diLepStat;
    Float_t oneLepStat;
    Float_t wStat;

    Float_t CR4CR5Unc;
    Float_t NJetsModellingUnc;

    Float_t RTopUnc;

    Float_t wSFUnc;
    Float_t wXSUnc;
    
    Float_t rareXSUnc;
    
    this->MtPeakStat( mtPeakStat);
    this->DiLepStat( diLepStat);
    this->OneLepStat( oneLepStat);
    this->WStat( wStat);

    this->CR4CR5Unc( CR4CR5Unc); 
    this->NJetsModellingUnc( 1.2, NJetsModellingUnc);
    
    this->RTopUnc( RTopUnc);
    
    this->WSFUnc( wSFUnc); 
    this->WXSUnc( wXSUnc);
    
    this->RareXSUnc( rareXSUnc);
    
    unc = sqrt( mtPeakStat * mtPeakStat + diLepStat * diLepStat + oneLepStat * oneLepStat + wStat * wStat +
		CR4CR5Unc * CR4CR5Unc + NJetsModellingUnc * NJetsModellingUnc + 
		RTopUnc * RTopUnc + 
		wSFUnc * wSFUnc + wXSUnc * wXSUnc +
		rareXSUnc * rareXSUnc);

    return 0;
  };

  int TotUnc( Float_t& unc, Float_t& percent){
    this->TotUnc( unc);

    Double_t bkg;
    Double_t error;

    this->BkgAndError( bkg, error);

    percent = unc / bkg * 100.;

    return 0;
  };
  
  TH1D* peakPreDatah;
  TH1D* peakPostDatah;
  TH1D* srDatah;

  TH1D* peakPreBkgh[4];
  TH1D* peakPreBkgCorrh[4];
  TH1D* peakPostBkgh[4];
  TH1D* peakPostBkgCorrh[4];
  TH1D* srBkgh[4];
  TH1D* srBkgCorrh[4];
  TH1D* peakPreAllBkgh;
  TH1D* peakPreAllBkgCorrh;
  TH1D* peakPostAllBkgh;
  TH1D* peakPostAllBkgCorrh;
  TH1D* srAllBkgh;
  TH1D* srAllBkgCorrh;  

  TH1D* SFPreh;
  TH1D* SFPosth;

  TH1D* RWh;
  TH1D* RToph;

  TH1D* tmph;

  Double_t CR4CR5;
  Double_t K;
};

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

  const int NSamples = 6;  
  const int NLep = 3;

  TString sample[NSamples];
  sample[0] = "Data";
  sample[1] = "DiLep";
  sample[2] = "SemiLep";
  sample[3] = "WJets";
  sample[4] = "Rare";
  sample[5] = "SingleTop";

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

  Float_t pre[NSamples + 1] = {};     Float_t preErr[NSamples + 1] = {};
  Float_t preCorr[NSamples + 1] = {};     Float_t preCorrErr[NSamples + 1] = {};

  Float_t post[NSamples + 1] = {};    Float_t postErr[NSamples + 1] = {};
  Float_t postCorr[NSamples + 1] = {};    Float_t postCorrErr[NSamples + 1] = {};

  Float_t SFPre = 0.;  Float_t SFPreErr = 0.;
  Float_t SFPost = 0.; Float_t SFPostErr = 0.;

  Float_t sr[NSamples + 1] = {};      Float_t srErr[NSamples + 1] = {};
  Float_t srCorr[NSamples + 1] = {};  Float_t srCorrErr[NSamples + 1] = {};
  
  Float_t MtPeakStat = 0.; Float_t MtPeakStatPercent = 0.;
  Float_t DiLepStat = 0.;  Float_t DiLepStatPercent = 0.;
  Float_t OneLepStat = 0.; Float_t OneLepStatPercent = 0.;
  Float_t WStat = 0.;      Float_t WStatPercent = 0.;

  Float_t NJetsUnc = 0.;   Float_t NJetsUncPercent = 0.;
  Float_t CR4CR5Unc = 0.;  Float_t CR4CR5UncPercent = 0.;

  Float_t RTopUnc = 0.;    Float_t RTopUncPercent = 0.;

  Float_t WXSUnc = 0.;     Float_t WXSUncPercent = 0.;
  Float_t WSFUnc = 0.;     Float_t WSFUncPercent = 0.;

  Float_t RareXSUnc = 0.;  Float_t RareXSUncPercent = 0.;

  Float_t TotUnc = 0.;     Float_t TotUncPercent = 0.;

  TTree* tree[3];  
  for(int ilep = 0; ilep < NLep; ilep++){
    tree[ilep] = new TTree( lep[ilep], lep[ilep]);

    for (int iSample = 0; iSample < NSamples; iSample++){
      tree[ilep]->Branch( TString("pre") + sample[iSample], &pre[iSample]);
      tree[ilep]->Branch( TString("pre") + sample[iSample] + TString("Err"), &preErr[iSample]);
      tree[ilep]->Branch( TString("pre") + sample[iSample] + TString("Corr"), &preCorr[iSample]);
      tree[ilep]->Branch( TString("pre") + sample[iSample] + TString("CorrErr"), &preCorrErr[iSample]);

      tree[ilep]->Branch( TString("post") + sample[iSample], &post[iSample]);
      tree[ilep]->Branch( TString("post") + sample[iSample] + TString("Err"), &postErr[iSample]);
      tree[ilep]->Branch( TString("post") + sample[iSample] + TString("Corr"), &postCorr[iSample]);
      tree[ilep]->Branch( TString("post") + sample[iSample] + TString("CorrErr"), &postCorrErr[iSample]);

      tree[ilep]->Branch( TString("sr") + sample[iSample], &sr[iSample]);
      tree[ilep]->Branch( TString("sr") + sample[iSample] + TString("Err") , &srErr[iSample]);
      tree[ilep]->Branch( TString("sr") + sample[iSample] + TString("Corr") , &srCorr[iSample]);
      tree[ilep]->Branch( TString("sr") + sample[iSample] + TString("CorrErr") , &srCorrErr[iSample]);   
    }   
    tree[ilep]->Branch( "preAllBkg", &pre[NSamples]);
    tree[ilep]->Branch( "preAllBkgErr", &preErr[NSamples]);
    tree[ilep]->Branch( "preAllBkgCorr", &preCorr[NSamples]);
    tree[ilep]->Branch( "preAllBkgCorrErr", &preCorrErr[NSamples]);

    tree[ilep]->Branch( "postAllBkg", &post[NSamples]);
    tree[ilep]->Branch( "postAllBkgErr", &postErr[NSamples]);
    tree[ilep]->Branch( "postAllBkgCorr", &postCorr[NSamples]);
    tree[ilep]->Branch( "postAllBkgCorrErr", &postCorrErr[NSamples]);

    tree[ilep]->Branch( "srAllBkg", &sr[NSamples]);
    tree[ilep]->Branch( "srAllBkg", &srErr[NSamples]);    
    tree[ilep]->Branch( "srAllBkgCorr", &srCorr[NSamples]);
    tree[ilep]->Branch( "srAllBkgCorrErr", &srCorrErr[NSamples]);

    tree[ilep]->Branch( "SFPre", &SFPre);
    tree[ilep]->Branch( "SFPreErr", &SFPreErr);
    tree[ilep]->Branch( "SFPost", &SFPost);
    tree[ilep]->Branch( "SFPostErr", &SFPostErr);

    tree[ilep]->Branch( "MtPeakStat", &MtPeakStat);
    tree[ilep]->Branch( "MtPeakStatPercent", &MtPeakStatPercent);
    tree[ilep]->Branch( "DiLepStat", &DiLepStat);
    tree[ilep]->Branch( "DiLepStatPercent", &DiLepStatPercent);
    tree[ilep]->Branch( "OneLepStat", &OneLepStat);
    tree[ilep]->Branch( "OneLepStatPercent", &OneLepStatPercent);
    tree[ilep]->Branch( "WStat", &WStat);
    tree[ilep]->Branch( "WStatPercent", &WStatPercent);

    tree[ilep]->Branch( "NJetsUnc", &NJetsUnc);
    tree[ilep]->Branch( "NJetsUncPercent", &NJetsUncPercent);
    tree[ilep]->Branch( "CR4CR5Unc", &CR4CR5Unc);
    tree[ilep]->Branch( "CR4CR5UncPercent", &CR4CR5UncPercent);

    tree[ilep]->Branch( "RTopUnc", &RTopUnc);
    tree[ilep]->Branch( "RTopUncPercent", &RTopUncPercent);

    tree[ilep]->Branch( "WXSUnc", &WXSUnc);
    tree[ilep]->Branch( "WXSUncPercent", &WXSUncPercent);
    tree[ilep]->Branch( "WSFUnc", &WSFUnc);
    tree[ilep]->Branch( "WSFUncPercent", &WSFUncPercent);

    tree[ilep]->Branch( "RareXSUnc", &RareXSUnc);
    tree[ilep]->Branch( "RareXSUncPercent", &RareXSUncPercent);

    tree[ilep]->Branch( "TotUnc", &TotUnc);
    tree[ilep]->Branch( "TotUncPercent", &TotUncPercent);
  }

  TH1D* peakPreh[NSamples];
  TH1D* peakPosth[NSamples];
  TH1D* tailh[NSamples];

  TH1D* RWh = new TH1D( "RW", "RW", 1, 0., 1.);
  TH1D* RToph = new TH1D( "RTop", "RTop", 1, 0., 1.);

  for ( int iSample = 0; iSample < NSamples; iSample++){
    peakPreh[iSample] = new TH1D( TString("peakPre") + sample[iSample], TString("peakPre") + sample[iSample], 1, 0., 1.);
    peakPosth[iSample] = new TH1D( TString("peakPost") + sample[iSample], TString("peakPost") + sample[iSample], 1, 0., 1.);
    tailh[iSample] = new TH1D( TString("tail") + sample[iSample], TString("tail") + sample[iSample], 1, 0., 1.);
  }
  
  TFile* cr1File = new TFile( "../CR/CR1.root", "READ");
  TFile* cr4File = new TFile( "../CR/CR4.root", "READ");
  TFile* cr5File = new TFile( "../CR/CR5.root", "READ");

  TTree* cr1Tree;
  TTree* cr4Tree;
  TTree* cr5Tree;

  int N = srTree->GetEntries();

  TFile* inFile[NSamples];
  TDirectory* inBaseDir[NSamples];
  TH1D* mth[NSamples];

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
    cr1File->GetObject( "ElAndMu", cr1Tree);
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

    Float_t CR5 = 0.;
    cr5Tree->SetBranchAddress( "SF", &CR5);

    for ( int iSR = 0; iSR < N; iSR++){ 
      
      cr1Tree->GetEntry(iSR);
      cr4Tree->GetEntry(iSR);
      cr5Tree->GetEntry(iSR);
  
      srTree->GetEntry(iSR);
      
      RWh->SetBinContent( 1, RW);
      RWh->SetBinError( 1, RWErr);

      RToph->SetBinContent( 1, RTop);
      RToph->SetBinError( 1, RTopErr);
      
      if ( fabs(1 - CR5) > 0.45)
	CR5 = 0.5;
      else if ( fabs(1 - CR5) > 0.3)
	CR5 = 0.35;
      else if ( fabs(1 - CR5) > 0.1)
	CR5 = 0.25;
      else CR5 = 0.15;

      Float_t CR4CR5 = sqrt( CR5 * CR5 + 0.1 * 0.1);

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
	peakPreh[iSample]->Clear();   peakPreh[iSample]->Reset();
	peakPosth[iSample]->Clear();   peakPosth[iSample]->Reset();
	tailh[iSample]->Clear(); tailh[iSample]->Reset();
      }

      ////////////////////
      // Reading input
      ////////////////////
      TString histoName;
      
      histoName = lep[ilep]; histoName += "-SearchRegionPreIsoTrackVeto/wTopPtReweight";	
      inBaseDir[1]->GetObject( histoName, wToph);
      inBaseDir[2]->GetObject( histoName, tmph);
      wToph->Add( tmph);     
      
      histoName = lep[ilep]; histoName += "-SearchRegionPreIsoTrackVeto/woTopPtReweight";	
      inBaseDir[1]->GetObject( histoName, woToph);
      inBaseDir[2]->GetObject( histoName, tmph);
      woToph->Add( tmph);

      for (int iSample = 0; iSample < NSamples; iSample++){
	histoName = lep[ilep]; histoName += "-SearchRegionPreIsoTrackVeto/Mt";
	inBaseDir[iSample]->GetObject(histoName, mth[iSample]);

	n1 = mth[0]->FindBin(50.);
	n2 = mth[0]->FindBin(79.99); 
	n3 = mth[0]->FindBin(mtCut);
	n4 = mth[0]->GetNbinsX() + 1;

	integral = mth[iSample]->IntegralAndError( n1, n2, error);
	peakPreh[iSample]->SetBinContent(1, integral);
	peakPreh[iSample]->SetBinError(1, error);
      }
      if (wToph->GetBinContent(1) > 0.){
	peakPreh[1]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
	peakPreh[2]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
      }
      peakPreh[2]->Add(peakPreh[5]);
      
      histoName = lep[ilep]; histoName += "-SearchRegionPostIsoTrackVeto/wTopPtReweight";	
      inBaseDir[1]->GetObject( histoName, wToph);
      inBaseDir[2]->GetObject( histoName, tmph);
      wToph->Add( tmph);     
      
      histoName = lep[ilep]; histoName += "-SearchRegionPostIsoTrackVeto/woTopPtReweight";	
      inBaseDir[1]->GetObject( histoName, woToph);
      inBaseDir[2]->GetObject( histoName, tmph);
      woToph->Add( tmph);
      
      for (int iSample = 0; iSample < NSamples; iSample++){
	histoName = lep[ilep]; histoName += "-SearchRegionPostIsoTrackVeto/Mt";
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

      if (wToph->GetBinContent(1) > 0.){
	peakPosth[1]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
	peakPosth[2]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));

	tailh[1]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
	tailh[2]->Scale( woToph->GetBinContent(1) / wToph->GetBinContent(1));
      }
      peakPosth[2]->Add(peakPosth[5]);
      tailh[2]->Add(tailh[5]);

      ////////////////////////////////////
      // Initializing bgk estimation
      ////////////////////////////////////
      
      BkgEstimation* bkg = new BkgEstimation();
      
      bkg->peakPreDatah = new TH1D( *peakPreh[0]);
      bkg->peakPostDatah = new TH1D( *peakPosth[0]);
      bkg->srDatah = new TH1D( *tailh[0]);

      for (int ibkg = 0; ibkg < 4; ibkg++){
	bkg->peakPreBkgh[ibkg] = new TH1D( *peakPreh[ibkg + 1]);
	bkg->peakPostBkgh[ibkg] = new TH1D( *peakPosth[ibkg + 1]);
	bkg->srBkgh[ibkg] = new TH1D( *tailh[ibkg + 1]);
      }
      
      bkg->RWh = new TH1D( *RWh);
      bkg->RToph = new TH1D( *RToph);     
      
      bkg->CR4CR5 = CR4CR5;

      Double_t b;
      Double_t e;
      bkg->BkgAndError( b, e);

      ///////////////////////////////////////////////
      // Fill outTree
      ///////////////////////////////////////////////

      pre[0] = bkg->peakPreDatah->GetBinContent(1);
      preErr[0] = bkg->peakPreDatah->GetBinError(1);
      preCorr[0] = pre[0];
      preCorrErr[0] = preErr[0];

      post[0] = bkg->peakPostDatah->GetBinContent(1);
      postErr[0] = bkg->peakPostDatah->GetBinError(1);
      postCorr[0] = post[0];
      postCorrErr[0] = postErr[0];
      
      sr[0] = bkg->srDatah->GetBinContent(1);
      srErr[0] = bkg->srDatah->GetBinError(1);      
      srCorr[0] = sr[0];
      srCorrErr[0] = srErr[0];

      for (int iSample = 1; iSample < NSamples; iSample++){
	pre[iSample] = bkg->peakPreBkgh[iSample-1]->GetBinContent(1);
	preErr[iSample] = bkg->peakPreBkgh[iSample-1]->GetBinError(1);
	preCorr[iSample] = bkg->peakPreBkgCorrh[iSample-1]->GetBinContent(1);
	preCorrErr[iSample] = bkg->peakPreBkgCorrh[iSample-1]->GetBinError(1);

	post[iSample] = bkg->peakPostBkgh[iSample-1]->GetBinContent(1);
	postErr[iSample] = bkg->peakPostBkgh[iSample-1]->GetBinError(1);
	postCorr[iSample] = bkg->peakPostBkgCorrh[iSample-1]->GetBinContent(1);
	postCorrErr[iSample] = bkg->peakPostBkgCorrh[iSample-1]->GetBinError(1);

	sr[iSample] = bkg->srBkgh[iSample-1]->GetBinContent(1);
	srErr[iSample] = bkg->srBkgh[iSample-1]->GetBinError(1);
	srCorr[iSample] = bkg->srBkgCorrh[iSample-1]->GetBinContent(1);
	srCorrErr[iSample] = bkg->srBkgCorrh[iSample-1]->GetBinError(1);
      }

      pre[NSamples] = bkg->peakPreAllBkgh->GetBinContent(1);
      preErr[NSamples] = bkg->peakPreAllBkgh->GetBinError(1);
      preCorr[NSamples] = bkg->peakPreAllBkgCorrh->GetBinContent(1);
      preCorrErr[NSamples] = bkg->peakPreAllBkgCorrh->GetBinError(1);

      post[NSamples] = bkg->peakPostAllBkgh->GetBinContent(1);
      postErr[NSamples] = bkg->peakPostAllBkgh->GetBinError(1);
      postCorr[NSamples] = bkg->peakPostAllBkgCorrh->GetBinContent(1);
      postCorrErr[NSamples] = bkg->peakPostAllBkgCorrh->GetBinError(1);

      sr[NSamples] = bkg->srAllBkgh->GetBinContent(1);
      srErr[NSamples] = bkg->srAllBkgh->GetBinError(1);
      srCorr[NSamples] = bkg->srAllBkgCorrh->GetBinContent(1);
      srCorrErr[NSamples] = bkg->srAllBkgCorrh->GetBinError(1);
      
      SFPre = bkg->SFPreh->GetBinContent(1);
      SFPreErr = bkg->SFPreh->GetBinError(1);
      
      SFPost = bkg->SFPosth->GetBinContent(1);
      SFPostErr = bkg->SFPosth->GetBinError(1);
      
      bkg->MtPeakStat( MtPeakStat, MtPeakStatPercent);
      bkg->DiLepStat( DiLepStat, DiLepStatPercent);
      bkg->OneLepStat( OneLepStat, OneLepStatPercent);
      bkg->WStat( WStat, WStatPercent);

      bkg->NJetsModellingUnc( K, NJetsUnc, NJetsUncPercent);
      bkg->CR4CR5Unc( CR4CR5Unc, CR4CR5UncPercent);

      bkg->RTopUnc( RTopUnc, RTopUncPercent);

      bkg->WXSUnc( WXSUnc, WXSUncPercent);
      bkg->WSFUnc( WSFUnc, WSFUncPercent);

      bkg->RareXSUnc( RareXSUnc, RareXSUncPercent);

      bkg->TotUnc( TotUnc, TotUncPercent);
      
      tree[ilep]->Fill();

      //////////////////////////////////
      // Pritn some output
      //////////////////////////////////
      
      cout<<ID<<"\t";
      /*
      bkg->SFPosth->Add(peakPostDatah);
      tmph->Add(peakPostBkgh[0]);
      tmph->Multiply(SFPreh);
      SFPosth->Add(tmph, -1.); 
      SFPosth->Add(peakPostBkgh[3], -1.);
   
      cout<<SFPosth->GetBinContent(1)<<"+-"<<SFPosth->GetBinError(1)<<"\t";
      
      tmph->Clear(); tmph->Reset();
      tmph->Add(peakPostBkgh[1]);
      tmph->Add(peakPostBkgh[2]);
      SFPosth->Divide(tmph);
      
      cout<<tmph->GetBinContent(1)<<"+-"<<tmph->GetBinError(1)<<"\t";
      cout<<SFPosth->GetBinContent(1)<<"+-"<<SFPosth->GetBinError(1)<<endl;
      */

      //cout<<iSR<<" "<<ilep<<"\t";
      cout<<SFPre<<"+-"<<SFPreErr<<"\t";
      cout<<SFPost<<"+-"<<SFPostErr<<"\t";
      cout<<sr[0]<<"\t";
      cout<<b<<"+-"<<TotUnc<<"("<<TotUncPercent<<")"<<"\t";
      for (int iSample = 1; iSample < NSamples-1; iSample++)
	cout<<"\t"<<srCorr[iSample]<<"+-"<<srCorrErr[iSample];
      cout<<endl;
      
      delete bkg;
    }      
  }

  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();

  for (int ilep = 0; ilep < NLep; ilep++){
    outFile->cd();
    tree[ilep]->Write();
  }

  return 0;
};
