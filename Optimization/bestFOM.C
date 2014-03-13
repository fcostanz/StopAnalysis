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

int bestFOM( int iSample = 0){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  const int NSamples = 7;  
  TString sample[NSamples];
  sample[0] = "T2tb-mStop175mLSP50";
  sample[1] = "T2tb-mStop200mLSP25";
  sample[2] = "T2tb-mStop375mLSP50";
  sample[3] = "T2tb-mStop250mLSP25";
  sample[4] = "T2tb-mStop325mLSP100";
  sample[5] = "T2tb-mStop450mLSP150";
  sample[6] = "T2tb-mStop550mLSP1";

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  

  TString signalFileName = "./"; signalFileName += sample[iSample]; signalFileName +="_optimization.root";
  TFile* inFile[5];
  inFile[0] = new TFile(signalFileName,"READ");
  inFile[1] = new TFile(TString("DiLep") + "_optimization.root","READ");
  inFile[2] = new TFile(TString("OneLep") + "_optimization.root","READ");
  inFile[3] = new TFile(TString("WJets") + "_optimization.root","READ");
  inFile[4] = new TFile(TString("Rare") + "_optimization.root","READ");

  TString* SR[5];
  Float_t nEvents[5] = {};
  Float_t nEvents2[5] = {};

  TTree* tree[5];
  for( int i = 0; i < 5; i++){
    tree[i]= (TTree*)inFile[i]->Get("Optimization");
    tree[i]->SetBranchAddress( "nEvents", &nEvents[i]);
    tree[i]->SetBranchAddress( "nEvents2", &nEvents2[i]);
    SR[i] = new TString("");
    tree[i]->SetBranchAddress( "sr", &SR[i]);
  }
  //===========================================
  if(pcp)cout<<"inputs set!"<<endl;
  int N = tree[0]->GetEntries();  cout<<"THERE ARE "<<N<<" EVENTS IN "<<signalFileName<<endl;

  

  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 

  for (int ievt=0;ievt<N;++ievt){
    //if (ievt%13454 == 0) {
    //  cout<<"Event number "<<ievt<<"\r"<<flush;
    // }
    
    Float_t bkg = 0.;
    Float_t sig = 0.;
    
    Float_t sigma2_sig = 0.;
    Float_t sigma2_bkg = 0.;

    Float_t epsilon_FOM = 0.;
    Float_t epsilon2_FOM = 0.;

    for ( int i = 0; i < 5; i++)
      tree[i]->GetEntry(ievt);
    
    sig = nEvents[0];
    sigma2_sig = nEvents2[0];
 
    for ( int i = 1; i < 5; i++){
      bkg += nEvents[i];
      sigma2_bkg += nEvents2[i];
    }

    Float_t FOM = 0;
    if (bkg > 0 && sig > 3.){
      FOM = sig / sqrt ( bkg + 0.15 * 0.15 * bkg * bkg);
      
      epsilon2_FOM = sigma2_sig / sig / sig + sigma2_bkg * ( 1 + 2 * 0.15 * 0.15 * bkg ) * ( 1 + 2 * 0.15 * 0.15 * bkg ) * FOM * FOM * FOM * FOM / 4. / sig / sig / sig / sig;

      if ( FOM > 2. && sqrt(epsilon2_FOM) < .2) {
	cout<<SR[1]->Data()<<": "<<FOM<<" +- " << sqrt( epsilon2_FOM ) * FOM;
	cout<<"; sig="<<sig<<"; bkg="<<bkg<<"; 1Lep="<<nEvents[2]/bkg*100<<"; rare="<<nEvents[4]/bkg * 100<<"%";
	cout<<endl;
      }
    }
  }

  /////////////////////////////////////////////////////
  //  WriteHisto
  /////////////////////////////////////////////////////
  
  /*  baseDir->cd();

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
  */

  //cout<<endl;
  //cout<<count<<endl;
   delete[] SR;
  
  for ( int i = 0; i < 5; i++){
    inFile[i]->Close();
  }
  delete[] inFile;
 
}
