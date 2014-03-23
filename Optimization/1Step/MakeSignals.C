#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

bool pcp = true;

int MakeSignals( int iSample = 0){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/home/fcostanz/Bonsai/Optimization/";

  const int NSignals = 15;    
  
  Float_t mStop0[NSignals] = { 150., 250., 325., 175., 250., 325., 375., 225., 300., 450., 550., 200., 400., 600., 650.};  
  Float_t mLSP0[NSignals] =  {   1., 100., 175.,   1.,  75., 150.,  50.,   1., 100., 150.,  50.,   1., 200., 250.,   1.};


  /////////////////////////////////////////////////////
  //  Variable Definition
  ///////////////////////////////////////////////////// 

  Float_t mStop = 0.;
  Float_t mLSP = 0.;
  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////
  TString inFileName = mainDir; inFileName += "T2tb.root";

  TFile* inFile = new TFile(inFileName,"READ");
  if (!inFile->IsOpen()){
    std::cout<<"not open"<<std::endl;
  }
  
  TTree* inTree;
  inTree= (TTree*)inFile->Get("NoSystematic/bonsai");
  //===========================================
  if(pcp)cout<<"inputs set!"<<endl;
  int N = inTree->GetEntries();  cout<<"THERE ARE "<<N<<" EVENTS IN "<<inFileName<<endl;
  
  inTree->SetBranchAddress("mStop",&mStop);
  inTree->SetBranchAddress("mLSP",&mLSP);

  /////////////////////////////////////////////////////
  //  Output Definition
  ///////////////////////////////////////////////////// 
  
  TString outFileName = "./Signals/T2tb-mStop"; outFileName += mStop0[iSample];
  outFileName += "mLSP"; outFileName += mLSP0[iSample]; outFileName += ".root";
  cout<<outFileName<<endl;
  
  TFile* outFile = new TFile( outFileName, "RECREATE"); 
  outFile->mkdir("NoSystematic");
  TDirectory* dir = outFile->GetDirectory("NoSystematic");
  dir->cd();
  
  inTree->LoadTree(0); //force 1st tree to be loaded
  TTree *outTree = inTree->GetTree()->CloneTree(0); 
  
  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 
  for (int ievt=0;ievt<N;++ievt){
    if (ievt%13454 == 0) {
      cout<<"Event number "<<ievt<<"\r"<<flush;
    }
    
    inTree->GetEntry(ievt);

    if ( fabs(mStop -  mStop0[iSample]) > 0.001) continue;
    if ( fabs(mLSP -  mLSP0[iSample]) > 0.001) continue;

    outTree->Fill();

  }

  dir->cd();
  outTree->Write();

  outFile->Close();    
  inFile->Close();
  
  return 0;
}
