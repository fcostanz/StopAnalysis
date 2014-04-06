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

int MtPeakReweighting(int iSR = 6){
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  const int NSamples = 6;  
  const int NControlRegions = 6;
  const int NLep = 2;
  const int NDir = NControlRegions * NLep;

  TString sample[NSamples];
  sample[0] = "Data";
  sample[1] = "DiLep";
  sample[2] = "OneLep";
  sample[3] = "WJets";
  sample[4] = "Rare";
  sample[5] = "DrellYan";

  bool lepFlag[NLep];
  TString lep[NLep];
  lep[0] = "El";
  lep[1] = "Mu";
  
  TString controlRegion[NControlRegions];
  controlRegion[0] = "Preselection";
  controlRegion[1] = "SearchRegionPreIsoTrackVeto";
  controlRegion[2] = "SearchRegionPostIsoTrackVeto";
  controlRegion[3] = "CR1";
  controlRegion[4] = "CR4";
  controlRegion[5] = "CR5";
  
  Float_t data = 0.;
  Float_t diLep = 0.;
  Float_t oneLep = 0.;
  Float_t wJets = 0.;
  Float_t rare = 0.;
  Float_t dy = 0.;

  Float_t SF[NSamples] = {};

  Int_t n1 = 0;
  Int_t n2 = 0;

  TString baseDirName = ""; baseDirName += iSR;
  TFile* inFile[NSamples];
  TDirectory* inBaseDir[NSamples];
  TFile* outFile[NSamples];
  TDirectory* outBaseDir[NSamples];
  TH1D* mth[NSamples];
  for (int iSample = 0; iSample < NSamples; iSample++){
    TString inFileName = "./MakeHistos/"+sample[iSample]+".root";
    inFile[iSample]= new TFile( inFileName,"READ");

    inBaseDir[iSample] =  inFile[iSample]->GetDirectory(baseDirName);
  }
  for (int iSample = 1; iSample < NSamples; iSample++){
    TString outFileName = "./MtPeakReweighting/"+sample[iSample]+".root";
    outFile[iSample]= new TFile( outFileName,"UPDATE");

    if (outFile[iSample]->GetDirectory(baseDirName)) outFile[iSample]->Delete(baseDirName + ";*");
    outFile[iSample]->mkdir( baseDirName);
    outBaseDir[iSample] = outFile[iSample]->GetDirectory(baseDirName);    
  }

  for (int ilep = 0; ilep < NLep; ilep++){
    /////////////////////////////////////////////////////
    // Reweighting histos Preselection
    /////////////////////////////////////////////////////
    for (int iSample = 1; iSample < NSamples; iSample++){
      SF[iSample] = 1.;
      TDirectory* inDir =  inBaseDir[iSample]->GetDirectory(lep[ilep] + "-Preselection");
      TDirectory* outDir =  outBaseDir[iSample]->mkdir(lep[ilep] + "-Preselection");
      outDir->cd();

      TIter next(inDir->GetListOfKeys());
      TKey *key;	
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();	

	h->Scale(SF[iSample]);	
	h->Write( h->GetName());
      }
    }


    /////////////////////////////////////////////////////
    // Reweighting histos SearchRegionPreIsoTrackVeto
    /////////////////////////////////////////////////////
    for (int iSample = 0; iSample < NSamples; iSample++){
      SF[iSample] = 1.;
      TString histoName = lep[ilep]; histoName += "-SearchRegionPreIsoTrackVeto/Mt";
      inBaseDir[iSample]->GetObject(histoName, mth[iSample]);
    }

    n1 = mth[0]->FindBin(50.);
    n2 = mth[0]->FindBin(79.99); 

    data = mth[0]->Integral( n1, n2);
    diLep = mth[1]->Integral( n1, n2);
    oneLep = mth[2]->Integral( n1, n2);
    wJets = mth[3]->Integral( n1, n2);
    rare = mth[4]->Integral( n1, n2);
    dy = mth[5]->Integral( n1, n2);

    if (diLep + oneLep + wJets > 0.) SF[1] = (data - rare) / (diLep + oneLep + wJets);
    SF[2] = SF[1];
    SF[3] = SF[1];
    
    for (int iSample = 1; iSample < NSamples; iSample++){
      TDirectory* inDir =  inBaseDir[iSample]->GetDirectory(lep[ilep] + "-SearchRegionPreIsoTrackVeto");
      TDirectory* outDir =  outBaseDir[iSample]->mkdir(lep[ilep] + "-SearchRegionPreIsoTrackVeto");
      outDir->cd();

      TIter next(inDir->GetListOfKeys());
      TKey *key;	
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();	

	h->Scale(SF[iSample]);	
	h->Write( h->GetName());
      }
    }

    ////////////////////////////////////////////////////////////
    // Reweighting histos SearchRegionPostIsoTrackVeto and CR5
    ////////////////////////////////////////////////////////////
    for (int iSample = 0; iSample < NSamples; iSample++){
      //SF[iSample] = 1.;
      TString histoName = lep[ilep]; histoName += "-SearchRegionPostIsoTrackVeto/Mt";
      inBaseDir[iSample]->GetObject(histoName, mth[iSample]);
    }

    data = mth[0]->Integral( n1, n2);
    diLep = mth[1]->Integral( n1, n2);
    oneLep = mth[2]->Integral( n1, n2);
    wJets = mth[3]->Integral( n1, n2);
    rare = mth[4]->Integral( n1, n2);
    dy = mth[5]->Integral( n1, n2);

    SF[0] = 1.;
    if (oneLep + wJets > 0.) SF[2] = (data - SF[1] * diLep - rare) / (oneLep + wJets);
    SF[3] = SF[2];
    SF[4] = 1.;
    SF[5] = 1.;
    SF[6] = 1.;
 
    for (int iSample = 1; iSample < NSamples; iSample++){
      TDirectory* inDir =  inBaseDir[iSample]->GetDirectory(lep[ilep] + "-SearchRegionPostIsoTrackVeto");
      TDirectory* outDir =  outBaseDir[iSample]->mkdir(lep[ilep] + "-SearchRegionPostIsoTrackVeto");
      outDir->cd();

      TIter next(inDir->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();	

	h->Scale(SF[iSample]);	
	h->Write( h->GetName());
      }
    }
    
    for (int iSample = 0; iSample < NSamples; iSample++){
      //SF[iSample] = 1.;
      TString histoName = lep[ilep]; histoName += "-CR5/Mt";
      inBaseDir[iSample]->GetObject(histoName, mth[iSample]);
    }

    data = mth[0]->Integral( n1, n2);
    diLep = mth[1]->Integral( n1, n2);
    oneLep = mth[2]->Integral( n1, n2);
    wJets = mth[3]->Integral( n1, n2);
    rare = mth[4]->Integral( n1, n2);
    dy = mth[5]->Integral( n1, n2);

    SF[0] = 1.;
    if (oneLep + wJets > 0.) SF[2] = (data - SF[1] * diLep - rare) / (oneLep + wJets);
    SF[3] = SF[2];
    SF[4] = 1.;
    SF[5] = 1.;
    SF[6] = 1.;

    for (int iSample = 1; iSample < NSamples; iSample++){
      TDirectory* inDir =  inBaseDir[iSample]->GetDirectory(lep[ilep] + "-CR5");
      TDirectory* outDir =  outBaseDir[iSample]->mkdir(lep[ilep] + "-CR5");
      outDir->cd();

      TIter next(inDir->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();	

	h->Scale(SF[iSample]);	
	h->Write( h->GetName());
      }
    }
  
    /////////////////////////////////////////////////////
    // Reweighting histos CR1
    /////////////////////////////////////////////////////
    for (int iSample = 0; iSample < NSamples; iSample++){
      SF[iSample] = 1.;
      TString histoName = lep[ilep]; histoName += "-CR1/Mt";
      inBaseDir[iSample]->GetObject(histoName, mth[iSample]);
    }

    data = mth[0]->Integral( n1, n2);
    diLep = mth[1]->Integral( n1, n2);
    oneLep = mth[2]->Integral( n1, n2);
    wJets = mth[3]->Integral( n1, n2);
    rare = mth[4]->Integral( n1, n2);
    dy =  mth[5]->Integral( n1, n2);

    if (oneLep + wJets > 0.) SF[2] = (data - diLep - rare) / (oneLep + wJets);
    SF[3] = SF[2];

    for (int iSample = 1; iSample < NSamples; iSample++){
      TDirectory* inDir =  inBaseDir[iSample]->GetDirectory(lep[ilep] + "-CR1");
      TDirectory* outDir =  outBaseDir[iSample]->mkdir(lep[ilep] + "-CR1");
      outDir->cd();

      TIter next(inDir->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();	

	h->Scale(SF[iSample]);	
	h->Write( h->GetName());
      }
    }

    /////////////////////////////////////////////////////
    // Reweighting histos CR4
    /////////////////////////////////////////////////////
    for (int iSample = 0; iSample < NSamples; iSample++){
      SF[iSample] = 1.;
      TString histoName = lep[ilep]; histoName += "-CR4/Mt";
      inBaseDir[iSample]->GetObject(histoName, mth[iSample]);
    }

    data = mth[0]->Integral( n1, n2);
    diLep = mth[1]->Integral( n1, n2);
    oneLep = mth[2]->Integral( n1, n2);
    wJets = mth[3]->Integral( n1, n2);
    rare = mth[4]->Integral( n1, n2);
    dy = mth[5]->Integral( n1, n2);

    //SF[1] = (data - oneLep - wJets - rare - dy) / (diLep);
    //cout<<SF[1]<<endl;

    for (int iSample = 1; iSample < NSamples; iSample++){
      TDirectory* inDir =  inBaseDir[iSample]->GetDirectory(lep[ilep] + "-CR4");
      TDirectory* outDir =  outBaseDir[iSample]->mkdir(lep[ilep] + "-CR4");
      outDir->cd();

      TIter next(inDir->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *h = (TH1*)key->ReadObj();	

	h->Scale(SF[iSample]);	
	h->Write( h->GetName());
      }
    }
  }

  for (int iSample = 1; iSample < NSamples; iSample++){
    outBaseDir[iSample]->Write();    
    outFile[iSample]->Close();
  }

  for (int iSample = 1; iSample < NSamples; iSample++){
    TString outFileName = "./MtPeakReweighting/"+sample[iSample]+".root";
    outFile[iSample]= new TFile( outFileName,"UPDATE");
    outBaseDir[iSample] = outFile[iSample]->GetDirectory(baseDirName);

    for ( int iControlRegion = 0; iControlRegion < NControlRegions; iControlRegion++ ){
      TDirectory* inDirLep[NLep];
      for (int ilep = 0; ilep < NLep; ilep++)
	inDirLep[ilep] = outBaseDir[iSample]->GetDirectory(lep[ilep] + "-" + controlRegion[iControlRegion]);

      TDirectory* outDirTot = outBaseDir[iSample]->mkdir("ElAndMu-" + controlRegion[iControlRegion]);
      outDirTot->cd();

      TIter next(inDirLep[0]->GetListOfKeys());
      TKey *key;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	TH1 *hLep[NLep];
	TH1 *htot;
	htot= (TH1*)key->ReadObj();
	for (int ilep = 1; ilep < NLep; ilep++){
	  inDirLep[ilep]->GetObject(htot->GetName(), hLep[ilep]);
	  htot->Add(hLep[ilep]);
	}
	htot->Write( htot->GetName());
      }
    }
  }

  for (int iSample = 0; iSample < NSamples; iSample++)
    inFile[iSample]->Close();

  for (int iSample = 1; iSample < NSamples; iSample++)
    outFile[iSample]->Close();

  return 0;
}
