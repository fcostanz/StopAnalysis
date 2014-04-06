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
using namespace std;

int Rebinning( TString dirName = "Preselection/El-Preselection"){
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  const int NSamples = 11;  
  const int NSignals = 4;
  const int NControlRegions = 6;
  const int NLeps = 3;
  const int NDirs = NControlRegions * NLeps;

  TString sample[NSamples];
  sample[0] = "Data";
  sample[1] = "DiLep";
  sample[2] = "OneLep";
  sample[3] = "WJets";
  sample[4] = "Rare";
  sample[5] = "QCD";
  sample[6] = "DrellYan";
  sample[7] = "T2tb-mStop175mLSP50";
  sample[8] = "T2tb-mStop200mLSP25";
  sample[9] = "T2tb-mStop325mLSP100";
  sample[10] = "T2tb-mStop550mLSP1";

  TFile* inFile[NSamples];
  TDirectory* inDir[NSamples];
  TFile* outFile[NSamples];
  TDirectory* outDir[NSamples];

  for (int iSample = 0; iSample < NSamples; iSample++){
    TString inFileName = "./MtPeakReweighting/"+sample[iSample]+".root";
    inFile[iSample]= new TFile( inFileName,"READ");

    inDir[iSample] =  inFile[iSample]->GetDirectory(dirName);
  }
  for (int iSample = 0; iSample < NSamples; iSample++){
    TString outFileName = "./Rebinning/"+sample[iSample]+".root";
    outFile[iSample]= new TFile( outFileName,"UPDATE");

    if (outFile[iSample]->GetDirectory(dirName)) {
      string str = dirName.Data();

      unsigned pos = str.rfind("/");         // position of "live" in str

      std::string upperDir = str.substr (0,pos);
      std::string baseDir = str.substr (pos+1);
      
      cout<<upperDir<<" "<<baseDir<<endl;


      TDirectory* tempo = outFile[iSample]->GetDirectory( TString(upperDir));
      tempo->Delete(TString(baseDir)+";*");
	//outFile[iSample]->Delete(dirName + ";*")
    }
    outFile[iSample]->mkdir( dirName);
    outDir[iSample] = outFile[iSample]->GetDirectory(dirName);    
  }

  TIter next(inDir[0]->GetListOfKeys());
  TKey *key;	
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    TH1 *datah = (TH1*)key->ReadObj();	

    Int_t k=1;

    if (datah->Integral() > 0) {
      while( k < datah->GetXaxis()->GetNbins()){
	int binmin = 0;
	float min = 99999.;
	if (true){//datah->GetXaxis()->GetNbins()%k == 0){
	  TH1 *tmp = datah->Rebin(k,"tmp");
	  
	  int x1 = 0;
	  int x2 = 0;
	  
	  for (int ibin = 1; ibin < tmp->GetXaxis()->GetNbins(); ibin++){
	    if(tmp->GetBinContent(ibin) > 0.) {
	      x1 = ibin;
	      break;
	    }
	  }
	  
	  for (int ibin = tmp->GetXaxis()->GetNbins(); ibin > 1; ibin--){
	    if(tmp->GetBinContent(ibin) > 0.) {
	      x2 = ibin;
	      break;
	    }
	  }
	  
	  for (int ibin = x1+1; ibin < x2; ibin++) {
	    if(min > tmp->GetBinContent(ibin)){
	      binmin = ibin;
	      min = tmp->GetBinContent(ibin);
	    }
	  }
	  if (min > 0) break;
	}
	k++;
      }
    }
      
    outDir[0]->cd();    
    TH1F *rebDatah = (TH1F*)datah->Rebin(k, datah->GetName());
    int nbins = datah->GetXaxis()->GetNbins();
    rebDatah->SetBinContent(1, datah->GetBinContent(0) + datah->GetBinContent(1));
    rebDatah->SetBinContent(nbins, datah->GetBinContent(nbins) + datah->GetBinContent(nbins+1));
  
    rebDatah->Write();
    
    for (int iSample = 1; iSample < NSamples; iSample++){
      TH1 *mch;
      inDir[iSample]->GetObject(datah->GetName(), mch);

      outDir[iSample]->cd();   
      TH1F *rebMch = (TH1F*)mch->Rebin(k, mch->GetName());
      rebMch->SetBinContent(1, mch->GetBinContent(0) + mch->GetBinContent(1));
      rebMch->SetBinContent(nbins, mch->GetBinContent(nbins) + mch->GetBinContent(nbins+1));

      rebMch->Write();
    }

  }

  for (int iSample = 0; iSample < NSamples; iSample++){
    inFile[iSample]->Close();
    //outFile[iSample]->Write();    
    outFile[iSample]->Close();
  }

  return 0;
}
