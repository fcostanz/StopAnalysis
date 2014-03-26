#include "TROOT.h"
#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
#include "TCollection.h"
#include "TKey.h"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

bool pcp = true;
/*
Double_t On(Double_t *x, Double_t *par)
{
  Float_t mStop = 175.;
  if (( x[0] - x[1]) > mStop) return 1.;
  else return 0.;
}

Double_t Off(Double_t *x, Double_t *par)
{
  Float_t mStop = 175.;
  if (( x[0] - x[1]) <= mStop) return 1.;
  else return 0.;  
}
*/
class SR{
public:
  SR(Float_t BR_ = 0.5){
    BR = BR_;
    
    id_SR = 0;
    Ntry = -1;    
    exclusionID = -1;

    mtCut = 0.;
    nJetCut = 0.;
    topnessCut = 0.;
    mt2wCut = 0.;
    yCut = 0.;
    dphiCut = 0.;
    drlblCut = 0.;
    drlbgCut = 0.;
    chi2Cut = 0.;
    metCut = 0.;
    m3Cut = 0.;
    centralityCut = 0.;
    mlbCut = 0.;

    tth = 0;
    tbh = 0;  
    bbh = 0;
    sig_toth = 0;

    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkg[ibkg] = 0.;
      bkg2[ibkg] = 0.;
    }

    sys[0] = 0.15;
    sys[1] = 0.30;
    sys[2] = 0.30;
    sys[3] = 0.50;

    Exclusion();
    ExclusionDown();
    ExclusionUp();
  };
  ~SR(){};

  SR(const SR& copy){
    BR = copy.BR;
    
    id_SR = copy.id_SR;
    Ntry = copy.Ntry;
    exclusionID = copy.exclusionID;

    mtCut = copy.mtCut;
    nJetCut = copy.nJetCut;
    topnessCut = copy.topnessCut;
    mt2wCut = copy.mt2wCut;
    yCut = copy.yCut;
    dphiCut = copy.dphiCut;
    drlblCut = copy.drlblCut;
    drlbgCut = copy.drlbgCut;
    chi2Cut = copy.chi2Cut;
    metCut = copy.metCut;
    m3Cut = copy.m3Cut;
    centralityCut = copy.centralityCut;
    mlbCut = copy.mlbCut;

    tth = new TH2F(*copy.tth);
    tbh = new TH2F(*copy.tbh);
    bbh = new TH2F(*copy.bbh);
    sig_toth = new TH2F(*copy.sig_toth);

    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkg[ibkg] = copy.bkg[ibkg];
      bkg2[ibkg] = copy.bkg2[ibkg];
      sys[ibkg] = copy.sys[ibkg];
    }
    
    Exclusion();
    ExclusionDown();
    ExclusionUp();
  };

  SR& operator= (const SR& copy){
    BR = copy.BR;
    
    id_SR = copy.id_SR;
    Ntry = copy.Ntry;
    exclusionID = copy.exclusionID;

    mtCut = copy.mtCut;
    nJetCut = copy.nJetCut;
    topnessCut = copy.topnessCut;
    mt2wCut = copy.mt2wCut;
    yCut = copy.yCut;
    dphiCut = copy.dphiCut;
    drlblCut = copy.drlblCut;
    drlbgCut = copy.drlbgCut;
    chi2Cut = copy.chi2Cut;
    metCut = copy.metCut;
    m3Cut = copy.m3Cut;
    centralityCut = copy.centralityCut;
    mlbCut = copy.mlbCut;

    tth = new TH2F(*copy.tth);
    tbh = new TH2F(*copy.tbh);
    bbh = new TH2F(*copy.bbh);
    sig_toth = new TH2F(*copy.sig_toth);

    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkg[ibkg] = copy.bkg[ibkg];
      bkg2[ibkg] = copy.bkg2[ibkg];
      sys[ibkg] = copy.sys[ibkg];
    }
    
    return *this;
  };

  Float_t BR;

  Int_t id_SR;
  Int_t Ntry;
  Int_t exclusionID;

  Float_t mtCut;
  Float_t nJetCut;
  Float_t topnessCut;
  Float_t mt2wCut;
  Float_t yCut;
  Float_t dphiCut;
  Float_t drlblCut;
  Float_t drlbgCut;
  Float_t chi2Cut;
  Float_t metCut;
  Float_t m3Cut;
  Float_t centralityCut;
  Float_t mlbCut;

  TH2F* tth;
  TH2F* tbh;
  TH2F* bbh;
  TH2F* sig_toth;

  TH2F* sigh;
  TH2F* FOMh;
  TH2F* exclusionh;
  TH2F* exclusionUph;
  TH2F* exclusionDownh;

  Float_t bkg[4];
  Float_t bkg2[4];
  Float_t sys[4];
  
  TH2F* sig() {
    if (tth == 0) return 0;
    Float_t a = 4. * BR * BR;
    Float_t b = 4. * BR * (1 - BR);
    Float_t c = 4. * (1 - BR) * (1 - BR);

    sigh = new TH2F(*tth);

    sigh->SetName("sig");
    sigh->SetTitle("sig");
    sigh->Scale(a);

    sigh->Add( tbh, b);
    sigh->Add( bbh, c);

    return sigh;
  }

  Float_t bkg_tot() const { 
    Float_t tmp = 0;
    for (int ibkg = 0; ibkg < 4; ibkg++) tmp += bkg[ibkg];
    return tmp;
  }  
  
  TH2F* FOM(){
    if (tth == 0) return 0;
    sig();
    FOMh = new TH2F(*sigh);
    FOMh->SetName("FOM");
    FOMh->SetTitle("FOM");

    Float_t den = 0;
    for (int ibkg = 0; ibkg < 4; ibkg++)
      den += bkg[ibkg] * ( 1. + sys[ibkg] * sys[ibkg] * bkg[ibkg]);

    if (den!=0){
      FOMh->Scale(1./sqrt(den));
      return FOMh;
    }
    else 
      return 0;
  }

  TH2F* Exclusion(){
    if (tth == 0) return 0;
    FOM();
    exclusionh = new TH2F(*FOMh);
    exclusionh->SetName("exclusion");
    exclusionh->SetTitle("exclusion");

    for (int ibin = 0; ibin < exclusionh->GetSize(); ibin++){
      if ( exclusionh->GetBinContent( ibin) < 3.0)
	exclusionh->SetBinContent( ibin, 0);
      else 
	exclusionh->SetBinContent( ibin, 1);
    }
    
    return exclusionh;
  }

  TH2F* ExclusionUp(){
    if (tth == 0) return 0;
    FOM();
    exclusionUph = new TH2F(*FOMh);
    exclusionUph->SetName("exclusionUp");
    exclusionUph->SetTitle("exclusionUp");

    for (int ibin = 0; ibin < exclusionUph->GetSize(); ibin++){
      if ( exclusionUph->GetBinContent( ibin) + exclusionUph->GetBinError( ibin)< 3.0)
	exclusionUph->SetBinContent( ibin, 0);
      else 
	exclusionUph->SetBinContent( ibin, 1);
    }
    
    return exclusionUph;
  }
  
  TH2F* ExclusionDown(){
    if (tth == 0) return 0;
    FOM();
    exclusionDownh = new TH2F(*FOMh);
    exclusionDownh->SetName("exclusionDown");
    exclusionDownh->SetTitle("exclusionDown");

    for (int ibin = 0; ibin < exclusionDownh->GetSize(); ibin++){
      if ( exclusionDownh->GetBinContent( ibin) - exclusionDownh->GetBinError( ibin)< 3.0)
	exclusionDownh->SetBinContent( ibin, 0);
      else 
	exclusionDownh->SetBinContent( ibin, 1);
    }
    
    return exclusionDownh;
  }

  Float_t IntegralOn() const {
    TH2F* histo;
    if (exclusionID == 0 && exclusionh) histo = exclusionh;
    else if (exclusionID == 1 && exclusionUph) histo = exclusionUph;
    else if (exclusionID == -1 && exclusionDownh) histo = exclusionDownh;
    else return -1.;

    Float_t mStop = 175.;
    
    Float_t integral = 0.;
    Int_t i = 0;
    Int_t j = 0;
    Int_t k = 0;
    Float_t x = 0.;
    Float_t y = 0.;

    for (int ibin = 0; ibin < histo->GetSize(); ibin++){
      if (histo->IsBinOverflow(ibin)) continue;
      if (histo->IsBinUnderflow(ibin)) continue;

      histo->GetBinXYZ(ibin, i, j, k);
      
      x = histo->GetXaxis()->GetBinCenter(i);
      y = histo->GetYaxis()->GetBinCenter(j);

      if ( (x - y) >  mStop ) integral += histo->GetBinContent(ibin);
    }
    return integral;
  }


  Float_t IntegralOff() const {
    TH2F* histo;
    if (exclusionID == 0 && exclusionh) histo = exclusionh;
    else if (exclusionID == 1 && exclusionUph) histo = exclusionUph;
    else if (exclusionID == -1 && exclusionDownh) histo = exclusionDownh;
    else return -1.;

    Float_t mStop = 176.;
    
    Float_t integral = 0.;
    Int_t i = 0;
    Int_t j = 0;
    Int_t k = 0;
    Float_t x = 0.;
    Float_t y = 0.;

    for (int ibin = 0; ibin < histo->GetSize(); ibin++){
      if (histo->IsBinOverflow(ibin)) continue;
      if (histo->IsBinUnderflow(ibin)) continue;

      histo->GetBinXYZ(ibin, i, j, k);
      
      x = histo->GetXaxis()->GetBinCenter(i);
      y = histo->GetYaxis()->GetBinCenter(j);

      if ( (x - y) <=  mStop ) integral += histo->GetBinContent(ibin);
    }
    return integral;
  }


  
  Int_t Print() const {
    cout<<"ID = "<<Ntry<<endl;
    
    cout<<"mtCut = "<<mtCut<<endl;
    cout<<"nJetCut = "<<nJetCut<<endl;
    cout<<"topnessCut = "<<topnessCut<<endl;
    cout<<"mt2wCut = "<<mt2wCut<<endl;
    cout<<"yCut = "<<yCut<<endl;
    cout<<"dphiCut = "<<dphiCut<<endl;
    cout<<"drlblCut = "<<drlblCut<<endl;
    cout<<"drlbgCut = "<<drlbgCut<<endl;
    cout<<"chi2Cut = "<<chi2Cut<<endl;
    cout<<"metCut = "<<metCut<<endl;
    cout<<"m3Cut = "<<m3Cut<<endl;
    cout<<"centralityCut = "<<centralityCut<<endl;
    cout<<"mlbCut = "<<mlbCut<<endl;

    cout<<"; dilep="<<bkg[0]<<"+-"<<sqrt(bkg2[0]);
    cout<<"; onelep="<<bkg[1]<<"+-"<<sqrt(bkg2[1]);
    cout<<"; wjets="<<bkg[2]<<"+-"<<sqrt(bkg2[2]);
    cout<<"; rare="<<bkg[3]<<"+-"<<sqrt(bkg2[3]);
    cout<<"; integral="<<IntegralOn() + IntegralOff();
    cout<<endl;    

    return 0;
  }

  Int_t Write(TDirectory* base) const {    
    if (!base->IsFolder()) return 1;
    if (!base->IsWritable()) return 2;
    
    TDirectory* holdDir = gDirectory->GetDirectory("");

    TString dirName(""); dirName += BR;
    if (!base->GetDirectory(dirName)) base->mkdir(dirName);
    base->cd(dirName);
    dirName = ""; dirName += id_SR;
    if (!gDirectory->GetDirectory(dirName)) gDirectory->mkdir(dirName);
    gDirectory->cd(dirName);

    TH2F* sig_effh;
    TH1F* bkgh = new TH1F("bkg","bkg", 4, 0., 4.);   
    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkgh->SetBinContent(ibkg + 1, bkg[ibkg]);
      bkgh->SetBinError(ibkg + 1, sqrt(bkg2[ibkg]));    
    }
    bkgh->Write();

    if (sigh) sigh->Write();
    if (FOMh) FOMh->Write();      
    if (exclusionh)     exclusionh->Write();
    if (exclusionUph)   exclusionUph->Write();
    if (exclusionDownh) exclusionDownh->Write();
    if (sig_toth) sig_toth->Write();

    if (sig_toth && sigh) {
      sig_effh = new TH2F(*sigh);
      sig_effh->SetName("sig_eff");
      sig_effh->SetTitle("sig_eff");
      sig_effh->Divide(sig_toth);
    }

    if (sig_effh) sig_effh->Write();

    holdDir->cd();
    return 0;
  }
};

bool sortByIntegralOn (const SR* i,const SR* j) { return (i->IntegralOn() > j->IntegralOn());};
bool sortByIntegralOff (const SR* i,const SR* j) { return (i->IntegralOff() > j->IntegralOff());};

class SRCollection{
public:
  
  SRCollection(){
    SRs = new std::vector<SR*>;

    nBest = 4;
  };
  ~SRCollection(){};

  SR* at(int index){
    if (SRs == 0) return 0;
    else if ((int) SRs->size() <= index) return 0;
    else return SRs->at(index);
  }
  
  int size() const{
    if (SRs == 0) return -1;
    return (int) SRs->size();
  }

  Int_t SortOn(Float_t BR = 0.5) {
    TH2F* on = new TH2F(*SRs->at(0)->exclusionDownh);
    on->Reset();

    Float_t mStop = 175.;
    
    Int_t i = 0;
    Int_t j = 0;
    Int_t k = 0;
    Float_t x = 0.;
    Float_t y = 0.;

    for (int ibin = 0; ibin < on->GetSize(); ibin++){
      on->SetBinContent( ibin, 0.);
      
      if (on->IsBinOverflow(ibin)) continue;
      if (on->IsBinUnderflow(ibin)) continue;

      on->GetBinXYZ(ibin, i, j, k);
      
      x = on->GetXaxis()->GetBinCenter(i);
      y = on->GetYaxis()->GetBinCenter(j);

      if ( (x - y) >  mStop ) on->SetBinContent( ibin, 1.);
    }

    if (SRs == 0) return -1;

    for (int isr=0; isr < (int) SRs->size(); isr++){
      SRs->at(isr)->BR = BR;
      SRs->at(isr)->exclusionID = -1;
      SRs->at(isr)->ExclusionDown();
    }
    
    sort(SRs->begin(), SRs->end(), sortByIntegralOn);
    excluded = new TH2F(*SRs->at(0)->exclusionDownh);
    excluded->Multiply(on);

    cout<<SRs->at(0)->BR<<endl;    
    cout<<SRs->at(0)->id_SR<<" "<<SRs->at(0)->IntegralOn()<<endl;          
    
    for (int isr=1; isr < nBest && isr< (int) SRs->size(); isr++){
      int bestSR = -1;
      Float_t bestIntegral = -1.;
      
      for (int jsr=1; jsr< (int) SRs->size(); jsr++){
	TH2F* add = new TH2F(*SRs->at(jsr)->exclusionDownh);
	TH2F* mult = new TH2F(*SRs->at(jsr)->exclusionDownh);
	mult->Multiply(excluded);
	add->Add( mult, -1);
	add->Multiply(on);

	Float_t integral = add->Integral();
	if ( bestIntegral < integral){
	  bestIntegral = integral;
	  bestSR = jsr;
	}
      }
      
      SR* tmp = SRs->at(isr);
      SRs->at(isr) = SRs->at(bestSR);
      SRs->at(bestSR) = tmp;

      TH2F* add = new TH2F(*SRs->at(isr)->exclusionDownh);
      TH2F* mult = new TH2F(*SRs->at(isr)->exclusionDownh);
      mult->Multiply(excluded);
      add->Add( mult, -1); 
      add->Multiply(on);
      cout<<SRs->at(isr)->id_SR<<" "<<add->Integral()<<endl;
      excluded->Add(add);
    }
    
    return 0;
  }
  
  Int_t AddSortOff(Float_t BR = 0.5, Int_t skip = 4) {
    if (SRs == 0) return -1;
    TH2F* off = new TH2F(*SRs->at(0)->exclusionDownh);
    off->Reset();
    
    Float_t mTop = 176.;
    
    Int_t i = 0;
    Int_t j = 0;
    Int_t k = 0;
    Float_t x = 0.;
    Float_t y = 0.;

    for (int ibin = 0; ibin < off->GetSize(); ibin++){
      off->SetBinContent( ibin, 0.);
      
      if (off->IsBinOverflow(ibin)) continue;
      if (off->IsBinUnderflow(ibin)) continue;

      off->GetBinXYZ(ibin, i, j, k);
      
      x = off->GetXaxis()->GetBinCenter(i);
      y = off->GetYaxis()->GetBinCenter(j);

      if ( (x - y) <  mTop ) off->SetBinContent( ibin, 1.);
    }
    
    for (int isr=0; isr < (int) SRs->size(); isr++){
      SRs->at(isr)->BR = BR;
      SRs->at(isr)->exclusionID = -1;
      SRs->at(isr)->ExclusionDown();
    }

    sort(SRs->begin()+skip, SRs->end(), sortByIntegralOff);
    excluded = new TH2F(*SRs->at(skip)->exclusionDownh);
    excluded->Multiply(off);
    cout<<SRs->at(skip)->id_SR<<" "<<SRs->at(skip)->IntegralOff()<<" "<<excluded->Integral()<<endl;       

    for (int isr=skip+1; isr < skip + nBest && isr< (int) SRs->size(); isr++){
      int bestSR = -1;
      Float_t bestIntegral = -1.;
      
      for (int jsr=skip + 1; jsr< (int) SRs->size(); jsr++){
	TH2F* add = new TH2F(*SRs->at(jsr)->exclusionDownh);
	TH2F* mult = new TH2F(*SRs->at(jsr)->exclusionDownh);
	mult->Multiply(excluded);
	add->Add( mult, -1);
	add->Multiply(off);

	Float_t integral = add->Integral();
	if ( bestIntegral < integral){
	  bestIntegral = integral;
	  bestSR = jsr;
	}
      }
      
      SR* tmp = SRs->at(isr);
      SRs->at(isr) = SRs->at(bestSR);
      SRs->at(bestSR) = tmp;

      TH2F* add = new TH2F(*SRs->at(isr)->exclusionDownh);
      TH2F* mult = new TH2F(*SRs->at(isr)->exclusionDownh);
      mult->Multiply(excluded);
      add->Add( mult, -1); 
      add->Multiply(off);
      cout<<SRs->at(isr)->id_SR<<" "<<add->Integral()<<endl;
      excluded->Add(add);
    }
    
    return 0;
  }

  
  std::vector<SR*>* SRs;
  TH2F* excluded;

  Int_t nBest;
};

int bestFOM(){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  const Int_t NBR = 4;
  Float_t BRs[NBR] = {0.25, 0.50, 0.75, 1.0};

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  

  TString inFileName = "./Optimization-2Step.root";
  TFile* inFile = new TFile( inFileName,"READ");
  TTree* inTree = (TTree*)inFile->Get("Optimization");

  SR* sr = new SR();
  SRCollection* SRs = new SRCollection();

  inTree->SetBranchAddress( "SR", &sr->id_SR);

  inTree->SetBranchAddress( "mtCut", &sr->mtCut);
  inTree->SetBranchAddress( "njetCut", &sr->nJetCut);
  inTree->SetBranchAddress( "topnessCut", &sr->topnessCut);
  inTree->SetBranchAddress( "mt2wCut", &sr->mt2wCut);
  inTree->SetBranchAddress( "yCut", &sr->yCut);
  inTree->SetBranchAddress( "dphiCut", &sr->dphiCut);
  inTree->SetBranchAddress( "drlblCut", &sr->drlblCut);
  inTree->SetBranchAddress( "drlbgCut", &sr->drlbgCut);
  inTree->SetBranchAddress( "chi2Cut", &sr->chi2Cut);
  inTree->SetBranchAddress( "metCut", &sr->metCut);
  inTree->SetBranchAddress( "m3Cut", &sr->m3Cut);
  inTree->SetBranchAddress( "centralityCut", &sr->centralityCut);
  inTree->SetBranchAddress( "mlbCut", &sr->mlbCut);

  inTree->SetBranchAddress( "tt", &sr->tth);
  inTree->SetBranchAddress( "tb", &sr->tbh);
  inTree->SetBranchAddress( "bb", &sr->bbh);
  inTree->SetBranchAddress( "sig_tot", &sr->sig_toth);

  inTree->SetBranchAddress(  "diLep", &sr->bkg[0]);
  inTree->SetBranchAddress( "oneLep", &sr->bkg[1]);
  inTree->SetBranchAddress(  "wJets", &sr->bkg[2]);
  inTree->SetBranchAddress(   "rare", &sr->bkg[3]);
  
  inTree->SetBranchAddress(  "diLep2", &sr->bkg2[0]);
  inTree->SetBranchAddress( "oneLep2", &sr->bkg2[1]);
  inTree->SetBranchAddress(  "wJets2", &sr->bkg2[2]);
  inTree->SetBranchAddress(   "rare2", &sr->bkg2[3]);  
  ///////////////////////////////////////////////////// 

  int NSR = inTree->GetEntries();
  /////////////////////////////////////////////////////
  //  SR Loop
  ///////////////////////////////////////////////////// 
  for (int isr=0; isr<NSR; isr++){
    inTree->GetEntry(isr);
    
    sr->Ntry = isr;
    sr->BR = BRs[1];
    sr->exclusionID = -1;

    SRs->SRs->push_back( new SR(*sr));
  }

  TFile* outFile = new TFile( "Optimization-Coverage.root","RECREATE");
  SRs->nBest = 3;
  SRs->SortOn(0.5);
  SRs->nBest = 3;
  SRs->AddSortOff(0.5,3);

  for (int ibr = 0; ibr < 3; ibr++){    
    for (int isr=0; isr< 6; isr++){
      SRs->at(isr)->BR = BRs[ibr];

      SRs->at(isr)->sig();
      SRs->at(isr)->FOM();
      SRs->at(isr)->Exclusion();
      SRs->at(isr)->ExclusionUp();
      SRs->at(isr)->ExclusionDown();

      SRs->at(isr)->Print();
      SRs->at(isr)->Write(outFile);
    }
  }

  SRs->nBest = 3;
  SRs->SortOn(1.0);
  SRs->nBest = 2;
  SRs->AddSortOff(1.0, 3);
  
  for (int ibr = 3; ibr < 4; ibr++){    
    for (int isr=0; isr< 5; isr++){
      SRs->at(isr)->BR = BRs[ibr];

      SRs->at(isr)->sig();
      SRs->at(isr)->FOM();
      SRs->at(isr)->Exclusion();
      SRs->at(isr)->ExclusionUp();
      SRs->at(isr)->ExclusionDown();

      SRs->at(isr)->Print();
      SRs->at(isr)->Write(outFile);
    }
  }
  
  inFile->Close();
  outFile->Close();

  return 0;
}
  
