#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

bool pcp = true;

class SR{
public:
  SR(Float_t BR_ = 0.5){
    BR = BR_;
    
    tt = 0.;
    tb = 0.;  
    bb = 0.;

    tt2 = 0.;
    tb2 = 0.;  
    bb2 = 0.;

    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkg[ibkg] = 0.;
      bkg2[ibkg] = 0.;
    }

    sys[0] = 0.15;
    sys[1] = 0.30;
    sys[2] = 0.30;
    sys[3] = 0.50;

    Ntry = -1;    

    name = new TString("");
  };
  ~SR(){};

  SR(const SR& copy){
    BR = copy.BR;
    
    tt = copy.tt;
    tb = copy.tb;
    bb = copy.bb;

    tt2 = copy.tt2;
    tb2 = copy.tb2;
    bb2 = copy.bb2;

    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkg[ibkg] = copy.bkg[ibkg];
      bkg2[ibkg] = copy.bkg2[ibkg];
      sys[ibkg] = copy.sys[ibkg];
    }

    Ntry = copy.Ntry;

    name = new TString(copy.name->Data());
  };

  SR& operator= (const SR& copy){
    BR = copy.BR;
    
    tt = copy.tt;
    tb = copy.tb;
    bb = copy.bb;

    tt2 = copy.tt2;
    tb2 = copy.tb2;
    bb2 = copy.bb2;

    for (int ibkg = 0; ibkg < 4; ibkg++){
      bkg[ibkg] = copy.bkg[ibkg];
      bkg2[ibkg] = copy.bkg2[ibkg];
      sys[ibkg] = copy.sys[ibkg];
    }

    Ntry = copy.Ntry;

    name = new TString(copy.name->Data());

    return *this;
  };

  Float_t BR;

  Float_t tt;
  Float_t tb;
  Float_t bb;

  Float_t tt2;
  Float_t tb2;
  Float_t bb2;

  Float_t bkg[4];
  Float_t bkg2[4];
  Float_t sys[4];
  
  Int_t Ntry;
  
  TString* name;

  Float_t sig() const {
    Float_t a = 4. * BR * BR;
    Float_t b = 4. * BR * (1. - BR);
    Float_t c = 4. * (1. - BR) * (1. - BR);

    return (a * tt + b * tb + c * bb);
  }
  Float_t sigma2_sig() const {
    Float_t a = 4. * BR * BR;
    Float_t b = 4. * BR * (1. - BR);
    Float_t c = 4. * (1. - BR) * (1. - BR);
    
    return (a * a * tt2 + b * b * tb2 + c * c * bb2);
  }

  Float_t bkg_tot() const { 
    Float_t tmp = 0;
    for (int ibkg = 0; ibkg < 4; ibkg++) tmp += bkg[ibkg];
    return tmp;
  }
  

  Float_t FOM() const {
    Float_t den = 0;
    for (int ibkg = 0; ibkg < 4; ibkg++)
      den += bkg[ibkg] * ( 1. + sys[ibkg] * sys[ibkg] * bkg[ibkg]);
    
    return sig()/sqrt(den);
  }

  Float_t epsilon2_FOM() const {
    Float_t tmp = sigma2_sig() / sig() / sig();

    for (int ibkg = 0; ibkg < 1; ibkg++)
      tmp += pow( (double)FOM() / sig(), 4.) / 4. * bkg2[ibkg] * pow( 1. + 2. * sys[ibkg] * sys[ibkg] * bkg[ibkg], 2.);
    
    return tmp;
  }

  void Print() const {
    cout<<"OUT "<<Ntry<<" SRs.push_back("<<Ntry<<"); "<<name->Data()<<": "<<FOM()<<" +- " << sqrt( epsilon2_FOM() ) * FOM();
    cout<<"; sig="<<sig()<<" +- "<<sqrt(sigma2_sig())<<"; bkg="<<bkg_tot()<<"; 1Lep="<<bkg[1]/bkg_tot()*100<<"%; rare="<<bkg[3]/bkg_tot() * 100<<"%";
    cout<<endl;

    return;
  }
};

bool sortByFOM (const SR& i,const SR& j) { return (i.FOM() > j.FOM());}

int bestFOM( int iSample = 0, Float_t BR = 0.50){
  
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;

  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/Optimization/";

  const int NSamples = 15;
  
  Float_t mStop0[NSamples] = { 150., 250., 325., 175., 250., 325., 375., 225., 300., 450., 550., 200., 400., 600., 650.};  
  Float_t mLSP0[NSamples] =  {   1., 100., 175.,   1.,  75., 150.,  50.,   1., 100., 150.,  50.,   1., 200., 250.,   1.};
  TString sample = "T2tb-mStop"; sample += mStop0[iSample]; sample += "mLSP"; sample += mLSP0[iSample];

  std::vector<SR> SRs;
  Int_t NSavedFOMs = 30;

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////  

  TString signalFileName = "./"; signalFileName += sample; signalFileName +="_optimization-1Step.root";
  TFile* inFile[5];
  inFile[0] = new TFile(signalFileName,"READ");
  inFile[1] = new TFile(TString("DiLep") + "_optimization-1Step.root","READ");
  inFile[2] = new TFile(TString("OneLep") + "_optimization-1Step.root","READ");
  inFile[3] = new TFile(TString("WJets") + "_optimization-1Step.root","READ");
  inFile[4] = new TFile(TString("Rare") + "_optimization-1Step.root","READ");

  SR* sr = new SR(BR);

  TTree* tree[5];
  for( int i = 0; i < 5; i++)
    tree[i]= (TTree*)inFile[i]->Get("Optimization");

  int count = 0;
  //cout<<"here "<<count<<endl; count++;

  Float_t mtCut = 0.;
  Float_t nJetCut = 0.;
  Float_t topnessCut = 0.;
  Float_t mt2wCut = 0.;
  Float_t yCut = 0.;
  Float_t dphiCut = 0.;
  Float_t chi2Cut = 0.;
  Float_t drlblCut = 0.;
  Float_t drlbgCut = 0.;
  Float_t metCut = 0.;
  Float_t m3Cut = 0.;
  Float_t centralityCut = 0.;
  Float_t mlbCut = 0.;


  tree[0]->SetBranchAddress( "tt", &sr->tt);
  tree[0]->SetBranchAddress( "tb", &sr->tb);
  tree[0]->SetBranchAddress( "bb", &sr->bb);
  tree[0]->SetBranchAddress( "tt2", &sr->tt2);
  tree[0]->SetBranchAddress( "tb2", &sr->tb2);
  tree[0]->SetBranchAddress( "bb2", &sr->bb2);
  tree[0]->SetBranchAddress( "sr", &sr->name);

  tree[0]->SetBranchAddress( "mtCut", &mtCut);
  tree[0]->SetBranchAddress( "nJetCut", &nJetCut);
  tree[0]->SetBranchAddress( "topnessCut", &topnessCut);
  tree[0]->SetBranchAddress( "mt2wCut", &mt2wCut);
  tree[0]->SetBranchAddress( "yCut", &yCut);
  tree[0]->SetBranchAddress( "dphiCut", &dphiCut);
  tree[0]->SetBranchAddress( "drlblCut", &drlblCut);
  tree[0]->SetBranchAddress( "drlbgCut", &drlbgCut);
  tree[0]->SetBranchAddress( "chi2Cut", &chi2Cut);
  tree[0]->SetBranchAddress( "metCut", &metCut);
  tree[0]->SetBranchAddress( "m3Cut", &m3Cut);
  tree[0]->SetBranchAddress( "centralityCut", &centralityCut);
  tree[0]->SetBranchAddress( "mlbCut", &mlbCut);

  for (int ibkg = 0; ibkg < 4; ibkg++){
    tree[ibkg+1]->SetBranchAddress( "tt", &sr->bkg[ibkg]);
    tree[ibkg+1]->SetBranchAddress( "tt2", &sr->bkg2[ibkg]);
  }
  //===========================================
  if(pcp)cout<<"inputs set!"<<endl;
  int N = tree[0]->GetEntries();

  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 
  for (int ievt=0;ievt<N;++ievt){
    for ( int i = 0; i < 5; i++)
      tree[i]->GetEntry(ievt);
 
    sr->Ntry = ievt;
    
    if ( BR <= 0.5)
      if (nJetCut > 4) continue;

    if(sr->sig()<4.) continue;
    if(sr->bkg_tot()<1.8) continue;
    if(sqrt(sr->epsilon2_FOM()) > 0.4) continue;
    if(sr->FOM() - sqrt(sr->epsilon2_FOM()) < 1.) continue;
    if(sr->bkg[1]/sr->bkg_tot() > 0.35) continue;
    if(sr->bkg[3]/sr->bkg_tot() > 0.25) continue;
    
    
    if (SRs.size() < NSavedFOMs) SRs.push_back(*sr);
    else if ( sr->FOM() > SRs.at(9).FOM()) SRs.at(9) = *sr;
    sort(SRs.begin(), SRs.end(), sortByFOM);
  }  

  tree[0]->GetEntry(SRs.at(0).Ntry);

  Float_t mtCut0 = mtCut;
  Float_t nJetCut0 = nJetCut;
  Float_t topnessCut0 = topnessCut;
  Float_t mt2wCut0 = mt2wCut;
  Float_t yCut0 = yCut;
  Float_t dphiCut0 = dphiCut;
  Float_t drlblCut0 = drlblCut;
  Float_t drlbgCut0 = drlbgCut;
  Float_t chi2Cut0 = chi2Cut;
  Float_t metCut0 = metCut;
  Float_t m3Cut0 = m3Cut;
  Float_t centralityCut0 =  centralityCut;
  Float_t mlbCut0 = mlbCut;

  TString sigFileName = mainDir; sigFileName += sample; sigFileName += ".root";
  
  TFile* sigFile = new TFile(sigFileName,"READ");
  if (!sigFile->IsOpen()){
    std::cout<<"not open"<<std::endl;
  }

  Float_t weight = 1.;
  Float_t lumi = 19500.;

  Float_t globalWeight = 0.;
  Float_t triggerWeight = 0.;
  Float_t puWeight = 0.;
  Float_t isrWeight = 0.;
  Float_t topPtWeight = 0.;
  
  Float_t charginos = 0.;

  Float_t mt = 0.;
  Float_t njets = 0.;
  Float_t topness = 0.;
  Float_t mt2w = 0.;
  Float_t y = 0.;
  Float_t dphimin = 0.;
  Float_t drlb1 = 0.;
  Float_t hadChi2 = 0.;
  Float_t phiCorrMet = 0.;
  Float_t m3 = 0.;
  Float_t centrality = 0.;
  Float_t mlb = 0.;

  Float_t mStop = 0.;
  Float_t mLSP = 0.;

  Float_t lRelIso = 0.;
  Char_t searchRegionFlag = false;

  TTree* sigTree;
  sigTree= (TTree*)sigFile->Get("NoSystematic/bonsai");
  int NEvents = sigTree->GetEntries();
    
  sigTree->SetBranchAddress( "GlobalWeight", &globalWeight);
  sigTree->SetBranchAddress( "TriggerWeight", &triggerWeight);
  sigTree->SetBranchAddress( "PUWeight", &puWeight);
  sigTree->SetBranchAddress( "isrWeight", &isrWeight);
  sigTree->SetBranchAddress( "topPtWeight", &topPtWeight);

  sigTree->SetBranchAddress( "Charginos", &charginos);

  sigTree->SetBranchAddress( "mt", &mt);
  sigTree->SetBranchAddress( "njets", &njets);
  sigTree->SetBranchAddress( "topness", &topness);
  sigTree->SetBranchAddress( "mt2w", &mt2w);
  sigTree->SetBranchAddress( "y", &y);
  sigTree->SetBranchAddress( "dphimin", &dphimin);
  sigTree->SetBranchAddress( "drlb1", &drlb1);
  sigTree->SetBranchAddress( "hadChi2", &hadChi2);
  sigTree->SetBranchAddress( "phiCorrMet", &phiCorrMet);
  sigTree->SetBranchAddress( "m3", &m3);
  sigTree->SetBranchAddress( "centrality", &centrality);
  sigTree->SetBranchAddress( "mlb", &mlb);

  sigTree->SetBranchAddress("mStop",&mStop);
  sigTree->SetBranchAddress("mLSP",&mLSP);

  sigTree->SetBranchAddress( "lRelIso", &lRelIso);
  sigTree->SetBranchAddress( "searchRegionPost",&searchRegionFlag);

  Float_t bestX=0.;
  Float_t sigma_bestX=0.;
  Int_t SR_bestX = 0;

  Float_t a = 4. * BR * BR;
  Float_t b = 4. * BR * (1 - BR);
  Float_t c = 4. * (1 - BR) * (1 - BR);

  for (int isr=1;isr<SRs.size();++isr){
    Float_t ttUnc = 0.;
    Float_t tbUnc = 0.;
    Float_t bbUnc = 0.;
    
    Float_t tt2Unc = 0.;
    Float_t tb2Unc = 0.;
    Float_t bb2Unc = 0.;

    if (SRs.at(isr).FOM() < (0.7 * SRs.at(0).FOM() ) ) continue;

    tree[0]->GetEntry(SRs.at(isr).Ntry);
    for (int ievt=0;ievt<NEvents;++ievt){
      sigTree->GetEntry(ievt);
      
      weight = globalWeight * triggerWeight * puWeight * isrWeight * lumi;
      
      if ( fabs(mStop -  mStop0[iSample]) > 0.001) continue;
      if ( fabs(mLSP -  mLSP0[iSample]) > 0.001) continue;
      
      if (lRelIso > 0.1) continue;
      if (!searchRegionFlag) continue;
      
      if (mt < mtCut) continue;
      if (njets < nJetCut) continue;
      if (topness < topnessCut) continue;
      if (mt2w < mt2wCut) continue;
      if (y < yCut) continue;
      if (dphimin < dphiCut) continue; 
      if (drlb1 > drlblCut) continue; 
      if (drlb1 < drlbgCut) continue; 
      if (hadChi2 > chi2Cut) continue;
      if (phiCorrMet < metCut) continue;
      if (m3 < m3Cut) continue;
      if (centrality < centralityCut) continue;
      if (mlb > mlbCut) continue;

      if (!(mt < mtCut0) && !(njets < nJetCut0) && !(topness < topnessCut0)
	  && !(mt2w < mt2wCut0) && !(y < yCut0) && !(dphimin < dphiCut0)
	  && !(drlb1 > drlblCut0) && !(drlb1 < drlbgCut0) && !(hadChi2 > chi2Cut0)
	  && !(phiCorrMet < metCut0) && !(m3 < m3Cut0) && !(centrality < centralityCut0)
	  && !(mlb > mlbCut0) ) continue;
      
      if (charginos == 0){
	ttUnc  += weight;
	tt2Unc += weight * weight;
      }
      if (charginos == 1){
	tbUnc  += weight;
	tb2Unc += weight * weight;
      }
      if (charginos == 2){
	bbUnc  += weight;
	bb2Unc += weight * weight;
      }
    }
    
    Float_t X = (a * ttUnc + b * tbUnc + c * bbUnc) / SRs.at(isr).sig();

    if (X > bestX && X > 0.1){
      bestX = X;
      SR_bestX = isr;
      sigma_bestX = pow(a / SRs.at(isr).sig(), 2) * tt2Unc;
      sigma_bestX += pow(b / SRs.at(isr).sig(), 2) * tb2Unc;
      sigma_bestX += pow(c / SRs.at(isr).sig(), 2) * bb2Unc;
      sigma_bestX += pow( (a * ttUnc + b * tbUnc + c * bbUnc) / SRs.at(isr).sig() / SRs.at(isr).sig(), 2) * SRs.at(isr).sigma2_sig();
      sigma_bestX = sqrt(sigma_bestX);
    }          
  }
  SRs.at(0).Print();
  SRs.at(SR_bestX).Print();
  cout<<bestX<<"+-"<<sigma_bestX<<endl;
  
  return 0;
}
