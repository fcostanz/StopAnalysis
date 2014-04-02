#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"
//#include "../src/TriggerEfficiencyProvider.cpp"
#include <iostream>
#include <iomanip>

bool pcp = true;
using namespace std;

int MakeTree( int iSample = 0){
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;
    
  TString mainDir = "/nfs/dust/cms/user/fcost/store/Bonsai/Optimization/";

  const int NBkgs = 4;
  const int NSignals = 15;
  const int NSamples = NBkgs + NSignals;
  
  TString sampleName[NBkgs + 1];
  sampleName[0] = "DiLep";
  sampleName[1] = "OneLep";
  sampleName[2] = "WJets";
  sampleName[3] = "Rare";
  sampleName[4] = "T2tb";

  Float_t mStop0[NSignals] = { 150., 250., 325., 175., 250., 325., 375., 225., 300., 450., 550., 200., 400., 600., 650.};  
  Float_t mLSP0[NSignals] =  {   1., 100., 175.,   1.,  75., 150.,  50.,   1., 100., 150.,  50.,   1., 200., 250.,   1.};

  if(iSample >= NBkgs) {
    sampleName[NBkgs] += "-mStop"; sampleName[NBkgs] += mStop0[iSample - NBkgs];
    sampleName[NBkgs] += "mLSP"; sampleName[NBkgs] += mLSP0[iSample - NBkgs];
  }

  const int NmtCut = 4;
  const int NnJetCut = 3;
  const int NtopnessCut = 1;
  const int Nmt2wCut = 4;
  const int NyCut = 4;
  const int NdphiCut = 3;
  const int NdrlblCut = 3;
  const int NdrlbgCut = 1;
  const int Nchi2Cut = 1;
  const int NmetCut = 3;
  const int Nm3Cut = 3;
  const int NcentralityCut = 3;
  const int NmlbCut = 3;
  const int NSR = NmtCut * NnJetCut * NtopnessCut * Nmt2wCut * NyCut * NdphiCut * NdrlblCut * NdrlbgCut * Nchi2Cut * NmetCut * Nm3Cut * NcentralityCut * NmlbCut;

  Double_t mtCut[NmtCut]           = {100., 120., 160., 220.};
  Double_t nJetCut[NnJetCut]       = { 3., 4., 5.};
  Double_t topnessCut[NtopnessCut] = {-15.};
  Double_t mt2wCut[Nmt2wCut]       = {0., 200., 250., 300.};
  Double_t yCut[NyCut]             = {0., 8., 10., 13.};
  Double_t dphiCut[NdphiCut]       = {0., 0.8, 1.};
  Double_t drlblCut[NdrlblCut]     = {5., 2., 1.5};
  //Double_t drlbgCut[NdrlbgCut]     = {0., 1.5, 2.};
  Double_t drlbgCut[NdrlbgCut]     = {0.};
  Double_t chi2Cut[Nchi2Cut]       = {999999999.};
  Double_t metCut[NmetCut]         = {150., 250., 350.};
  Double_t m3Cut[Nm3Cut]           = {0., 350., 500.};
  Double_t centralityCut[NcentralityCut]  = {-1., 0.5, 0.7};
  Double_t mlbCut[NmlbCut]         = { 999999999, 120., 100.};

  /*Float_t tt[NmtCut][NnJetCut][NtopnessCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][NdrlbgCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t tb[NmtCut][NnJetCut][NtopnessCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][NdrlbgCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t bb[NmtCut][NnJetCut][NtopnessCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][NdrlbgCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t tt2[NmtCut][NnJetCut][NtopnessCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][NdrlbgCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut]= {};
  Float_t tb2[NmtCut][NnJetCut][NtopnessCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][NdrlbgCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t bb2[NmtCut][NnJetCut][NtopnessCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][NdrlbgCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};*/

  Float_t tt[NmtCut][NnJetCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t tb[NmtCut][NnJetCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t bb[NmtCut][NnJetCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t tt2[NmtCut][NnJetCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut]= {};
  Float_t tb2[NmtCut][NnJetCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};
  Float_t bb2[NmtCut][NnJetCut][Nmt2wCut][NyCut][NdphiCut][NdrlblCut][Nchi2Cut][NmetCut][Nm3Cut][NcentralityCut][NmlbCut] = {};

  for ( int imt = 0; imt < NmtCut ; imt++){
    for ( int ijet = 0; ijet < NnJetCut ; ijet++){
      for ( int itopness = 0; itopness < NtopnessCut ; itopness++){
	for ( int imt2w = 0; imt2w < Nmt2wCut ; imt2w++){
	  for ( int iy = 0; iy < NyCut ; iy++){
	    for ( int idphi = 0; idphi < NdphiCut ; idphi++){
	      for ( int idrlbl = 0; idrlbl < NdrlblCut ; idrlbl++){
		for ( int idrlbg = 0; idrlbg < NdrlbgCut ; idrlbg++){
		  for ( int ichi2 = 0; ichi2 < Nchi2Cut ; ichi2++){
		    for ( int imet = 0; imet < NmetCut ; imet++){
		      for ( int im3 = 0; im3 < Nm3Cut ; im3++){
			for ( int icentrality = 0; icentrality < NcentralityCut ; icentrality++){
			  for ( int imlb = 0; imlb < NmlbCut ; imlb++){
			    /*tt[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] = 0.;
			    tb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] = 0.;
			    bb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] = 0.;
			    tt2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] = 0.;
			    tb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] = 0.;
			    bb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] = 0.;*/

			    tt[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] = 0.;
			    tb[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] = 0.;
			    bb[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] = 0.;
			    tt2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] = 0.;
			    tb2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] = 0.;
			    bb2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] = 0.;

			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  Float_t lumi=19500.;
  Float_t weight = 1.;
  /////////////////////////////////////////////////////
  //  Variable Definition
  ///////////////////////////////////////////////////// 

  Float_t globalWeight = 0.;
  Float_t FE = 0.;
  Float_t triggerWeight = 0.;
  Float_t puWeight = 0.;
  Float_t isrWeight = 0.;
  Float_t topPtWeight = 0.;

  Float_t lepFromTop = 0.;
  Float_t charginos = 0.;

  Float_t njets = 0.;
  Float_t jet1 = 0.;
  Float_t jet2 = 0.;
  Float_t jet3 = 0.;
  Float_t jet4 = 0.;

  Float_t nbjets = 0.;
  Float_t bjet1 = 0.;
  Float_t bjetHighestDisc = 0.;
  
  Float_t lPt = 0.;
  Float_t lEta = 0.;
  Float_t lRelIso = 0.;

  Float_t isoTrack = 0.;
  Float_t tauVeto = 0.;
  
  Float_t rawmet = 0.;
  Float_t typeImet = 0.;
  Float_t phiCorrMet = 0.;
  
  Float_t ht = 0.;
  Float_t htRatio = 0.;
  Float_t meff = 0.;
  Float_t y = 0.;
  
  Float_t mt = 0.;
  Float_t mlb1 = 0.;
  Float_t mlb = 0.;
  Float_t m3b = 0.;
  Float_t m3 = 0.;
  Float_t centrality = 0.;
  Float_t mt2w = 0.;
  Float_t hadChi2 = 0.;
  Float_t topness = 0.;

  Float_t dphimin = 0.;
  Float_t drlb1 = 0.;
  Float_t drlbmin = 0.;
  
  Float_t mStop = 0.;
  Float_t mLSP = 0.;

  Int_t pdgIdLep1 = 0;
  Char_t kinRegionFlag = false;
  Char_t searchRegionFlag = false;

  /////////////////////////////////////////////////////
  //  Input Definition
  /////////////////////////////////////////////////////
  TString inFileName = mainDir; 
  if(iSample < NBkgs) inFileName += sampleName[iSample];
  else inFileName += sampleName[NBkgs];
  inFileName +=".root";
  cout<<inFileName<<endl;

  TFile* inFile = new TFile(inFileName,"READ");
  if (!inFile->IsOpen()){
    std::cout<<"not open"<<std::endl;
  }

  TTree* tree;
  tree= (TTree*)inFile->Get("NoSystematic/bonsai");
  //===========================================
  if(pcp)cout<<"inputs set!"<<endl;
  int N = tree->GetEntries();  cout<<"THERE ARE "<<N<<" EVENTS IN "<<inFileName<<endl;
    
  tree->SetBranchAddress( "GlobalWeight", &globalWeight);
  tree->SetBranchAddress( "FE", &FE);
  tree->SetBranchAddress( "TriggerWeight", &triggerWeight);
  tree->SetBranchAddress( "PUWeight", &puWeight);
  tree->SetBranchAddress( "isrWeight", &isrWeight);
  tree->SetBranchAddress( "topPtWeight", &topPtWeight);

  tree->SetBranchAddress( "LepFromTop", &lepFromTop);
  tree->SetBranchAddress( "Charginos", &charginos);

  tree->SetBranchAddress( "njets", &njets);
  tree->SetBranchAddress( "jet1", &jet1);
  tree->SetBranchAddress( "jet2", &jet2);
  tree->SetBranchAddress( "jet3", &jet3);
  tree->SetBranchAddress( "jet4", &jet4);

  tree->SetBranchAddress( "nbjets", &nbjets);
  tree->SetBranchAddress( "bjet1", &bjet1);
  tree->SetBranchAddress( "bjetHighestDisc", &bjetHighestDisc);

  tree->SetBranchAddress( "lPt", &lPt);
  tree->SetBranchAddress( "lEta", &lEta);
  tree->SetBranchAddress( "lRelIso", &lRelIso);

  tree->SetBranchAddress( "phiCorrMet", &phiCorrMet);
  
  tree->SetBranchAddress( "ht", &ht);
  tree->SetBranchAddress( "htRatio", &htRatio);
  tree->SetBranchAddress( "meff", &meff);
  tree->SetBranchAddress( "y", &y);

  tree->SetBranchAddress( "mt", &mt);
  tree->SetBranchAddress( "mlb1", &mlb1);
  tree->SetBranchAddress( "mlb", &mlb);
  tree->SetBranchAddress( "m3b", &m3b);
  tree->SetBranchAddress( "m3", &m3);
  tree->SetBranchAddress( "centrality", &centrality);
  tree->SetBranchAddress( "mt2w", &mt2w);
  tree->SetBranchAddress( "hadChi2", &hadChi2);
  tree->SetBranchAddress( "topness", &topness);

  tree->SetBranchAddress( "dphimin", &dphimin);
  tree->SetBranchAddress( "drlb1", &drlb1);
  tree->SetBranchAddress( "drlbmin", &drlbmin);

  tree->SetBranchAddress("mStop",&mStop);
  tree->SetBranchAddress("mLSP",&mLSP);

  tree->SetBranchAddress("pdgIdLep1",&pdgIdLep1);
  tree->SetBranchAddress("kinRegion",&kinRegionFlag);
  tree->SetBranchAddress("searchRegionPost",&searchRegionFlag);


  /////////////////////////////////////////////////////
  //  Output Definition
  ///////////////////////////////////////////////////// 

  TString outFileName = "./";  
  if(iSample < NBkgs) outFileName += sampleName[iSample];
  else outFileName += sampleName[NBkgs];
  outFileName +="_optimization-1Step.root";
  cout<<outFileName<<endl;

  TFile* outFile = new TFile( outFileName, "RECREATE"); 
  TTree* outTree = new TTree("Optimization", "Optimization subTree");

  /////////////////////////////////////////////////////
  //  Branch Definition
  ///////////////////////////////////////////////////// 

  Float_t mtCut0 = 0.;
  Float_t nJetCut0 = 0.;
  Float_t topnessCut0 = 0.;
  Float_t mt2wCut0 = 0.;
  Float_t yCut0 = 0.;
  Float_t dphiCut0 = 0.;
  Float_t drlblCut0 = 0.;
  Float_t drlbgCut0 = 0.;
  Float_t chi2Cut0 = 0.;
  Float_t metCut0 = 0.;
  Float_t m3Cut0 = 0.;
  Float_t centralityCut0 = 0.;
  Float_t mlbCut0 = 0.;

  TString SR0("");

  Float_t tt0 = 0.;
  Float_t tb0 = 0.;
  Float_t bb0 = 0.;
  Float_t tt20 = 0.;
  Float_t tb20 = 0.;
  Float_t bb20 = 0.;

  outTree->Branch("mtCut", &mtCut0, "mtCut/F");
  outTree->Branch("nJetCut", &nJetCut0, "nJetCut/F");
  outTree->Branch("topnessCut", &topnessCut0, "topnessCut/F");
  outTree->Branch("mt2wCut", &mt2wCut0, "mt2wCut/F");
  outTree->Branch("yCut", &yCut0, "yCut/F");
  outTree->Branch("dphiCut", &dphiCut0, "dphiCut/F");
  outTree->Branch("drlblCut", &drlblCut0, "drlblCut/F");
  outTree->Branch("drlbgCut", &drlbgCut0, "drlbgCut/F");
  outTree->Branch("chi2Cut", &chi2Cut0, "chi2Cut/F");
  outTree->Branch("metCut", &metCut0, "metCut/F");
  outTree->Branch("m3Cut", &m3Cut0, "m3Cut/F");
  outTree->Branch("centralityCut", &centralityCut0, "centralityCut/F");
  outTree->Branch("mlbCut", &mlbCut0, "mlbCut/F"); 
  outTree->Branch("sr", "TString", &SR0);
  
  outTree->Branch("tt", &tt0, "tt/F");
  outTree->Branch("tb", &tb0, "tb/F");
  outTree->Branch("bb", &bb0, "bb/F");
  outTree->Branch("tt2", &tt20, "tt2/F");
  outTree->Branch("tb2", &tb20, "tb2/F");
  outTree->Branch("bb2", &bb20, "bb2/F");

  /////////////////////////////////////////////////////
  //  Event Loop
  ///////////////////////////////////////////////////// 
  for (int ievt=0;ievt<N;++ievt){
    if (ievt%13454 == 0) {
      cout<<"Event number "<<ievt<<"\r"<<flush;
    }
    
    tree->GetEntry(ievt);

    weight = globalWeight * triggerWeight * puWeight * lumi;
    if (iSample < 2) weight *= topPtWeight;
    if (iSample > 3 ) {
      if ( fabs(mStop -  mStop0[iSample - NBkgs]) > 0.001) continue;
      if ( fabs(mLSP -  mLSP0[iSample - NBkgs]) > 0.001) continue;
      if (charginos == 0 && FE > 0.999) continue;
      weight *= isrWeight;
    }
    
    if (lRelIso > 0.1) continue;
    if (!searchRegionFlag) continue;

    for ( int imt = 0; imt < NmtCut ; imt++){
      if (mt < mtCut[imt]) continue;

      for ( int ijet = 0; ijet < NnJetCut ; ijet++){
	if (njets < nJetCut[ijet]) continue;

	for ( int itopness = 0; itopness < NtopnessCut ; itopness++){
	  if (topness < topnessCut[itopness]) continue;
	  
	  for ( int imt2w = 0; imt2w < Nmt2wCut ; imt2w++){
	    if (mt2w < mt2wCut[imt2w]) continue;

	    for ( int iy = 0; iy < NyCut ; iy++){
	      if (y < yCut[iy]) continue;
	      
	      for ( int idphi = 0; idphi < NdphiCut ; idphi++){
		if (dphimin < dphiCut[idphi]) continue;

		for ( int idrlbl = 0; idrlbl < NdrlblCut ; idrlbl++){
		  if (drlb1 > drlblCut[idrlbl]) continue;

		  for ( int idrlbg = 0; idrlbg < NdrlbgCut ; idrlbg++){
		    if (drlb1 < drlbgCut[idrlbg]) continue;
		    
		    for ( int ichi2 = 0; ichi2 < Nchi2Cut ; ichi2++){
		      if (hadChi2 > chi2Cut[ichi2]) continue;
		      
		      for ( int imet = 0; imet < NmetCut ; imet++){	     
			if (phiCorrMet < metCut[imet]) continue;
			
			for ( int im3 = 0; im3 < Nm3Cut ; im3++){
			  if (m3 < m3Cut[im3]) continue;

			  for ( int icentrality = 0; icentrality < NcentralityCut ; icentrality++){
			    if (centrality < centralityCut[icentrality]) continue;

			    for ( int imlb = 0; imlb < NmlbCut ; imlb++){
			      if (mlb > mlbCut[imlb]) continue;

			      /* if (charginos == 0){
				tt[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] += weight;
				tt2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] += weight * weight;
			      }
			      if (charginos == 1){
				tb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] += weight;
				tb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] += weight * weight;
			      }
			      if (charginos == 2){
				bb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] += weight;
				bb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb] += weight * weight;
			      }*/

			      if (charginos == 0){
				tt[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] += weight;
				tt2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] += weight * weight;
			      }
			      if (charginos == 1){
				tb[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] += weight;
				tb2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] += weight * weight;
			      }
			      if (charginos == 2){
				bb[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] += weight;
				bb2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb] += weight * weight;
			      }

			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  } 
	}
      }
    }
  }

  for ( int imt = 0; imt < NmtCut ; imt++){
    for ( int ijet = 0; ijet < NnJetCut ; ijet++){
      for ( int itopness = 0; itopness < NtopnessCut ; itopness++){
	for ( int imt2w = 0; imt2w < Nmt2wCut ; imt2w++){
	  for ( int iy = 0; iy < NyCut ; iy++){
	    for ( int idphi = 0; idphi < NdphiCut ; idphi++){
	      for ( int idrlbl = 0; idrlbl < NdrlblCut ; idrlbl++){
		for ( int idrlbg = 0; idrlbg < NdrlbgCut ; idrlbg++){
		  for ( int ichi2 = 0; ichi2 < Nchi2Cut ; ichi2++){
		    for ( int imet = 0; imet < NmetCut ; imet++){	     
		      for ( int im3 = 0; im3 < Nm3Cut ; im3++){			
			for ( int icentrality = 0; icentrality < NcentralityCut ; icentrality++){
			  for ( int imlb = 0; imlb < NmlbCut ; imlb++){

			    mtCut0 = mtCut[imt];
			    nJetCut0 = nJetCut[ijet];    
			    topnessCut0 = topnessCut[itopness];
			    mt2wCut0 = mt2wCut[imt2w];
			    yCut0 = yCut[iy];
			    dphiCut0 =  dphiCut[idphi];
			    drlblCut0 =  drlblCut[idrlbl];
			    drlbgCut0 =  drlbgCut[idrlbg];
			    chi2Cut0 =  chi2Cut[ichi2];
			    metCut0 = metCut[imet];
			    m3Cut0 = m3Cut[im3];   
			    centralityCut0 = centralityCut[icentrality];
			    mlbCut0 = mlbCut[imlb];
			    
			    TString tmp = "MT-"; tmp += mtCut[imt]; tmp += "_";
			    tmp += "Njet-"; tmp += nJetCut[ijet]; tmp += "_";
			    tmp += "Topness-"; tmp += topnessCut[itopness]; tmp += "_";
			    tmp += "Mt2w-"; tmp += mt2wCut[imt2w]; tmp += "_";
			    tmp += "Y-"; tmp += yCut[iy]; tmp += "_";
			    tmp += "Dphi-"; tmp += dphiCut[idphi]; tmp += "_";
			    tmp += "Drlbl-"; tmp += drlblCut[idrlbl]; tmp += "_";
			    tmp += "Drlbg-"; tmp += drlbgCut[idrlbg]; tmp += "_";
			    tmp += "Chi2-"; tmp += chi2Cut[ichi2]; tmp += "_";
			    tmp += "Met-"; tmp += metCut[imet];
			    tmp += "M3-"; tmp += m3Cut[im3];
			    tmp += "Centrality-"; tmp += centralityCut[icentrality];
			    tmp += "Mlb-"; tmp += mlbCut[imlb];
		      
			    SR0 = tmp;
			    
			    /*tt0 = tt[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb];
			    tb0 = tb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb];
			    bb0 = bb[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb];
			    tt20 = tt2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb];
			    tb20 = tb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb];
			    bb20 = bb2[imt][ijet][itopness][imt2w][iy][idphi][idrlbl][idrlbg][ichi2][imet][im3][icentrality][imlb];*/
			    
			    tt0 = tt[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb];
			    tb0 = tb[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb];
			    bb0 = bb[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb];
			    tt20 = tt2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb];
			    tb20 = tb2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb];
			    bb20 = bb2[imt][ijet][imt2w][iy][idphi][idrlbl][ichi2][imet][im3][icentrality][imlb];

			    outTree->Fill();
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  outFile->cd();
  outTree->Write();

  outFile->Close();    
  inFile->Close();
  
  return 0;
}

