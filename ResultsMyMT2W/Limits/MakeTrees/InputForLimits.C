#include "TH2F.h"
#include "TTree.h"
#include "TFile.h"
#include "TMath.h"

#include <iostream>
#include <iomanip>
#include <vector>

bool pcp = true;
using namespace std;

class Systematics{
public:
  Systematics(){};

  int Calc(){
    Float_t a = 4. * BR * BR;
    Float_t b = 4. * BR * (1 - BR);
    Float_t c = 4. * (1 - BR) * (1 - BR);

    sigh = new TH2F( *h[0]);
    sigh->SetName("sig");
    sigh->SetTitle("sig");    
    sigh->Scale(a);
    sigh->Add( h[1], b);
    sigh->Add( h[2], c);
    
    effh = new TH2F( *sigh);
    effh->SetName("eff");
    effh->SetTitle("eff");  
    effh->Divide(sig_toth);

    sigLh = new TH2F( *Lh);
    sigLh->SetName("sigL");
    sigLh->SetTitle("sigL");    
    sigLh->Scale(a);

    effLh = new TH2F( *sigLh);
    effLh->SetName("effL");
    effLh->SetTitle("effL");  
    effLh->Divide(sig_toth);

    sigRh = new TH2F( *Rh);
    sigRh->SetName("sigR");
    sigRh->SetTitle("sigR");    
    sigRh->Scale(a);

    effRh = new TH2F( *sigRh);
    effRh->SetName("effR");
    effRh->SetTitle("effR");  
    effRh->Divide(sig_toth);

    TString histoname;
    for ( int ishift = 0; ishift < 2; ishift++){
      shiftSigh[ishift] = new TH2F( *shifth[ishift][0]);
      histoname = "shiftSig"; histoname += ishift;
      shiftSigh[ishift]->SetName( histoname);
      shiftSigh[ishift]->SetTitle( histoname);
      shiftSigh[ishift]->Scale(a);
      shiftSigh[ishift]->Add( shifth[ishift][1], b);
      shiftSigh[ishift]->Add( shifth[ishift][2], c);
    }

    sysh = new TH2F( *sigh); 
    sysh->Reset();
    sysh->SetName("sys");
    sysh->SetTitle("sys");

    Float_t cen = 0.;    
    Float_t sys = 0.;

    Float_t up = 0.;
    Float_t down = 0.;

    Float_t min = 0.;
    Float_t max = 0.;

    for (int ibin = 0; ibin < sigh->GetSize(); ibin++){
      cen = sigh->GetBinContent(ibin);
      up = shiftSigh[0]->GetBinContent(ibin);
      down = shiftSigh[1]->GetBinContent(ibin);

      if ( up < down) {
	min = up;
	max = down;
      }
      else {
	min = down;
	max = up;
      }

      sys = sqrt( (cen - (max + min) / 2.) * (cen - (max + min) / 2.) + (max - min) * (max - min) / 4.);
      
      sysh->SetBinContent( ibin, sys);
    }
  
    return 0;
  }
  
  Float_t BR;
  
  TH2F* h[3];
  TH2F* Lh;
  TH2F* Rh;

  TH2F* sig_toth;
  TH2F* sigh;
  TH2F* sigLh;
  TH2F* sigRh;
  TH2F* effh;
  TH2F* effLh;
  TH2F* effRh;

  TH2F* shifth[2][3];
  TH2F* shiftSigh[2];

  TH2F* sysh;
};
  

int InputForLimits(){
  TH1::SetDefaultSumw2(true);
  if(pcp)cout<<"going to set inputs"<<endl;

  Int_t NBR = 4;
  Float_t BR[] = { 1., 0.75, 0.5, 0.25};

  TFile* bkgFile = new TFile( "../../BkgPrediction/BkgPrediction.root", "READ"); 
  TTree* bkgTree;
  bkgFile->GetObject( "ElAndMu", bkgTree);
  
  Float_t bkg = 0.;
  Float_t bkgTotUnc = 0.;

  Float_t obs = 0.;

  bkgTree->SetBranchAddress( "srData", &obs);
  bkgTree->SetBranchAddress( "srAllBkgCorr", &bkg);
  bkgTree->SetBranchAddress( "TotUnc", &bkgTotUnc);

  TFile* sigFile = new TFile( "../../SignalSystematics/SignalSys.root", "READ"); 

  std::vector<std::vector<TString> > sysColl;
  std::vector<TString> sys;
  sys.push_back(TString("JES_Up"));
  sys.push_back(TString("JES_Down"));
  sysColl.push_back(sys);

  sys.clear();
  sys.push_back(TString("BTagReweight_UpBC"));
  sys.push_back(TString("BTagReweight_DownBC"));
  sysColl.push_back(sys);

  sys.clear();
  sys.push_back(TString("BTagReweight_UpLight"));
  sys.push_back(TString("BTagReweight_DownLight"));
  sysColl.push_back(sys);

  std::vector<TString> decayMode;
  decayMode.push_back(TString("tt"));
  decayMode.push_back(TString("tb"));
  decayMode.push_back(TString("bb"));
		      
  Systematics* systematics[3];
  TString dirname;
  TString histoname;

  TDirectory* srDir;
  TDirectory* histoDir;

  TFile* outFile = new TFile( "InputForLimits.root", "RECREATE");
  TDirectory* outBRDir;
  TDirectory* outSRDir;

  TH1F* datah = new TH1F( "data", "data", 1, 0., 1.);
  TH1F* bkgh = new TH1F( "bkg", "bkg", 1, 0., 1.);

  TH2F* sigh;
  TH2F* sig_toth;
  TH2F* effh;
  
  TH2F* sigLh;
  TH2F* effLh;

  TH2F* sigRh;
  TH2F* effRh;
  
  TH2F* jesh;
  TH2F* btagBCh;
  TH2F* btagLighth;
  TH2F* btagh;
  
  TH2F* sysh;

  TH2F* unch;

  TH2F* jesPercenth;
  TH2F* btagBCPercenth;
  TH2F* btagLightPercenth;
  TH2F* btagPercenth;
  
  TH2F* sysPercenth;

  TH2F* uncPercenth;

  Float_t sig = 0.;
  Float_t stat = 0.;

  Float_t jes = 0.;
  Float_t bc = 0.;
  Float_t light = 0.;
  
  Float_t unc = 0.;

  int N = bkgTree->GetEntries();

  for ( int ibr = 0; ibr < NBR; ibr++){
    dirname = ""; dirname += BR[ibr];
    outFile->mkdir(dirname);
    outBRDir = outFile->GetDirectory(dirname);

    for ( int iSR = 0; iSR < 8; iSR++){
      bkgTree->GetEntry(iSR);      
      datah->SetBinContent( 1, obs);
      bkgh->SetBinContent( 1, bkg);
      bkgh->SetBinError( 1, bkgTotUnc);

      dirname = ""; dirname += iSR;
      outBRDir->mkdir(dirname);
      outSRDir = outBRDir->GetDirectory(dirname);
      
      dirname = ""; dirname += iSR; dirname += ".root";      
      srDir = sigFile->GetDirectory( dirname);
      for ( int isys = 0; isys < (int) sysColl.size(); isys++){
	systematics[isys] = new Systematics();
	systematics[isys]->BR = BR[ibr];

	histoDir = srDir->GetDirectory( "NoSystematic");
	for ( int idecay = 0; idecay < (int) decayMode.size(); idecay++){
	  histoDir->GetObject( decayMode.at(idecay), systematics[isys]->h[idecay]);
	  histoDir->GetObject( decayMode.at(idecay) + "l", systematics[isys]->Lh);
	  histoDir->GetObject( decayMode.at(idecay) + "r", systematics[isys]->Rh);
	}
	histoDir->GetObject( decayMode.at(0) + "l", systematics[isys]->Lh);
	histoDir->GetObject( decayMode.at(0) + "r", systematics[isys]->Rh);
	histoDir->GetObject( "sig_tot", systematics[isys]->sig_toth);

	for ( int ishift = 0; ishift < 2; ishift++){
	  histoDir = srDir->GetDirectory(sysColl.at(isys).at(ishift));
	  for ( int idecay = 0; idecay < (int) decayMode.size(); idecay++){
	    histoDir->GetObject( decayMode.at(idecay), systematics[isys]->shifth[ishift][idecay]);
	  }
	}
	systematics[isys]->Calc();
      }
      sig_toth = new TH2F( *systematics[0]->sig_toth);
      sig_toth->SetName("sig_tot");
      sig_toth->SetTitle("sig_tot");    

      sigh = new TH2F( *systematics[0]->sigh);
      sigh->SetName("sig");
      sigh->SetTitle("sig");

      effh = new TH2F( *systematics[0]->effh);
      effh->SetName("eff");
      effh->SetTitle("eff");

      sigLh = new TH2F( *systematics[0]->sigRh);
      sigLh->SetName("sigL");
      sigLh->SetTitle("sigL");

      effLh = new TH2F( *systematics[0]->effLh);
      effLh->SetName("effL");
      effLh->SetTitle("effL");

      sigRh = new TH2F( *systematics[0]->sigRh);
      sigRh->SetName("sigR");
      sigRh->SetTitle("sigR");

      effRh = new TH2F( *systematics[0]->effRh);
      effRh->SetName("effR");
      effRh->SetTitle("effR");

      jesh = new TH2F( *systematics[0]->sysh);
      jesh->SetName("jes");
      jesh->SetTitle("jes");

      btagBCh = new TH2F( *systematics[1]->sysh);
      btagBCh->SetName("btagBC");
      btagBCh->SetTitle("btagBC");
      
      btagLighth = new TH2F( *systematics[2]->sysh);
      btagLighth->SetName("btagLight");
      btagLighth->SetTitle("btagLight");
      
      btagh = new TH2F( *btagBCh);
      btagh->Reset();
      btagh->SetName("btag");
      btagh->SetTitle("btag");
 
      sysh = new TH2F( *jesh);
      sysh->Reset();
      sysh->SetName("sys");
      sysh->SetTitle("sys");

      unch = new TH2F( *jesh);
      unch->Reset();
      unch->SetName("unc");
      unch->SetTitle("unc");


      jesPercenth = new TH2F( *systematics[0]->sysh);
      jesPercenth->SetName("jesPercent");
      jesPercenth->SetTitle("jesPercent");
      jesPercenth->Divide( sigh);
      jesPercenth->Scale( 100.);      

      btagBCPercenth = new TH2F( *systematics[1]->sysh);
      btagBCPercenth->SetName("btagBCPercent");
      btagBCPercenth->SetTitle("btagBCPercent");
      btagBCPercenth->Divide( sigh);
      btagBCPercenth->Scale( 100.);   

      btagLightPercenth = new TH2F( *systematics[2]->sysh);
      btagLightPercenth->SetName("btagLightPercent");
      btagLightPercenth->SetTitle("btagLightPercent");
      btagLightPercenth->Divide( sigh);
      btagLightPercenth->Scale( 100.);   

      btagPercenth = new TH2F( *btagBCh);
      btagPercenth->Reset();
      btagPercenth->SetName("btagPercent");
      btagPercenth->SetTitle("btagPercent");
 
      sysPercenth = new TH2F( *jesh);
      sysPercenth->Reset();
      sysPercenth->SetName("sysPercent");
      sysPercenth->SetTitle("sysPercent");

      uncPercenth = new TH2F( *jesh);
      uncPercenth->Reset();
      uncPercenth->SetName("uncPercent");
      uncPercenth->SetTitle("uncPercent");

      for (int ibin = 0; ibin < sigh->GetSize(); ibin++){	
	sig = sigh->GetBinContent( ibin);
	stat = sigh->GetBinError( ibin);

	jes = jesh->GetBinContent(ibin);
	bc = btagBCh->GetBinContent(ibin);
	light = btagLighth->GetBinContent(ibin);
    
	unc = sqrt( bc * bc + light * light);
	btagh->SetBinContent( ibin, unc);
	btagPercenth->SetBinContent( ibin, unc / sig * 100.);
	
	unc = sqrt( jes * jes + bc * bc + light * light +
		    sig * sig * (0.044 * 0.044 + // Lumi
				 0.03  * 0.03  + // Trigger
				 0.05  * 0.05    // Lep Id
				 )
		    );       
	sysh->SetBinContent( ibin, unc);
	sysPercenth->SetBinContent( ibin, unc / sig * 100.);
	
	unc = sqrt( jes * jes + bc * bc + light * light + stat * stat +  
		    sig * sig * (0.044 * 0.044 + // Lumi
				 0.03  * 0.03  + // Trigger
				 0.05  * 0.05    // Lep Id
				 )
		    );	
	unch->SetBinContent( ibin, unc);
	uncPercenth->SetBinContent( ibin, unc / sig * 100.);
      }
      outSRDir->cd();

      datah->Write();
      bkgh->Write();

      sig_toth->Write();
      sigh->Write();
      effh->Write();
      
      sigLh->Write();
      effLh->Write();

      sigRh->Write();
      effRh->Write();

      jesh->Write();
      btagBCh->Write(); 
      btagLighth->Write();
      btagh->Write();
      
      sysh->Write();
      unch->Write();

      jesPercenth->Write();
      btagBCPercenth->Write(); 
      btagLightPercenth->Write();
      btagPercenth->Write();
      
      sysPercenth->Write();
      uncPercenth->Write();
    }
  }

  delete systematics[0];
  delete systematics[1];
  delete systematics[2];

  outFile->Close();
  sigFile->Close();
  bkgFile->Close();

  return 0;
}
