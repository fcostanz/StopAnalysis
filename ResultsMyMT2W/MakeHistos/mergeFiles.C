// Author: Rene Brun

#include "TROOT.h"
#include "TKey.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

  //Example of script showing how to copy all objects (including directories)
  //from a source file.
  //For each input file, a new directory is created in the current directory 
  //with the name of the source file.
  //After execution of:
  // root > .x copyFiles.C
  //the file result.root will contain 4 subdirectories:
  // "tot100.root", "hsimple.root", "hs1.root","hs2.root"
      
void CopyDir(TDirectory *source) {
   //copy all objects and subdirs of directory source as a subdir of the current directory   
   source->ls();
   TDirectory *savdir = gDirectory;
   TDirectory *adir = savdir->mkdir(source->GetName());
   adir->cd();
   //loop on all entries of this directory
   TKey *key;
   TIter nextkey(source->GetListOfKeys());
   while ((key = (TKey*)nextkey())) {
      const char *classname = key->GetClassName();
      TClass *cl = gROOT->GetClass(classname);
      if (!cl) continue;
      if (cl->InheritsFrom(TDirectory::Class())) {
         source->cd(key->GetName());
         TDirectory *subdir = gDirectory;
         adir->cd();
         CopyDir(subdir);
         adir->cd();
      } else if (cl->InheritsFrom(TTree::Class())) {
         TTree *T = (TTree*)source->Get(key->GetName());
         adir->cd();
         TTree *newT = T->CloneTree(-1,"fast");
         newT->Write();
      } else {
         source->cd();
         TObject *obj = key->ReadObj();
         adir->cd();
         obj->Write();
         delete obj;
     }
  }
  adir->SaveSelf(kTRUE);
  savdir->cd();
}
void CopyFile(const char *fname) {
   //Copy all objects and subdirs of file fname as a subdir of the current directory
   TDirectory *target = gDirectory;
   TFile *f = TFile::Open(fname);
   if (!f || f->IsZombie()) {
      printf("Cannot copy file: %s\n",fname);
      target->cd();
      return;
   }
   target->cd();
   CopyDir(f);
   delete f;
   target->cd();
}  
void mergeFiles(int iSample = 0) {
  
  const int NSamples = 12;  

  TString sample[NSamples];
  sample[0] = "SingleMu";
  sample[1] = "SingleElectron";
  sample[2] = "DiLep";
  sample[3] = "OneLep";
  sample[4] = "WJets";
  sample[5] = "Rare";
  sample[6] = "SemiLep";
  sample[7] = "SingleTop";
  sample[8] = "QCD";
  sample[8] = "DrellYan";
  sample[8] = "T2tb-mStop175mLSP50";
  sample[9] = "T2tb-mStop200mLSP25";
  sample[10] = "T2tb-mStop325mLSP100";
  sample[11] = "T2tb-mStop550mLSP1";

  //Int_t SR[] = {31, 82, 32, 39, 61, 46, 22, 35, 51, 14, 65, 75, 1, 72, 70, 12, 9, 27};
  //Int_t NSR = 18;

  TString fileName = sample[iSample]; fileName += ".root";

  TFile *f = new TFile( "./" + fileName,"recreate");

  /*for (int iSR = 0; iSR < NSR; iSR++){
    TString inFileName = ""; inFileName += SR[iSR]; inFileName += "/"; inFileName += fileName;
    CopyFile(inFileName);
    }*/

  for (int iSR = 0; iSR < 14; iSR++){
    TString inFileName = ""; inFileName += iSR; inFileName += "/"; inFileName += fileName;
    CopyFile(inFileName);
  }

  delete f;
}
