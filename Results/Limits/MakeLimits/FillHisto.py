#!/usr/bin/env python

import math
import os, sys
from subprocess import call
from ROOT import ROOT, gDirectory, gROOT, TFile, TDirectory, TH1F, TH2F, TIter, TKey

gROOT.SetBatch()


def GetKeys( self, dir = "" ):
        self.cd(dir)
        return [key for key in gDirectory.GetListOfKeys()]
TDirectory.GetKeys = GetKeys


def makeConfig( mStop, mLSP, dir):
    datah = TH1F()
    bkgh = TH1F()
    sigh = TH2F()
    unch = TH2F()

    dir.GetObject( "data", datah)
    dir.GetObject( "bkg", bkgh)
    dir.GetObject( "sig", sigh)
    dir.GetObject( "unc", unch)
    
    data = datah.GetBinContent(1)

    bkg = bkgh.GetBinContent(1)
    bkgUnc = bkgh.GetBinError(1)    

    sig = sigh.GetBinContent(sigh.FindBin(float(mStop),float(mLSP)))
    sigUnc = unch.GetBinContent(sigh.FindBin(float(mStop),float(mLSP)))

    if (sig < 1 ):
    	    return 1    
    
    sigC= 1.
    bkgC= 1.
    
    if (sig > 0):  sigC= 1 + sigUnc / sig
    if (bkg > 0):  bkgC= 1 + bkgUnc / bkg

    f = open('template.card', 'r')
    h = open('my.card', 'w')

    for line in f:
	    line=line.replace( "$OBSERVATION", str(data))
	    
	    line=line.replace( "$SIG_", str(sig))
	    line=line.replace( "$BKG_", str(bkg))
	    
	    line=line.replace( "$SIGC_", str(sigC))
	    line=line.replace( "$BKGC_", str(bkgC))

	    h.write(line)

    h.close()
    f.close()

    return 0

mStop = sys.argv[1]
mLSP = sys.argv[2]

inFile = TFile("InputForLimits.root")
outFile = TFile("Exclusion.root", "recreate")

brList = inFile.GetKeys()
for br in brList:
	cl = gROOT.GetClass(br.GetClassName());
	if (cl.InheritsFrom("TDirectory")):
		inDir = br.ReadObj()

		outFile.mkdir(br.GetName())
		brDir=outFile.GetDirectory(br.GetName())
		
		srList = inDir.GetKeys()
		for sr in srList:
			cl = gROOT.GetClass(sr.GetClassName());
			if (cl.InheritsFrom("TDirectory")):

				print br.GetName() + "/" + sr.GetName()
				
				inDir = sr.ReadObj()				
				brDir.mkdir(sr.GetName())
				srDir = brDir.GetDirectory(sr.GetName())
				srDir.cd()

				sigh = TH2F()
				inDir.GetObject( "sig", sigh)				

				rm2sh = TH2F(sigh)
				rm2sh.SetName("rm2s")
				rm2sh.SetTitle("rm2s")
				rm2sh.Reset()

				rm1sh = TH2F(rm2sh)
				rm1sh.SetName("rm1s")
				rm1sh.SetTitle("rm1s")
			
				rmeanh = TH2F(rm2sh)
				rmeanh.SetName("rmean")
				rmeanh.SetTitle("rmean")
			
				rp1sh = TH2F(rm2sh)
				rp1sh.SetName("rp1s")
				rp1sh.SetTitle("rp1s")
				
				rp2sh = TH2F(rm2sh)
				rp2sh.SetName("rp2s")
				rp2sh.SetTitle("rp2s")				

				obsh = TH2F(rm2sh)
				obsh.SetName("obs")
				obsh.SetTitle("obs")	
				
				if (not makeConfig( mStop, mLSP, inDir)):

					os.system("./expected.sh")
				
					resultsf = TFile("results.root", "read")
				
					resultsh=TH1F()
					resultsf.GetObject("results",resultsh)

					rm2sh.SetBinContent(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinContent(1))
					rm1sh.SetBinContent(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinContent(2))
					rmeanh.SetBinContent(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinContent(3))
					rp1sh.SetBinContent(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinContent(4))
					rp2sh.SetBinContent(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinContent(5))

					resultsf.Close()

					os.system("./observed.sh")
				
					resultsf = TFile("results.root", "read")
				
					resultsh=TH1F()
					resultsf.GetObject("results",resultsh)

					obsh.SetBinContent(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinContent(6))
					obsh.SetBinError(sigh.FindBin(float(mStop),float(mLSP)), resultsh.GetBinError(6))

					resultsf.Close()					
					
				srDir.cd()
				rm2sh.Write()
				rm1sh.Write()
				rmeanh.Write()
				rp1sh.Write()
				rp2sh.Write()
				obsh.Write()
