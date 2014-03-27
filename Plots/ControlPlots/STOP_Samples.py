#! /usr/bin/env pyth
import sys, os
import commands
import subprocess
from time import asctime
import random
import copy
#from SelectionSamples import SelectionSample
#from SampleClasses import SampleTask,Sample,FakeRatioSample,SelectionSample
from SampleClasses import Sample
from BreakDownInputPath import *
from GetDate import *

#==========================================
##  Data
#==========================================
def Data_STOP_NoTail(options):
    
    Samp=Sample()
    Samp.isData=True
    #
    Samp.AddRootFile("Data/NoSub/Data_NoSub_STOP_NoTail.root")

    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def SinleMu_STOP_NoTail(options):
    
    Samp=Sample()
    Samp.isData=True
    #
    Samp.AddRootFile("SingleMu/NoSub/SingleMu_NoSub_STOP_NoTail.root")

    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def SingleElectron_STOP_NoTail(options):
    
    Samp=Sample()
    Samp.isData=True
    #
    Samp.AddRootFile("SingleElectron/NoSub/SingleElectron_NoSub_STOP_NoTail.root")

    Samp.ReweightAndMerge(options)
    return Samp
#==========================================


#==========================================
##  MC Background
#==========================================
def DiLep_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    
    Samp.AddRootFile("DiLep/NoSub/DiLep_NoSub_STOP_NoTail.root")

    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def OneLep_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    
    Samp.AddRootFile("OneLep/NoSub/OneLep_NoSub_STOP_NoTail.root")

    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def WJets_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("WJets/NoSub/WJets_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#=========================================
def Rare_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("Rare/NoSub/Rare_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def QCD_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("QCD/NoSub/QCD_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def DrellYan_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("DrellYan/NoSub/DrellYan_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================

#==========================================
##  Signal Samples
#==========================================
def T2tb_mStop175mLSP50br50_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop175mLSP50br50/NoSub/T2tb-mStop175mLSP50br50_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop175mLSP50br100_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop175mLSP50br100/NoSub/T2tb-mStop175mLSP50br100_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop200mLSP25br50_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop200mLSP25br50/NoSub/T2tb-mStop200mLSP25br50_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop200mLSP25br100_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop200mLSP25br100/NoSub/T2tb-mStop200mLSP25br100_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop325mLSP100br50_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop325mLSP100br50/NoSub/T2tb-mStop325mLSP100br50_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop325mLSP100br100_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop325mLSP100br100/NoSub/T2tb-mStop325mLSP100br100_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop550mLSP1br50_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop550mLSP1br50/NoSub/T2tb-mStop550mLSP1br50_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================
def T2tb_mStop550mLSP1br100_STOP_NoTail(options):

    Samp=Sample()
    Samp.isData=False
    #
    Samp.AddRootFile("T2tb-mStop550mLSP1br100/NoSub/T2tb-mStop550mLSP1br100_NoSub_STOP_NoTail.root")
    #
    Samp.ReweightAndMerge(options)
    return Samp
#==========================================

