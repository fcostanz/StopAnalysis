#
import sys, os
import commands
from ROOT import gROOT, gDirectory,TCanvas, TF1,TH1F,TH1D,TH1I,TFile,TKey,TString,THStack,TList,TLegend,TPaveText,TIter
from FileMerger import MergeFiles
from BreakDownInputPath import *
from copy import *
#from Plotter import *
from Plot_DataBkgSigPlot import *
from Plot_SetRootObjectsOptions import GetColor
from Plot_PlottingStyles import *
from STOP_Samples import *
from CheckDirectory import *
from CreateHTML import *

if __name__=='__main__':

    gROOT.SetBatch()

    #
    #
    #
    datasample=[]
    #datasample.append("SingleMu")
    #datasample.append("SingleElectron")
    datasample.append("Combined")
    #
    #
    for lep in datasample:
        
         DataOptions={}
         DataOptions['RW']=False
         DataOptions['Merge']=True
         #import the data sample============================
         if lep=="SingleMu":
             Data=SingleMu_STOPmu_NoTail(DataOptions)
             Data.SampleLumi=19500
         elif lep=="SingleElectron":
             Data=SingleElectron_STOPel_NoTail(DataOptions)
             Data.SampleLumi=19500
         elif lep=="Combined":
             Data=Data_STOP_NoTail(DataOptions)
             Data.SampleLumi=19500

             
         #
         #
         #===================================================

         
         #==========DIRECTORY FOR THE PLOTS=================
         #if this dir does not exist, it will be created
         topPlotsDir=Data.SampleName+'/ControlPlots/'+Data.Estimation+'/'+Data.Tail

         #remove the control plots there:
         print "going to remove the controlplots dir"
         rm_command='rm -rf '+topPlotsDir
         print "the command will be", rm_command

         #answer=raw_input('proceed? Answer \'yes\' if yes:  ')
         answer='yes'
         if answer=='yes':
             out=commands.getoutput(rm_command)
         else:
              print 'the directory is not removed'
         #
         #
         htmlOptions={}
         #
         #
         #
         #
         #======IMPORT THE MC SAMPLES================
         MCsamples={}

         MCOptions={}
         MCOptions['RW']=False
         MCOptions['Merge']=True
         MCOptions['dataLumi']=[Data.SampleLumi]

         if lep=="mu" or lep=="singleMu":
             #TTJETS
             TTJets=TTJetsPOWHEG_STOPmu_NoTail(MCOptions)
             MCsamples['04TTJets']=[TTJets]
             #WJETS
             WJets=WJetsNJetsBinned_STOPmu_NoTail(MCOptions)
             MCsamples['03WJets']=[WJets]
             #SingleTop
             SingleTop=SingleTop_STOPmu_NoTail(MCOptions)
             MCsamples['02SingleTop']=[SingleTop]
             #Rare
             Rare=Rare_STOPmu_NoTail(MCOptions)
             MCsamples['01Rare']=[Rare]

             htmlOptions['sample']=Data.SampleName
             #htmlOptions['PtCut']=20.0
             #htmlOptions['METCut']=60.0

         elif lep=="el" or lep=="singleElectron":
             #TTJETS
             TTJets=TTJetsPOWHEG_STOPel_NoTail(MCOptions)
             MCsamples['04TTJets']=[TTJets]
             #WJETS
             WJets=WJetsNJetsBinned_STOPel_NoTail(MCOptions)
             MCsamples['03WJets']=[WJets]
             #SingleTop
             SingleTop=SingleTop_STOPel_NoTail(MCOptions)
             MCsamples['02SingleTop']=[SingleTop]
             #Rare
             Rare=Rare_STOPel_NoTail(MCOptions)
             MCsamples['01Rare']=[Rare]
         elif lep=="Combined":
             #DiLep
             DiLep=DiLep_STOP_NoTail(MCOptions)
             MCsamples['06DiLep']=[DiLep]
             #OneLep
             OneLep=OneLep_STOP_NoTail(MCOptions)
             MCsamples['05OneLep']=[OneLep]
             #WJETS
             WJets=WJets_STOP_NoTail(MCOptions)
             MCsamples['04WJets']=[WJets]
             #Rare
             Rare=Rare_STOP_NoTail(MCOptions)
             MCsamples['03Rare']=[Rare]
             #QCD
             QCD=QCD_STOP_NoTail(MCOptions)
             MCsamples['02QCD']=[QCD]
             #DrellYan
             DrellYan=DrellYan_STOP_NoTail(MCOptions)
             MCsamples['01DrellYan']=[DrellYan]

             htmlOptions['sample']=Data.SampleName
             #htmlOptions['PtCut']=30.0
             #htmlOptions['METCut']=100.0    
             
         #======IMPORT THE Signal SAMPLES=================
         #
         Sigsamples={}
         SigOptions={}
         SigOptions['RW']=False
         SigOptions['True']=True
         SigOptions['dataLumi']=[Data.SampleLumi]
         ## if lep=='mu' or lep=="singleMu":
##              #T2tt stop=400, LSP=50
##              t400n50=T2tt_t400n50_RA4bmu_NoTail(SigOptions)
##              Sigsamples['2T2tt_t400n50']=[t400n50]
##              #T2tt stop=500, LSP=100
##              t500n100=T2tt_t500n100_RA4bmu_NoTail(SigOptions)
##              Sigsamples['1T2tt_t500n100']=[t500n100]
##          elif lep=='el' or lep=="singleEl":
##              #T2tt stop=400, LSP=50
##              t400n50=T2tt_t400n50_RA4bel_NoTail(SigOptions)
##              Sigsamples['2T2tt_t400n50']=[t400n50]
##              #T2tt stop=500, LSP=100
##              t500n100=T2tt_t500n100_RA4bel_NoTail(SigOptions)
##              Sigsamples['1T2tt_t500n100']=[t500n100]
         if lep=="Combined":
             #T2tb stop=550, LSP=1, br=50
             t550n1b50=T2tb_mStop550mLSP1br50_STOP_NoTail(SigOptions)
             Sigsamples['2T2tb-mStop550mLSP1br50']=[t550n1b50]
             #T2tb stop=325, LSP=100, br=50
             t325n100b50=T2tb_mStop325mLSP100br50_STOP_NoTail(SigOptions)
             Sigsamples['1T2tb-mStop325mLSP100br50']=[t325n100b50]
              
         #===================================================
         #
         #
         #
         #OPEN THE ROOT FILES===============================
         datafile=TFile(Data.RootFile,"READ")
         #
         #
         #MC bkg
         MCsamples['06DiLep'].append(TFile(DiLep.RootFile,"READ"))
         MCsamples["05OneLep"].append(TFile(OneLep.RootFile,"READ"))
         MCsamples['04WJets'].append(TFile(WJets.RootFile,"READ"))
         MCsamples['03Rare'].append(TFile(Rare.RootFile,"READ")) 
         MCsamples['02QCD'].append(TFile(QCD.RootFile,"READ"))
         MCsamples['01DrellYan'].append(TFile(DrellYan.RootFile,"READ"))         
         #MC Signal
         Sigsamples['2T2tb-mStop550mLSP1br50'].append(TFile(str(t550n1b50.RootFile),"READ"))
         Sigsamples['1T2tb-mStop325mLSP100br50'].append(TFile(str(t325n100b50.RootFile),"READ"))
         #Create the mcfiles dictionary for convencience
         #maps the name of the sample to the TFile
         mcfiles={}         
         for name,lists in MCsamples.iteritems():
             mcfiles[name]=lists[1]
         #
         print 'data ',Data.RootFile
         print 'DiLep ',DiLep.RootFile
         print 'OneLep ',OneLep.RootFile
         print 'WJets ',WJets.RootFile
         print 'Rare ',Rare.RootFile
         print 'QCD ',QCD.RootFile
         print 'DrellYan ',DrellYan.RootFile

         
         print 'T2tb: m_stop = 550; m_LSP = 1',t550n1b50.RootFile
         print 'T2tb: m_stop = 325; m_LSP = 100',t325n100b50.RootFile
         
         #raw_input()

         #
         #
         #===================================================
         #
         #
         #
         #
         #
         #==================================================
         #CUTNAMES, which are the directories with the plots
         #==================================================    
         cutnames=[]

         for SR in ["26", "23", "57", "41", "24", "17", "42", "53", "66", "52"]:
             for lep in ["El", "Mu", "ElAndMu"]:
                 #for CR in ["Preselection", "SearchRegionPreIsoTrackVeto", "SearchRegionPostIsoTrackVeto", "CR1", "CR4", "CR5"]:
                 for CR in ["SearchRegionPostIsoTrackVeto", "CR1", "CR5"]:   
                     cutnames.append(SR+"/"+lep+"-"+CR)



         #====Form the empty dictionaries for the histograms
         dataDict={}
         #the key is the name of the sample, 
         datakey=Data.SampleName
         #the itemx is a list
         dataDict[datakey]=[]
         #the first entry of the list will have to be the histogram pointer
         #I put 'nil' but I will overwrite it with the pointer later on.
         dataDict[datakey].append('nil')
         #the second one is a dictionary of properties, than CAN be empty
         myoptdict={}
         dataDict[datakey].append(myoptdict)
         #
         #
         #
         #====Form the empty dictionaries for the MC histograms
         #CONSTRUCT THE EMPTY BACKGROUND DICTIONARY
         mcDict={}
         for sampleskey,sampleslist in MCsamples.iteritems():

             #for the MC, it has to be the same key
             #as for the MCsample dictionary.
             mckey=sampleskey
             #same thing, I put 'nil' in here but will put the right
             #pointer afterwards
             mcDict[mckey]=['nil']
             myoptdict={}
             #
             mcDict[mckey].append(copy.deepcopy(myoptdict))
         #   
         #CONSTRUCT THE EMPTY SIGNAL DICTIONARY, OPTIONAL
         sigDict={}
         for sampleskey, sampleslist in Sigsamples.iteritems():

             #
             #
             sigDict[sampleskey]=['nil']
             myoptdict={}
             #
             sigDict[sampleskey].append(myoptdict)


         #==============================================
         #LOOP OVER THE DIRECTORIES AND HISTOGRAMS
         #==============================================
         for dir in cutnames:

             # 
             #ITERATE OVER ALL THE HISTOS IN THIS DIRECTORY
             #
             myDataDirName=""+dir
             myDataDir=datafile.Get(myDataDirName)
             #
             if (str(myDataDir).find("nil") != -1):
                 raise NameError("the directory ",dir, "was not found")
             #
             nextkey=TIter(myDataDir.GetListOfKeys());   
             key=nextkey();        
             #
             #
             #==============ITERATE OVER HISTOS================
             while (str(key).find("nil") == -1): #this means that is not a null pointer

                 #===============GET THE DATA HISTO=============================
                 datahist = myDataDir.Get(key.GetName())
                 if (datahist.GetEntries() == 0):
                     key=nextkey()
                     continue
                 #UPDATE THE DICTIONARY FOR THE DATA
                 for datakey,dataList in dataDict.iteritems():
                     #PUTTING IN A HISTO CLONE
                     dataList[0]=copy.deepcopy(datahist.Clone())
                     #
                     #Get the style of the histo from the
                     #standard style and put it in the dictionary
                     dataList[1]=Style_DataHistogram()
                     #

                 #==============GET THE MONTE CARLO HISTOGRAMS=================
                 #=============================================================                 
                 # Bkg
                 #=============================================================
                 for mckey,mclist in MCsamples.iteritems():
                     #key is a string with the name of the sample
                     #mclist[0] is the datasample
                     #mclist[1] is the opened root file
                     mcfile=mclist[1]
                     #
                     #
                     #print str(myDataDir.GetPath())
                     dataPath=str(myDataDir.GetPath())
                     dataPath_index=dataPath.rfind(':')
                     dataPath=dataPath[dataPath_index+2 : ]
                     #
                     #print dataPath+"/"+key.GetName()
                     dummyhist = mcfile.Get(dataPath+"/"+key.GetName())
                     #
                     #
                     if str(dummyhist).find('nil') != -1:
                         print key.GetName(), 'was not found ',keyfile
                         raise NameError('ERROR, the histo was not found')
                     #====UPDATE THE MC Bkg DICTIONARY WITH THE HISTO POINTER CLONE
                     mcDict[mckey][0]=copy.deepcopy(dummyhist.Clone())                     
                     #SET THE DICTIONARY OPTIONS FROM THE STANDARDS
                     #GetColor is defined in SetRootObjectsOptions
                     mcDict[mckey][1]['lineColor']=GetColor(mckey)
                     mcDict[mckey][1]['fillColor']=GetColor(mckey)
                     mcDict[mckey][1]['NoErrors']=True;
                     #
                     
                 #=============================================================                 
                 # Signal
                 #=============================================================
                     
                 for sigkey,siglist in Sigsamples.iteritems():

                     sigfile=siglist[1]

                     dataPath=str(myDataDir.GetPath())
                     dataPath_index=dataPath.rfind(':')
                     dataPath=dataPath[dataPath_index+2 : ]
                     #
                     dummyhist = sigfile.Get(dataPath+"/"+key.GetName())
                     #
                     #
                     if str(dummyhist).find('nil') != -1:
                         print key.GetName(), 'was not found ', sigkey
                         raise NameError('ERROR, the histo was not found')
                     #===UPDATE THE MC Signal DICTIONARY WITH THE HISTO POINTER CLONE
                     sigDict[sigkey][0]=copy.deepcopy(dummyhist.Clone())                     
                     #SET THE DICTIONARY OPTIONS FROM THE STANDARDS
                     #GetColor is defined in SetRootObjectsOptions
                     sigDict[sigkey][0].SetDirectory(0)
                     sigDict[sigkey][1]=Style_SignalHistogram(sigkey)
                     sigDict[sigkey][1]['fillColor']=GetColor(sigkey)

                 #
                 #
                 #NOW CONSTRUCT THE PLOT WITH THESE HISTOGRAMS
                 #
                 #
                 #===FIRST CONSTUCT THE PROPERTIES OF THE PLOTS
                 #get some standard properties
                 PlotProperties=Style_ControlPlots(Data.SampleLumi/1000.)
                 PlotProperties['xtitle_ratioplot']=GetUnits(datahist.GetName())
                 #
                 #if you want the MC not stacked, uncomment the next line
                 #PlotProperties['drawingOptions_stack']='NoStack'
                 #
                 #
                 #
                 #
                 
                 #
                 #SET THE NAME OF THE PLOT AS THE WHOLE PATH OF THE PLOT
                 #some crazy operations on strings to have the right name
                 #essentially a combination of string slicing with the operator[a:b]
                 #and the method 'find()'
                 c=myDataDir.GetPath()
                 d=c[c.find(':')+1 : ] #this means that d equals a slice of c from where ":"
                 #is found +1 position,
                 #until the end of the string.
                 currentdir=d[ c.find('"') +2 : ]
                 currentdir=topPlotsDir+'/'+currentdir
                 #print 'currentdir =',currentdir
                 #raw_input()
                 
                 #outPlotName=currentdir+'/'+datahist.GetName()
                 #outPlotName=outPlotName.replace(Data.SampleName+'/ControlPlots/','')
                 #outPlotName=outPlotName.replace('/','_')
                 outPlotName=datahist.GetName()
                 PlotProperties['outPlotName']=outPlotName
                 #print 'outPlotName is ', outPlotName
                 #raw_input()
                 #PlotProperties['outPlotName']=outPlotName                 
                 #======================================
                 #
                 #
                 #
                 #
                 #
                 #      CREATE THE PLOT
                 #
                 #
                 #
                 #
                 #
                 #
                 #========================================
                 #Lets create the input list for the class that we are using
                 InputList=['StackWithRatio',PlotProperties,dataDict,mcDict,sigDict]
                 newPlot=DataBkgSigPlot(InputList)
                 #
                 #      DRAW IT
                 newPlot.Draw()
                 #
                 #
                 #     PUT THE PLOTS IN THE CORRESPONDING DIR:
                 CheckDirectory(currentdir)
                 #
                 commands.getoutput('mv '+outPlotName+'.*'+' '+currentdir)                 
                 #
                 #Run the html template
                 plotdescription=datahist.GetTitle()
                 #
                 #
                 #CREATE AN HTML FILE
                 CreateHTML(outPlotName,plotdescription,currentdir)
                 #        TO THE NEXT PLOT
                 myDataDir.cd()
                 #
                 key=nextkey()




#END
