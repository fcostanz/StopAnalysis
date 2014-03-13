import sys, os
import commands
import optparse
import datetime

###########################################################################

__version__ = "0.0"
__author__ = "Francesco Costanza (francesco.costanza@desy.de)"


###########################################################################

class copyToWeb(object):

    def __init__(self):
        
	self.version=__version__
	
	cmd_line_opts = sys.argv[1:]
	#print cmd_line_opts
	self.cmd_line_opts = cmd_line_opts

        #default value command-line options
        
        self.tag=''
        self.timeTag=datetime.datetime.today().strftime("%Y-%m-%d-%Hh%Mm")
        
        self.DS="SingleMu"
        self.controlPlotsDir="./"+self.DS+"/ControlPlots/"
        self.webDir="/afs/desy.de/user/f/fcost/www/RA4b/"
        self.storageDir="/scratch/hh/dust/naf/cms/user/costanza/ControlPlots/"
        
##         self.cuts=["Before_CutFlow","Triggers","Jet_Cuts","Wide_Electrons","HT","MET","HT400","MET200","btag1"]
##         self.cutsLabels={"Before_CutFlow" : "Before any cut",
##                          "Triggers"       : "After triggers",
##                          "Jet_Cuts"       : "After jets cuts",
##                          "Wide_Electrons" : "After single lepton requirements",
##                          "HT"             : "After HT cut",
##                          "MET"            : "After MET cut",
##                          "HT400"          : "HT > 400",
##                          "MET200"         : "MET > 200",
##                          "btag1"          : "at least one btagged jet" }
        self.cuts=["1+BJets_MET0",
                   "1+BJets_MET50",
                   "1+BJets_MET100",
                   "1+BJets_MET150"
                   ]
        self.cutsLabels={"1+BJets_MET0" : "1 BJet, MET > 0",
                         "1+BJets_MET50" : "1 BJet, MET > 50",
                         "1+BJets_MET100" : "1 BJet, MET > 100",
                         "1+BJets_MET150" : "1 BJet,MET > 150"
                         }
        self.main_index=''' '''

        self.versions_index='''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>PLOT COLLECTION</title>
<link href="../styles/site.css" rel="stylesheet">
<style>

</style>
<script src="../styles/jquery-1.6.3.min.js"></script>


</head>

<!--==========================================-->
<!-- BODY                                     -->
<!--==========================================-->

<body>

  <div class="main">
<h1>Control Plots for the sample $DATASET </h1>
  </div>

<!--ADDHERE-->

</body>
</head>
'''

        self.cuts_index='''<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<title>PLOT COLLECTION</title>
<link href="../../../../styles/site.css" rel="stylesheet">
<style>

</style>
<script src="../../../../styles/jquery-1.6.3.min.js"></script>


</head>

<!--==========================================-->
<!-- BODY                                     -->
<!--==========================================-->

<body>

  <div class="main">
<h1>Control Plots for the sample $DATASET </h1>
  </div>

<!--ADDHERE-->

</body>
</head>
'''
    def parse_cmd_line_options(self):

        parser = optparse.OptionParser()        
        self.option_parser = parser
        
        # Option to specify the name (or a regexp) of the dataset(s)
        # to be used.
        parser.add_option("", "--tag",
                          help='Name of the tag to add to the end of the dir name (default="")',
                          action="callback",
                          callback=self.option_handler_input_tag,
                          type="string",
                          metavar="TAG")
        parser.add_option("", "--dataset",
                          help='Name of the dataset you want to load the plots from (default="SingleMu")',
                          action="callback",
                          callback=self.option_handler_input_dataset,
                          type="string",
                          metavar="DS")
        parser.add_option("", "--controlPlotsDir",
                          help='Path to the dir where the control Plots are corrently stored (default="'+self.controlPlotsDir+'")',
                          action="callback",
                          callback=self.option_handler_input_controlPlotsDir,
                          type="string",
                          metavar="CONTROLPLOTSDIR")
        parser.add_option("", "--webDir",
                          help='Path to the directory where the dir tree for the web page starts (default="/afs/desy.de/user/f/fcost/www/")',
                          action="callback",
                          callback=self.option_handler_input_webDir,
                          type="string",
                          metavar="WEBDIR")
        parser.add_option("", "--storageDir",
                          help='Path to the directory where to store all plots (not only the ones for the web page)  \n (default="/scratch/hh/dust/naf/cms/user/costanza/ControlPlots")',
                          action="callback",
                          callback=self.option_handler_input_webDir,
                          type="string",
                          metavar="STORAGEDIR")
        
	parser.set_defaults()
        (self.options, self.args) = parser.parse_args(self.cmd_line_opts)
        
    def option_handler_input_tag(self, option, opt_str, value, parser):
        self.tag=value
        
    def option_handler_input_dataset(self, option, opt_str, value, parser):
        self.DS=value
        self.controlPlotsDir="./"+self.DS+"/ControlPlots/"
        
    def option_handler_input_controlPlotsDir(self, option, opt_str, value, parser):
        self.controlPlotsDir=value
        
    def option_handler_input_webDir(self, option, opt_str, value, parser):
        self.webDir=value

    def option_handler_input_storageDir(self, option, opt_str, value, parser):
        self.storageDir=value
        
    def run(self):
        self.parse_cmd_line_options()


        ##create dir structure for the storage
        
        self.storage=self.storageDir.replace("//","/")
        
        self.storageDir+="/"+self.DS+"/"
        self.storage=self.storageDir.replace("//","/")
        
        self.estimation=commands.getoutput("ls "+self.controlPlotsDir)
        self.storageDir+=self.estimation+"/"
        self.controlPlotsDir+=self.estimation+"/"
        
        self.tail=commands.getoutput("ls "+self.controlPlotsDir)
        self.storageDir+=self.tail+"/"
        self.controlPlotsDir+=self.tail+"/ControlPlots"
        
        self.fullTag=self.timeTag
        if (not self.tag == "" ):
            self.fullTag+="_"+self.tag
            
        self.storageDir+="/"+self.fullTag
        self.storageDir=self.storageDir.replace("//","/")
        os.makedirs( self.storageDir)
        
        command="cp -r "+self.controlPlotsDir+"/* "+self.storageDir+"/."
        command=command.replace("//","/")
        print command
        out=commands.getoutput(command)
        #print out


        self.webDir+=self.DS+"/"
	print self.webDir
        if (not os.path.exists(self.webDir)):
            os.makedirs(self.webDir)
        #create html files:

        if (not os.path.isfile(self.webDir+'/versions_index.html')):
            f=open(self.webDir+'/versions_index.html','w')
            f.write(self.versions_index)
            f.close()
        
        newFile=''
        f=open(self.webDir+'/versions_index.html','r')
        for line in f.readlines():
            if line.find("$DATASET") > 0:
                line=line.replace("$DATASET",self.DS)
            if line.find("ADDHERE") > 0:
                line+='<p><a href= "./'+self.estimation+'/'+self.tail+'/'+self.fullTag+'/cuts_index.html">'+self.DS+'  '+self.estimation+'  '+self.tail+'  '+self.fullTag+'</a></p>\n'
            newFile+=line
        f.close()
        f=open(self.webDir+'/versions_index.html','w')
        f.write(newFile)
        f.close()

        ##create dir structure for the web page
        
        self.webDir+=self.estimation+"/"
        self.webDir+=self.tail+"/"
        self.webDir+=self.fullTag+"/"
        os.makedirs(self.webDir)

        if (not os.path.isfile(self.webDir+'/cuts_index.html')):
            f=open(self.webDir+'/cuts_index.html','w')
            f.write(self.cuts_index)
            f.close()
                
        newFile=''
        f=open(self.webDir+'/cuts_index.html','r')
        for line in f.readlines():
            if line.find("$DATASET") > 0:
                line=line.replace("$DATASET",self.DS)
            if line.find("ADDHERE") > 0:
                for cut in self.cuts:
                    line+='<p><a href= "./'+cut+'/cplots_index.html">'+self.cutsLabels[cut]+'</a></p>\n'
            newFile+=line
        f.close()
        f=open(self.webDir+'/cuts_index.html','w')
        f.write(newFile)
        f.close()

        ##cp plots
        for cut in self.cuts:
            os.makedirs(self.webDir+"/"+cut)
            command="cp "+self.controlPlotsDir+"/"+cut+"/*.gif "+self.webDir+"/"+cut+"/."
            out=commands.getoutput(command)
            command="cp myDummy2.html "+self.webDir+"/"+cut+"/cplots_index.html"           
            
            #command="cp "+self.controlPlotsDir+"/"+cut+"/cplots_"+cut+".html "+self.webDir+"/"+cut+"/cplots_index.html"
            out=commands.getoutput(command)
        
            

        
        



        
###########################################################################
## Main entry point.
###########################################################################

if __name__ == "__main__":
    "Main entry point for copyToWeb."

    copyToWeb().run()

    # Done.

###########################################################################
