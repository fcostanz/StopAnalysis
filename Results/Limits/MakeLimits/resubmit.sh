data="Optimization-Coverage.root"
template="template.card"
script="FillHisto.py"
basedir="/nfs/dust/cms/user/fcost/workdir/StopAnalysis/Results/ExpectedLimits/MakeLimits/"

for lock in `ls */my.lock`
do
  echo ${lock}
  dir=`echo ${lock} | sed 's/\/my.lock//g'`
  cd ${dir}
  
  tmp=`echo ${dir} | sed 's/mStop//g' | sed 's/mLSP//g'`
  
  scriptName=script${tmp}.sh
  qsub -cwd -V -l h_rt=24:59:00 -l h_vmem=2000M -P cms ./${scriptName}

  cd ${basedir}
done
