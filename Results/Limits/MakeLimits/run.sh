data="InputForLimits.root"
template="template.card"
script="FillHisto.py"
basedir="/nfs/dust/cms/user/fcost/workdir/StopAnalysis/Results/Limits/MakeLimits/"

for ((iStop=0; iStop<29; iStop++)) #29
  do
  mStop=`expr $iStop \* 25`
  mStop=`expr $mStop + 100`
  
  dirName=mStop${mStop}-mLSP1
  rm -rf ${dirName}
  mkdir ${dirName}

  cp ${data} ${dirName}/.
  cp ${template} ${dirName}/.
  cp lands.exe expected.sh observed.sh macro.C fitRvsCLs_C.so ${dirName}/.
  cp ${script} ${dirName}/.

  cd ${dirName}

  scriptName=script${mStop}-1.sh

  rm -f ${scriptName}
  touch ${scriptName}
    
  echo "#!/bin/zsh" >> ${scriptName}
  echo pwd >> ${scriptName}
  dir=`pwd`
  echo "python "${script}" "${mStop}" 1" >> ${scriptName}
  echo "rm my.lock" >> ${scriptName}
  
  touch my.lock

  qsub -cwd -V -l h_rt=5:59:00 -l h_vmem=2000M -P cms ./${scriptName}
  cd ${basedir}

  for ((iLSP=1; iLSP<15; iLSP++)) #15
    do
    mLSP=`expr $iLSP \* 25`

    deltaM=`expr $mStop - $mLSP`
    if [[ ${deltaM} -lt 90 ]]
	then
	continue
    fi

    dirName=mStop${mStop}-mLSP${mLSP}
    rm -rf ${dirName}
    mkdir ${dirName}

    cp ${data} ${dirName}/.
    cp ${template} ${dirName}/.
    cp lands.exe lands.sh macro.C fitRvsCLs_C.so ${dirName}/.
    cp ${script} ${dirName}/.
    
    cd ${dirName}

    scriptName=script${mStop}-${mLSP}.sh

    rm -f ${scriptName}
    touch ${scriptName}
    
    echo "#!/bin/zsh" >> ${scriptName}
    echo pwd >> ${scriptName}
    dir=`pwd`
    echo "python "${script}" "${mStop}" "${mLSP} >> ${scriptName}
    echo "rm my.lock" >> ${scriptName}

    touch my.lock

    qsub -cwd -V -l h_rt=5:59:00 -l h_vmem=2000M -P cms ./${scriptName}
    cd ${basedir}    
  done
done

