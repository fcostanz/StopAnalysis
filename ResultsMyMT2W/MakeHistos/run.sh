#!bin/zsh
root -l -q './MakeHistos.C+( 2, 0)'

for SR in 0 1 2 3 4 5 6 7 8 9 10 11
  do
  mkdir $SR
  for ((i=0; i<8; i++))
    do
    rm -f script${i}-${SR}.sh
    touch script${i}-${SR}.sh
    
    echo "#!/bin/zsh" >> script${i}-${SR}.sh
    echo pwd >> script${i}-${SR}.sh
    dir=`pwd`
    echo "root -l -q '${dir}/MakeHistos.C+($i, $SR)'" >> script${i}-${SR}.sh
    
    qsub -cwd -V -l h_rt=2:30:00 -l h_vmem=1700M -P af-cms ./script${i}-${SR}.sh
  done
done

#cp MakeHistos/*.root MtPeakReweighting/.

#for iSR in 26 23 57 41 24 17 42 53 66 52
#do
#    root -l -q 'MtPeakReweighting.C+('$iSR')'
#done

#for iSR in 26 23 57 41 24 17 42 53 66 52
#do    
#    for lep in ElAndMu
#    do
#	for Cut in SearchRegionPostIsoTrackVeto CR1 CR4 CR5
#	do
#	    root -l -q 'Rebinning.C+("'${iSR}/${lep}-${Cut}'")'
#	done
#    done
#done


#. ./interface.sh

#python newControlPlots.py
