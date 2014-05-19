#!bin/zsh
for ((iSample=0; iSample<8; iSample++))
  do
  root -l -q 'mergeFiles.C+('$iSample')'
done

hadd -f Data.root SingleElectron.root SingleMu.root

#cp MakeHistos/*.root MtPeakReweighting/.
#for iSR in 31 82 32 39 61 46 22 35 51 14 65 75 1 72 70 12 9 27
#do
#    root -l -q 'MtPeakReweighting.C+('$iSR')'
#done

#for iSR in 31 82 32 39 61 46 22 35 51 14 65 75 1 72 70 12 9 27
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
