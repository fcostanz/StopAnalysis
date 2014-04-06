#!bin/zsh
for ((i=0; i<11; i++))
do
    for iSR in 26 23 57 41 24 17 42 53 66 52
    do
	root -l -q 'MakeHistos.C+('$i','$iSR')'
    done
done

cp MakeHistos/*.root MtPeakReweighting/.

for iSR in 26 23 57 41 24 17 42 53 66 52
do
    root -l -q 'MtPeakReweighting.C+('$iSR')'
done

for iSR in 26 23 57 41 24 17 42 53 66 52
do    
    for lep in ElAndMu
    do
	for Cut in SearchRegionPostIsoTrackVeto CR1 CR4 CR5
	do
	    root -l -q 'Rebinning.C+("'${iSR}/${lep}-${Cut}'")'
	done
    done
done


. ./interface.sh

python newControlPlots.py
