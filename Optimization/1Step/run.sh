#!bin/zsh
#for ((i=0; i<11; i++))
#do
#    root -l -q 'MakeTree.C+('$i')'
#done

rm -f SRs.txt
touch SRs.txt

for ((i=0; i<7; i++))
do
    for BR in "0.25" "0.50" "0.75" "1.0"
    do
	root -l -q 'bestFOM.C+('$i','$BR')' >> SRs.txt  
    done
done

