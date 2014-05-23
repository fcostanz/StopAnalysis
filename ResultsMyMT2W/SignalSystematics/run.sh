for ((SR=0; SR<9; SR++))
  do
  rm -f script${SR}.sh
  touch script${SR}.sh

  echo "#!/bin/zsh" >> script${SR}.sh
  echo pwd >> script${SR}.sh
  dir=`pwd`
  echo "root -l -q '${dir}/MakeHistos.C+($SR)'" >> script${SR}.sh

  qsub -cwd -V -l h_rt=2:30:00 -l h_vmem=1700M -P cms ./script${SR}.sh
done
