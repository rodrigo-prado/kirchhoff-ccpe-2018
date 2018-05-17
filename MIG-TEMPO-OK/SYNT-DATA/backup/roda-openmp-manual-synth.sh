firstcdp=1
lastcdp=318
nvelcdp=318
offmin=0
offmax=1260
infile=synt-off.su
outfile=saida-synt-openmp-manual.su
vfile=velxz.bin
dx=10
intoff=20
prog=~/dev/kirchhoff/CodigoKirchhoff/OpenMP_Manual/suktmig2d_OpenMP_manual

export OMP_NUM_THREADS=6

$prog <$infile dx=$dx intoff=$intoff firstcdp=$firstcdp lastcdp=$lastcdp nvelcdp=$nvelcdp offmin=$offmin offmax=$offmax vfile=$vfile > $outfile

