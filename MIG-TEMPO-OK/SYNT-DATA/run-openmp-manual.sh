if [ -z "$1" ] ; then
        nfc=16
else
        nfc=$1
fi
if [ -z "$2" ] ; then
        fwidth=5
else
        fwidth=$2
fi
if [ -z "$3" ] ; then
        nt=6
else
        nt=$3
fi
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
#prog=~/dev/kirchhoff/CodigoKirchhoff/OpenMP_Manual/suktmig2d_OpenMP_manual
prog=./bin/suktmig2d_OpenMP_Manual

export OMP_NUM_THREADS=$nt

$prog <$infile nfc=$nfc fwidth=$fwidth dx=$dx intoff=$intoff firstcdp=$firstcdp lastcdp=$lastcdp nvelcdp=$nvelcdp offmin=$offmin offmax=$offmax vfile=$vfile > $outfile

