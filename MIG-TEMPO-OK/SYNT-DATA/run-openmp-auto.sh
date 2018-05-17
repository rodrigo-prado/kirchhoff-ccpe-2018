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
outfile=saida-synt-openmp-auto.su
vfile=velxz.bin
dx=10
intoff=20
#prog=~/dev/kirchhoff/CodigoKirchhoff/OpenMP/suktmig2d_OpenMP
prog=./bin/suktmig2d_OpenMP_Auto

export OMP_NUM_THREADS=$nt

$prog <$infile dx=$dx intoff=$intoff firstcdp=$firstcdp lastcdp=$lastcdp nvelcdp=$nvelcdp offmin=$offmin offmax=$offmax vfile=$vfile > $outfile

