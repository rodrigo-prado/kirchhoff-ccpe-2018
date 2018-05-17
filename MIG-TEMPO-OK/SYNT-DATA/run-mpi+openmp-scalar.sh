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
        np=3
else
        np=$(($3+1))
fi
if [ -z "$4" ] ; then
	nt=3
else
	nt=$4
fi
firstcdp=1
lastcdp=318
nvelcdp=318
offmin=0
offmax=1260
infile=synt-off.su
outfile=saida-synt-mpi+openmp-scalar.su
input_file=./$infile
output_file=./$outfile
vfile=velxz.bin
dx=10
intoff=20
#prog=~/dev/kirchhoff/CodigoKirchhoff/MPI+OpenMP/suktmig2d_MPI+OpenMP
prog=./bin/suktmig2d_MPI+OpenMP_Scalar

export OMP_NUM_THREADS=$nt

/usr/bin/mpirun --mca btl ^openib -np $np $prog nfc=$nfc fwidth=$fwidth dx=$dx intoff=$intoff firstcdp=$firstcdp lastcdp=$lastcdp nvelcdp=$nvelcdp offmin=$offmin offmax=$offmax vfile=$vfile input_file=$input_file output_file=$output_file
