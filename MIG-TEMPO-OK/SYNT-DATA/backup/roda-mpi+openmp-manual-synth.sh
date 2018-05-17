firstcdp=1
lastcdp=318
nvelcdp=318
offmin=0
offmax=1260
infile=synt-off.su
outfile=saida-synt-mpi-manual.su
vfile=velxz.bin
dx=10
intoff=20
prog=~/dev/kirchhoff/CodigoKirchhoff/MPI+OpenMP+AVX_Manual/suktmig2d_MPI+OpenMP_manual

export OMP_NUM_THREADS=2

/usr/bin/mpirun -np 5 $prog dx=$dx intoff=$intoff firstcdp=$firstcdp lastcdp=$lastcdp nvelcdp=$nvelcdp offmin=$offmin offmax=$offmax vfile=$vfile 
