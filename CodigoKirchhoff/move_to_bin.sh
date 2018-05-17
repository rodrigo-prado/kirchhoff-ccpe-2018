#!/bin/bash

HOME_DIR=~/dev/kirchhoff-project/CodigoKirchhoff
DEST_DIR=~/dev/kirchhoff-project/MIG-TEMPO-OK/SYNT-DATA/bin
PROJ_DIR=Sequential/

for i in {1..4};
do
	case $i in
	"1")
		PROJ_DIR=Sequential
		cp -f $HOME_DIR/${PROJ_DIR}_Auto_Optimized/suktmig2d_${PROJ_DIR}_Auto_Optimized $DEST_DIR/
	;;
	"2")
		PROJ_DIR=OpenMP
	;;
	"3")
		PROJ_DIR=MPI	
		cp -f $HOME_DIR/${PROJ_DIR}_Auto_Optimized/suktmig2d_${PROJ_DIR}_Auto_Optimized $DEST_DIR/
	;;
	"4")
		PROJ_DIR=MPI+OpenMP
	;;
	esac
	cp -f $HOME_DIR/$PROJ_DIR/suktmig2d_${PROJ_DIR}_Scalar $DEST_DIR/
	cp -f $HOME_DIR/$PROJ_DIR/suktmig2d_${PROJ_DIR}_Auto $DEST_DIR/
	cp -f $HOME_DIR/${PROJ_DIR}_Manual/suktmig2d_${PROJ_DIR}_Manual $DEST_DIR/
done

