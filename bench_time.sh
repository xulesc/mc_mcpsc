#!/bin/sh

PRG="mcpsc"
TMP="$HOME/local_wspace/mc_mcpsc/tmp/"
## CK34
PDB_DIR="$HOME/workspace_git/data/pdb_chew_kedem/"
CMP_DIR="$HOME/workspace_git/data/contact_maps_chew_kedem/"
## RS 119
#PDB_DIR="$HOME/workspace_git/data/pdb_rost_sander/"
#CMP_DIR="$HOME/workspace_git/data/contact_maps_rost_sander/"
LOG="log"
CMD="./$PRG $TMP $PDB_DIR $CMP_DIR"

echo "#CMD=$CMD"
for N in `seq 1 8`; do
    export OMP_NUM_THREADS=$N; $CMD > $LOG.$N
    grep 'Total time (msec)' $LOG.$N | awk '{print $NF}' | python -c 'import sys;nc=sys.argv[1];(usm,tm,ce,ignore)=sys.stdin.read().split("\n"); print "%s %s %s %s %d" %(nc,usm,tm,ce,(int(usm)+int(tm)+int(ce)))' $N
done

    
    