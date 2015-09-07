#!/bin/sh

PRG="mcpsc"
TMP="$HOME/workspace_git/mc_mcpsc/tmp/"
## CK34
PDB_DIR="$HOME/workspace_git/data/pdb_chew_kedem/"
CMP_DIR="$HOME/workspace_git/data/contact_maps_chew_kedem/"
DOM_FILES="$HOME/workspace_git/data/test_dom_files"
JOB_PAIRS="$HOME/workspace_git/data/test_jobs"
CM_FILES="$HOME/workspace_git/data/test_cm_files"
## RS 119
#PDB_DIR="$HOME/workspace_git/data/pdb_rost_sander/"
#CMP_DIR="$HOME/workspace_git/data/contact_maps_rost_sander/"
LOG="log"
CMD="./$PRG $TMP $PDB_DIR $CMP_DIR $DOM_FILES $JOB_PAIRS $CM_FILES"

#echo "#CMD=$CMD"
#$CMD
#for N in `seq 1 8`; do
#    export OMP_NUM_THREADS=$N; $CMD > $LOG.$N
#    grep 'Total time (msec)' $LOG.$N | awk '{print $NF}' | python -c 'import sys;nc=sys.argv[1];(usm,tm,ce,ignore)=sys.stdin.read().split("\n"); print "%s %s %s %s %d" %(nc,usm,tm,ce,(int(usm)+int(tm)+int(ce)))' $N
#done
#./$PRG /home/anuj/workspace/mc_mcpsc/tmp/ /home/anuj/Downloads/pdb40d/ /home/anuj/Downloads/contact_maps_pdb40d/ ./test_dom_files ./test_jobs ./test_cm_files > pdb40d.psc.log
#export OMP_NUM_THREADS=1
#./$PRG /home/anuj/workspace/mc_mcpsc/tmp/ /home/anuj/Downloads/pdb40d/ /home/anuj/Downloads/contact_maps_pdb40d/ ./dom_list ./job_pairs ./cm_list > pdb40d.psc.log
./$PRG /home/anuj/workspace/mc_mcpsc/tmp/ /home/anuj/Downloads/astral_40p/ /home/anuj/Downloads/contact_maps_astral_40p/ ./dom_list1 ./job_pairs1 ./cm_list1 > pdb40d.psc.log.1


    
    
