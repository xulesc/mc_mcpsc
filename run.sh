#!/bin/sh

PDB_DIR='../../../pdb_chew_kedem'
TMP_DIR='tmp'
TEST_PDB_DIR='../../../pdb_test'

PDB_DIR_CK='../pdb_chew_kedem/'
CM_DIR_CK='../contact_maps_chew_kedem/'
PDB_DIR_RS='../pdb_rost_sander/'
CM_DIR_RS='../contact_maps_rost_sander/'
LOG_CK='tmp.log.pctest.ck'
LOG_RS='tmp.log.pctest.rs'
UNIPROT_DUMP_DIR='uniprot_dump'
PDB_DOM_LIST='pdb_dom_list'
UNIPROT_URL='http://www.uniprot.org/uniprot/?query='

BLAST_PRG='/home/xule/downloads/software/ncbi-blast-2.2.28+-src/c++/GCC480-Debug64/bin/blastp'

run_ce()
{
	./ceskel - $PDB_DIR/$1 - $PDB_DIR/$2 - $TMP_DIR	
}

ce_batch_dir()
{
	for file1 in `ls $PDB_DIR`; do
		for file2 in `ls $PDB_DIR`; do
			echo "$file1 vs $file2"
			run_ce $file1 $file2
		done;
	done;
}

ce_batch_dir2()
{
	./ceskel $TMP_DIR $PDB_DIR
}

mcpsc_pc_test()
{
	./mcpsc $1 $2 $3> $4
	grep -E 'CE|TMa|USM ' $4 | less -S	
}

run_pairwise_blast()
{
	perf=`$BLAST_PRG -threshold 5 -subject $1 -query $2 | ./get_blast_scores.py`
	echo "$1 $2 $perf" >> $3
}

blast_pdb()
{
	echo "protein1 protein2 score expect identity positives gaps" > blast_log_$2
	for file1 in `ls $1*.seq`; do
		for file2 in `ls $1*.seq`; do
			run_pairwise_blast $file1 $file2 blast_log_$2		
		done;
	done;
}


get_uniprot_entries_for_pdb()
{
	[ -e $1 ] || mkdir $1
	
	[ -e delme ] || rm delme
	cat *.lq.dat | awk '{print $1}' > delme
	cat *.lq.dat | awk '{print $2}' >> delme
	cat *.uq.dat | awk '{print $1}' >> delme
	cat *.uq.dat | awk '{print $2}' >> delme	
	sort -u delme > $PDB_DOM_LIST
	rm delme
	
	for dom in `cat $PDB_DOM_LIST`; do
		pdb_name=`echo $dom | awk -F'.' '{print $1}'`
		wget $UNIPROT_URL$pdb_name -O $UNIPROT_DUMP_DIR/$pdb_name
	done;
}

get_uniprot_summary_data_for_pdb()
{
	../make_data.py $PDB_DOM_LIST $UNIPROT_DUMP_DIR
}

get_uniprot_features_for_pdb()
{
	../make_data.py $PDB_DOM_LIST.uniprot.dat
}

gen_features_file()
{
	../make_data.py $1 $2 $PDB_DOM_LIST.uniprot.features.dat > full_features.$1
}

##### processing #####
#
##Data from log
mcpsc_pc_test $TMP_DIR $PDB_DIR_CK $CM_DIR_CK $LOG_CK
#mcpsc_pc_test $TMP_DIR $PDB_DIR_RS $CM_DIR_RS $LOG_RS
##Blast a pair if sequences
#run_pairwise_blast ../pdb_chew_kedem/1aa9.seq ../pdb_chew_kedem/1ash.seq
##Make sequence files for all domains 
#blast_pdb $PDB_DIR_CK ck
#blast_pdb $PDB_DIR_RS rs
##Download results of pdb id search on uniprot
#get_uniprot_entries_for_pdb $UNIPROT_DUMP_DIR
##Extract uniprot ids for pdb files
#get_uniprot_summary_data_for_pdb
#get_uniprot_features_for_pdb
##Make features files
#gen_features_file tm.log.pctest.ck.lq.dat blast_log_ck 
#gen_features_file tm.log.pctest.ck.uq.dat blast_log_ck
#gen_features_file tm.log.pctest.rs.lq.dat blast_log_rs
#gen_features_file tm.log.pctest.rs.uq.dat blast_log_rs
#gen_features_file ce.log.pctest.ck.lq.dat blast_log_ck
#gen_features_file ce.log.pctest.ck.uq.dat blast_log_ck
#gen_features_file ce.log.pctest.rs.lq.dat blast_log_rs 
#gen_features_file ce.log.pctest.rs.uq.dat blast_log_rs
##
#mv full_features.* full_feature_data
#sed -e "s/^/'0',/g" full_feature_data/full_features.tm.log.pctest.ck.lq.dat > full_feature_data/full_features.tm.ck.lq
#sed -e "s/^/'1',/g" full_feature_data/full_features.tm.log.pctest.ck.uq.dat > full_feature_data/full_features.tm.ck.uq
#cat full_feature_data/full_features.tm.ck.lq full_feature_data/full_features.tm.ck.uq > full_features.tm.ck

##



