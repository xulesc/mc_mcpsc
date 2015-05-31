///////////////////////////////////////////////////////////////////////////
//  Combinatorial Extension Algorithm for Protein Structure Alignment    //
//                                                                       //
//  Authors: I.Shindyalov & P.Bourne                                     //
///////////////////////////////////////////////////////////////////////////
/*
 
 Copyright (c)  1997-2000   The Regents of the University of California
 All Rights Reserved
 
 Permission to use, copy, modify and distribute any part of this CE
 software for educational, research and non-profit purposes, without fee,
 and without a written agreement is hereby granted, provided that the above
 copyright notice, this paragraph and the following three paragraphs appear
 in all copies.
 
 Those desiring to incorporate this CE Software into commercial products
 or use for commercial purposes should contact the Technology Transfer
 Office, University of California, San Diego, 9500 Gilman Drive, La Jolla,
 CA 92093-0910, Ph: (619) 534-5815, FAX: (619) 534-7345.
 
 IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
 DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING
 LOST PROFITS, ARISING OUT OF THE USE OF THIS CE SOFTWARE, EVEN IF THE
 UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH
 DAMAGE.
 
 THE CE SOFTWARE PROVIDED HEREIN IS ON AN "AS IS" BASIS, AND THE
 UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
 SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  THE UNIVERSITY OF
 CALIFORNIA MAKES NO REPRESENTATIONS AND EXTENDS NO WARRANTIES OF ANY KIND,
 EITHER IMPLIED OR EXPRESS, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT
 THE USE OF THE CE SOFTWARE WILL NOT INFRINGE ANY PATENT, TRADEMARK OR
 OTHER RIGHTS.
 
 */
#include "comm_block.h"
#include "ce.h"
#include "usm.hpp"
#include "tmalign_extern.h"
#include "pom/ipdb.h"
#include "rckskel.h"
#include <sys/timeb.h>
#include <stdio.h>
//
#include <fstream>
#include <sstream>
#include <string>
//
// little trick to allow the application to be called "RCCE_APP" under
// OpenMP, and "main" otherwise

//#define VERBOSE_PRINTS
//

#define CE_ALGO_TYPE 1
#define TMALIGN_ALGO_TYPE 0
////////////////////////////////////////////////////////////////////
// tmalign variables
char *seq_amino_acids, *amino_acid_chain_seq;
int *atom_index, *seq_length, zero = 0;
float *seq_coords;
//
typedef struct {
	int algo_type;
	protein_data *protein1, *protein2;
	int protein1_idx, protein2_idx;
} job_descriptor;

class SortData {
public:
	SortData(int p, int i) {
		product = p;
		index = i;
	}	
	int product;
	int index;
	bool operator<(const SortData &rhs) const { 
		return product > rhs.product;  // descending
	}	
};

int *job_indexes;
job_descriptor *job_pool;

void tmalign_load_data(int line_count, char **filenames, char *dir_name);
int tmalign_read_pdb_structure(char *fname, float *prot_backbone, char *ss,
		int *start, int *lenarray, char *seq);
int get_file_count(char *dir_name);
void ce_load_data(CE *ce, char *db_tmp_path, char *pdb_dir_name, int file_count,
		char **filenames, protein_data *proteins_data);
void read_params(char *db_tmp_path, char *pdb_dir_name, char **argv);

MasterToClientTransferBlock client_in_data;
ClientToMasterTransferBlock master_in_data;

int diff_timeb(timeb t1, timeb t2) {
	return (1000 * t2.time + t2.millitm) - (1000 * t1.time + t1.millitm);
}
//

#ifdef PC_TEST
void do_target_pc(int argc, char **argv);
#endif

#ifndef PC_TEST
void master_send_job_data(int job_idx, int ue_id);
int client_receive_job(timeb t1);
int master_receive_result(int id);
int check_ready(int id);
void terminate(int id);
void do_target_scc(int argc, char **argv);
#endif
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
template<class A> void NewArray(A *** array, int Narray1, int Narray2) {
	*array = new A*[Narray1];
	for (int i = 0; i < Narray1; i++)
		*(*array + i) = new A[Narray2];
}
//

void fill_dataset_details(char *dir_name, char **filenames) {
        DIR *mydir = opendir(dir_name);
        struct dirent *entry = NULL;
        int ccount = 0;

        // read file name and store
        while ((entry = readdir(mydir))) {
                if (strcmp(entry->d_name, ".") == 0
                                || strcmp(entry->d_name, "..") == 0) {
                        continue;
                }
                char *name;
                setText(&name, entry->d_name);
                filenames[ccount++] = name;
        }
        closedir(mydir);
}

int get_min_index(int size, int *data) {
	int lowest = 0;
	for(int i = 0; i < size; i++)
		if( data[i] <= data[lowest])
			lowest = i;
	return lowest;
}

std::vector<SortData> get_partitions(int split_size, int algo_type, int job_count, bool sort) {
	int part_sum[split_size];
	// init all sums
	for(int i = 0; i < split_size; i++) {
		part_sum[i] = 0;
	}
	// put job index and length product in vector
	std::vector<SortData> data;
	for(int i = 0; i < job_count; i++) {
		if(job_pool[i].algo_type == algo_type) {
			int product;
			if(algo_type == CE_ALGO_TYPE) {
				product = job_pool[i].protein1->nSe * job_pool[i].protein2->nSe;
			} else {
				product = seq_length[job_pool[i].protein1_idx] * seq_length[job_pool[i].protein2_idx];
			}
			data.push_back(SortData(i, product));
		}
	}
	// sort vector (reverse)
	if(sort)
		std::sort(data.begin(), data.end());
	// greedy partition
	//std::vector<int> partitions[split_size];
	//for(auto datum : data) {
	//	int index = get_min_index(split_size, part_sum);
	//	partitions[index].push_back(datum.index);
	//	part_sum[index] += datum.product;
	//}
	return data;
}

void dispatch(int ue_count, int job_index) {
	while(1) {
		for(int i = 1; i < ue_count; i++) {
			// if node is accepting send data to it
		}
	}
}

int main(int argc, char **argv) {
#ifdef PC_TEST
	do_target_pc(argc, argv);
	return 0;
#endif
#ifndef PC_TEST
	do_target_scc(argc, argv);
#endif
	return (0);
}

#ifndef PC_TEST
void terminate(int id) {
	if (id == MASTER_ID) {
		printf("not terminating master.\n");
		return;
	}
	printf("terminate ue %d\n", id);
	client_in_data.die = RCKSKEL_TRUE;
	RCCE_send((char *) &client_in_data, sizeof(MasterToClientTransferBlock),
			id);
}

void master_send_job_data(int job_idx, int ue_id) {
	printf("sending job: %d to node: %d, p1: %s, p2: %s\n", job_idx, ue_id,
			job_pool[job_idx].protein1->name, job_pool[job_idx].protein2->name);
	job_descriptor job = job_pool[job_idx];

	// send stay alive
	client_in_data.die = RCKSKEL_FALSE;
	client_in_data.job_index = job_idx;
	client_in_data.algo_type = job.algo_type;
	RCCE_send((char *) &client_in_data, sizeof(MasterToClientTransferBlock),
			ue_id);

	if (job.algo_type == CE_ALGO_TYPE) {
		protein_data *protein1 = job.protein1, *protein2 = job.protein2;
		// send name length
		int p1_length = strlen(protein1->name);
		int p2_length = strlen(protein2->name);
		RCCE_send((char *) &p1_length, sizeof(int), ue_id);
		RCCE_send((char *) &p2_length, sizeof(int), ue_id);
		// send names
		//printf("sending: %s\t%s\n",protein1->name, protein2->name);
		RCCE_send((char *) protein1->name, p1_length * sizeof(char), ue_id);
		RCCE_send((char *) protein2->name, p2_length * sizeof(char), ue_id);
		// send lengths
		int nse1 = protein1->nSe, nse2 = protein2->nSe;
		RCCE_send((char *) &nse1, sizeof(int), ue_id);
		RCCE_send((char *) &nse2, sizeof(int), ue_id);
		// send sequences
		//printf("sending: %s\t%s\n",protein1->seq, protein2->seq);
		RCCE_send((char *) protein1->seq, nse1 * sizeof(char), ue_id);
		RCCE_send((char *) protein2->seq, nse2 * sizeof(char), ue_id);
		// send iEnp
		int iEnp1 = protein1->iEnp1, iEnp2 = protein2->iEnp1;
		RCCE_send((char *) &iEnp1, sizeof(int), ue_id);
		RCCE_send((char *) &iEnp2, sizeof(int), ue_id);
		// send structures
		RCCE_send((char *) protein1->ca, nse1 * sizeof(XYZ), ue_id);
		RCCE_send((char *) protein2->ca, nse2 * sizeof(XYZ), ue_id);
	} else {
		int dom1_index = job.protein1_idx, dom2_index = job.protein2_idx;
		// amino_acid_chain_seq
#ifdef VERBOSE_PRINTS
		printf("Master: sending amino_acid_chain_seq\n");
#endif
		RCCE_send((char *) &amino_acid_chain_seq[5001 * dom1_index],
				5001 * sizeof(char), ue_id);
		RCCE_send((char *) &amino_acid_chain_seq[5001 * dom2_index],
				5001 * sizeof(char), ue_id);
		// atom_index
#ifdef VERBOSE_PRINTS
		printf("Master: sending atom_index\n");
#endif
		RCCE_send((char *) &atom_index[MAX_ATOMS * dom1_index],
				MAX_ATOMS * sizeof(int), ue_id);
		RCCE_send((char *) &atom_index[MAX_ATOMS * dom2_index],
				MAX_ATOMS * sizeof(int), ue_id);
		// seq1_amino_acids
#ifdef VERBOSE_PRINTS
		printf("Master: sending seq_amino_acids\n");
#endif
		RCCE_send((char *) &seq_amino_acids[3 * MAX_ATOMS * dom1_index],
				3 * MAX_ATOMS * sizeof(char), ue_id);
		RCCE_send((char *) &seq_amino_acids[3 * MAX_ATOMS * dom2_index],
				3 * MAX_ATOMS * sizeof(char), ue_id);
		// seq_coords
#ifdef VERBOSE_PRINTS
		printf("Master: sending seq_coords\n");
#endif
		RCCE_send((char *) &seq_coords[dom1_index * 3 * MAX_ATOMS],
				3 * MAX_ATOMS * sizeof(float), ue_id);
		RCCE_send((char *) &seq_coords[dom2_index * 3 * MAX_ATOMS],
				3 * MAX_ATOMS * sizeof(float), ue_id);
		// lengths
#ifdef VERBOSE_PRINTS
		printf("Master: sending seq_length\n");
		printf("l1: %d\n",seq_length[dom1_index]);
		printf("l2: %d\n",seq_length[dom2_index]);
#endif
		RCCE_send((char *) &seq_length[dom1_index], sizeof(int), ue_id);
		RCCE_send((char *) &seq_length[dom2_index], sizeof(int), ue_id);
#ifdef VERBOSE_PRINTS
		printf("sent all\n");
#endif
	}
}

int check_ready(int id) {
	int test_flag;
	READY_T ready_t;
	RCCE_recv_test((char *) &ready_t, sizeof(READY_T), id, &test_flag);
	if (test_flag == 0)
		return RCKSKEL_FALSE;
	return (ready_t.ready == 1) ? RCKSKEL_TRUE : RCKSKEL_FALSE;
}
int master_receive_result(int id) {
	int test;
	RCCE_recv_test((char *) &master_in_data,
			sizeof(ClientToMasterTransferBlock), id, &test);
	if (test == RCKSKEL_TRUE) {
		printf("Algo: %d, Protein1: %s, Protein2: %s, len1: %d, len2: %d, "
				"Rmsd: %.2fA, TM1: %.2f, TM2: %.2f, "
				"id: %d, sec_job_start: %ld, sec_result_send: %ld, "
				"processing_time_msec: %ld, data_collect_time_msec: %ld, "
				"idle_time_msec: %ld\n", master_in_data.algo_type,
				job_pool[master_in_data.job_index].protein1->name,
				job_pool[master_in_data.job_index].protein2->name,
				job_pool[master_in_data.job_index].protein1->nSe,
				job_pool[master_in_data.job_index].protein2->nSe,
				master_in_data.rmsd, master_in_data.tm1, master_in_data.tm2, id,
				master_in_data.seconds_at_job_start,
				master_in_data.seconds_at_result_send,
				master_in_data.processing_time,
				master_in_data.data_collect_time, master_in_data.idle_time);
	}
	return test;
}
int client_receive_job(timeb t1) {
	timeb t;

	RCCE_recv((char *) &client_in_data, sizeof(MasterToClientTransferBlock),
			MASTER_ID);
	ftime(&t);
	master_in_data.seconds_at_job_start = t.time;
	master_in_data.idle_time = diff_timeb(t1, t);
	// early exit if server is asking to leave
	if (client_in_data.die == RCKSKEL_TRUE) {
		return RCKSKEL_TRUE;
	}

	if (client_in_data.algo_type == CE_ALGO_TYPE) {
		int n1, n2, nse1, nse2, iEnp1, iEnp2;
		struct timeb tp1, tp2;

		ftime(&tp1);
		// receive name lengths
		RCCE_recv((char *) &n1, sizeof(int), MASTER_ID);
		RCCE_recv((char *) &n2, sizeof(int), MASTER_ID);
		// receive names
		char *name1 = (char*) malloc(n1 * sizeof(char));
		char *name2 = (char*) malloc(n2 * sizeof(char));
		RCCE_recv((char *) name1, n1 * sizeof(char), MASTER_ID);
		RCCE_recv((char *) name2, n2 * sizeof(char), MASTER_ID);
		//printf("received: %s\t%s\n",name1, name2);
		// receive structure lengths
		RCCE_recv((char *) &nse1, sizeof(int), MASTER_ID);
		RCCE_recv((char *) &nse2, sizeof(int), MASTER_ID);
		// receive sequences
		char *seq1 = (char*) malloc(nse1 * sizeof(char));
		char *seq2 = (char*) malloc(nse2 * sizeof(char));
		RCCE_recv((char *) seq1, nse1 * sizeof(char), MASTER_ID);
		RCCE_recv((char *) seq2, nse2 * sizeof(char), MASTER_ID);
		//printf("received: %s\t%s\n",seq1, seq2);
		// receive iEnps
		RCCE_recv((char *) &iEnp1, sizeof(int), MASTER_ID);
		RCCE_recv((char *) &iEnp2, sizeof(int), MASTER_ID);
		// receive structure data
		XYZ *ca1 = (XYZ *) malloc(nse1 * sizeof(XYZ));
		XYZ *ca2 = (XYZ *) malloc(nse2 * sizeof(XYZ));
		RCCE_recv((char *) ca1, nse1 * sizeof(XYZ), MASTER_ID);
		RCCE_recv((char *) ca2, nse2 * sizeof(XYZ), MASTER_ID);
		//printf("processing\n");
		ftime(&tp2);
		master_in_data.data_collect_time = diff_timeb(tp1, tp2);

		ftime(&tp1);

		double d_[20];
		int isPrint = 0;
		int lcmp;

		int *align_se1 = new int[nse1 + nse2];
		int *align_se2 = new int[nse1 + nse2];

		int aln_len, gaps;
		float rmsd;
		//double z = 0;
		double z = ce_1(name1, name2, ca1, ca2, nse1, nse2, align_se1,
				align_se2, lcmp, 8, 3.0, 4.0, d_, isPrint, &aln_len, &gaps,
				&rmsd);
		// print scores and alignments
		isPrint = 1;

		if (lcmp > 0) {
			int lsim = 0, lali = 0;
			for (int l = 0; l < lcmp; l++)
				if (align_se1[l] != -1 && align_se2[l] != -1) {
					if (seq1[align_se1[l]] == seq2[align_se2[l]])
						lsim++;
					lali++;
				}

			int lstep = 70;

			for (int l = 0; l < lcmp; l += lstep) {
				for (int ie = 0; ie < 2; ie++) {
					int *align_se = (ie == 0 ? align_se1 : align_se2);
					for (int l_ = l; l_ < l + lstep && l_ < lcmp; l_++)
						if (align_se[l_] != -1) {
							break;
						}
				}
			}

		}

		// free stuff
		delete[] ca1;
		delete[] ca2;
		delete[] seq1;
		delete[] seq2;
		delete[] name1;
		delete[] name2;
		if (nse1 + nse2 > 0) {
			delete[] align_se1;
			delete[] align_se2;
		}
		ftime(&tp2);

		master_in_data.aln_len = aln_len;
		master_in_data.job_index = client_in_data.job_index;
		master_in_data.len1 = nse1;
		master_in_data.len2 = nse2;
		master_in_data.z = z;
		master_in_data.rmsd = rmsd;
		master_in_data.gaps = gaps;
		master_in_data.processing_time = diff_timeb(tp1, tp2);
		master_in_data.algo_type = CE_ALGO_TYPE;
	} else {
		int i;
		char seq1_amino_acids[3 * MAX_ATOMS], seq2_amino_acids[3 * MAX_ATOMS];
		char amino_acid_chain_seq1[5001], amino_acid_chain_seq2[5001];
		backbone_ *backbone_1 = (backbone_ *) malloc(sizeof(backbone_));
		int atom_indices1[MAX_ATOMS], atom_indices2[MAX_ATOMS];
		int len1, len2;
		struct timeb tp1, tp2;
		ftime(&tp1);
		// now collect the batch of data
		RCCE_recv((char *) &amino_acid_chain_seq1[0], 5001 * sizeof(char),
				MASTER_ID);
		RCCE_recv((char *) &amino_acid_chain_seq2[0], 5001 * sizeof(char),
				MASTER_ID);
#ifdef VERBOSE_PRINTS
		printf("got sequences\n");
#endif
		// atom_index
		RCCE_recv((char *) &atom_indices1[0], MAX_ATOMS * sizeof(int),
				MASTER_ID);
		RCCE_recv((char *) &atom_indices2[0], MAX_ATOMS * sizeof(int),
				MASTER_ID);
#ifdef VERBOSE_PRINTS
		printf("got atom indices\n");
#endif
		// seq_amino_acids
		RCCE_recv((char *) &seq1_amino_acids[0], 3 * MAX_ATOMS * sizeof(char),
				MASTER_ID);
		RCCE_recv((char *) &seq2_amino_acids[0], 3 * MAX_ATOMS * sizeof(char),
				MASTER_ID);
#ifdef VERBOSE_PRINTS
		printf("got sequences\n");
#endif
		// seq_coords
		RCCE_recv((char *) &backbone_1->xa[0], 3 * MAX_ATOMS * sizeof(float),
				MASTER_ID);
		RCCE_recv((char *) &backbone_1->xa[3 * MAX_ATOMS],
				3 * MAX_ATOMS * sizeof(float), MASTER_ID);
#ifdef VERBOSE_PRINTS
		printf("got coordinates\n");
#endif
		// lengths
		RCCE_recv((char *) &len1, sizeof(int), MASTER_ID);
		RCCE_recv((char *) &len2, sizeof(int), MASTER_ID);
#ifdef VERBOSE_PRINTS
		printf("got lengths\n");
		printf("all data finished\n");
		printf("\n\n");
		printf("data received on UE\n");
		printf("%c",amino_acid_chain_seq1[i]);
		printf("\n");
#endif
		ftime(&tp2);
		master_in_data.data_collect_time = diff_timeb(tp1, tp2);

		//do tmalign
		ftime(&tp1);
		main_process(seq1_amino_acids, seq2_amino_acids, amino_acid_chain_seq1,
				amino_acid_chain_seq2, &master_in_data, backbone_1,
				atom_indices1, atom_indices2, len1, len2);
		ftime(&tp2);
#ifdef VERBOSE_PRINTS
		printf("*****tmalign_time1 = %ld.%d\n",tp1.time,tp1.millitm);
		printf("*****tmalign_time2 = %ld.%d\n",tp2.time,tp2.millitm);
//		printf("\nresults: %f\t%f\t%f\t%d\t%d\n",master_in_data.tm1,
//				master_in_data.tm2,master_in_data.rmsd,master_in_data.start_time,master_in_data.end_time);
		printf("finished TMalign\n");
#endif
		master_in_data.job_index = client_in_data.job_index;
		master_in_data.processing_time = diff_timeb(tp1, tp2);
		master_in_data.algo_type = TMALIGN_ALGO_TYPE;
		delete[] backbone_1;
	}

	ftime(&t);
	master_in_data.seconds_at_result_send = t.time;
	return RCKSKEL_FALSE;
}
#endif
////////////////////////////////////////////////////////////////////
void CE::scratch_align_ent(char *db_tmp_path, char *pdb_dir_name) {
	// master node specific work

	int file_count = get_file_count(pdb_dir_name);
	char * filenames[file_count];
	fill_dataset_details(pdb_dir_name, filenames);
	printf("filecount: %d\n", file_count);

	protein_data proteins_data[file_count];
	parse_dataset_files_make_entries(pdb_dir_name, filenames, file_count,
			db_tmp_path, proteins_data);

	// mapable processing
	for (int i = 0; i < file_count; i++) {
		for (int j = 0; j < file_count; j++) {
			printf("%s vs %s\n", proteins_data[i].name, proteins_data[j].name);
			printf(
					"\nStructure Alignment Calculator, version 1.02, last modified: Jun 15, 2001.\n\n");
			// slave node specific work
			align_pair(&proteins_data[i], &proteins_data[j], 8, 3.0, 4.0);
		}
	}

	// free stuff
}

////////////////////////////////////////////////////////////////////
void CE::align_pair(protein_data *protein1, protein_data *protein2,
		int winSize_, double rmsdThr, double rmsdThrJoin) {
	double d_[20];
	int isPrint = 2;
	int lcmp;

	int *align_se1 = new int[protein1->nSe + protein2->nSe];
	int *align_se2 = new int[protein1->nSe + protein2->nSe];

	char *seq1 = protein1->seq, *seq2 = protein2->seq;

	//NOTE: this is the real meat of the work.
	int aln_len, gaps;
	float rmsd;
	ce_1(protein1->name, protein2->name, protein1->ca, protein2->ca,
			protein1->nSe, protein2->nSe, align_se1, align_se2, lcmp, 8, 3.0,
			4.0, d_, isPrint, &aln_len, &gaps, &rmsd);

	// print scores and alignments
	isPrint = 1;
	if (lcmp > 0) {
		int lsim = 0, lali = 0;
		for (int l = 0; l < lcmp; l++)
			if (align_se1[l] != -1 && align_se2[l] != -1) {
				if (seq1[align_se1[l]] == seq2[align_se2[l]])
					lsim++;
				lali++;
			}
		printf("Sequence identities = %.1f%%", lsim * 100.0 / lali);

		int lstep = 70;

		for (int l = 0; l < lcmp; l += lstep) {
			printf("\n");
			for (int ie = 0; ie < 2; ie++) {

				char *seq = (ie == 0 ? seq1 : seq2);

				int *align_se = (ie == 0 ? align_se1 : align_se2);

				printf("\n%8.8s ", (ie == 0 ? "Chain 1:" : "Chain 2:"));
				int ip = -1;
				for (int l_ = l; l_ < l + lstep && l_ < lcmp; l_++)
					if (align_se[l_] != -1) {
						ip = align_se[l_] + 1;
						break;
					}
				if (ip != -1)
					printf("%4d ", ip);
				else
					printf("     ");

				for (int l_ = l; l_ < l + lstep && l_ < lcmp; l_++)
					printf("%c",
							align_se[l_] == -1 ? '-' : (seq[align_se[l_]]));

			}
		}
		printf("\n");

		printf("\n     X2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",
				d_[0], d_[1], d_[2], d_[9]);
		printf("     Y2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",
				d_[3], d_[4], d_[5], d_[10]);
		printf("     Z2 = (%9.6f)*X1 + (%9.6f)*Y1 + (%9.6f)*Z1 + (%12.6f)\n",
				d_[6], d_[7], d_[8], d_[11]);

	}

	// free stuff
	if (protein1->nSe + protein2->nSe > 0) {
		delete[] align_se1;
		delete[] align_se2;
	}
}
////////////////////////////////////////////////////////////////////
void CE::fill_dataset_details(char *dir_name, char **filenames) {
	DIR *mydir = opendir(dir_name);
	struct dirent *entry = NULL;
	int ccount = 0;

	// read file name and store
	while ((entry = readdir(mydir))) {
		if (strcmp(entry->d_name, ".") == 0
				|| strcmp(entry->d_name, "..") == 0) {
			continue;
		}
		char *name;
		setText(&name, entry->d_name);
		filenames[ccount++] = name;
	}
	closedir(mydir);
}
////////////////////////////////////////////////////////////////////
void CE::parse_dataset_files_make_entries(char *data_dir, char** filenames,
		int file_count, char *db_tmp_path, protein_data *proteins_data) {
	char cmd[] = "scratch";
	char *parms[6];
	int parm_count;

	// initialize scratch
	parms[1] = cmd;
	char slash_cmd[] = "/";
	for (int i = 0; i < file_count; i++) {
		char *mkdir_cmd1;
		char rm_cmd[] = "rm -rf ";
		setText(&mkdir_cmd1, rm_cmd);
		addText(&mkdir_cmd1, db_tmp_path);
		addText(&mkdir_cmd1, slash_cmd);
		addText(&mkdir_cmd1, filenames[i]);
		system(mkdir_cmd1);

		char *mkdir_cmd;
		char mk_cmd[] = "mkdir -p ";
		setText(&mkdir_cmd, mk_cmd);
		addText(&mkdir_cmd, db_tmp_path);
		addText(&mkdir_cmd, slash_cmd);
		addText(&mkdir_cmd, filenames[i]);
		system(mkdir_cmd);
		// db path
		char *db_path;
		setText(&db_path, db_tmp_path);
		addText(&db_path, slash_cmd);
		addText(&db_path, filenames[i]);
		parms[2] = db_path;
		// get qualified path to file
		char *file_path;
		setText(&file_path, data_dir);
		addText(&file_path, slash_cmd);
		addText(&file_path, filenames[i]);
		//
		parms[3] = file_path;
		parm_count = 4;
		alt_entry(parm_count, parms);
	}

	for (int i = 0; i < file_count; i++) {
		//
		DB db;
		char *db_path;
		setText(&db_path, db_tmp_path);
		addText(&db_path, slash_cmd);
		addText(&db_path, filenames[i]);
		db.setPath(db_path);

		char fname1[] = "name.enp", fname2[] = "n_se.enp", fname3[] = "c_a.enp",
				fname4[] = "seq.enp", fname5[] = "code3.mon", fname6[] =
						"se.enc", fname7[] = "i_enc.enp", fname8[] = "id.com",
				fname9[] = "compnd.com", fname10[] = "i_com.enp", fname11[] =
						"i_enp.com";
		Property name_enp(fname1), nse_enp(fname2), ca_enp(fname3), seq_enp(
				fname4), code3_mon(fname5), se_enc(fname6), i_enc_enp(fname7),
				id_com(fname8), comp_enp(fname9), i_com_enp(fname10), i_enp_com(
						fname11);

		protein_data protein;
		char ent[] = "USR1:";
		char *com;
		setText(&com, ent);
		com[4] = '\0';
		protein.iCom = id_com.find(com);
		protein.iEnp1 = *i_enp_com.item4(protein.iCom);
		protein.iEnp2 = *i_enp_com.item4(protein.iCom + 1);
		if (protein.iEnp1 == -1 || protein.iEnp2 == -1) {
			printf("Chain %s not found\n", filenames[i]);
			continue;
		}
		//
		char *name;
		setText(&name, filenames[i]);
		char col_cmd[] = ":";
		addText(&name, col_cmd);
		addText(&name, name_enp.item1(protein.iEnp1) + 5);
		//
		protein.name = name;
		protein.ent = ent;
		protein.nSe = *nse_enp.item2(protein.iEnp1);
		//protein.se = se_enc.item2(*i_enc_enp.item4(protein.iEnp1), 1);
		protein.seq = seq_enp.item1(protein.iEnp1, 1);
		flt4 *raw_struct = ca_enp.itemf(protein.iEnp1);
		protein.raw_struct = raw_struct;
		protein.ca = arrayToXYZ(raw_struct, protein.nSe);
		proteins_data[i] = protein;
	}
}
////////////////////////////////////////////////////////////////////
int get_file_count(char *dir_name) {
	///TODO: can we get the count in a better way?
	// get file count
	int file_count = 0;
	DIR *mydir = opendir(dir_name);
	struct dirent *entry = NULL;
	while ((entry = readdir(mydir)))
		file_count += 1;
	closedir(mydir);
	file_count -= 2;
	return file_count;
}
////////////////////////////////////////////////////////////////////
void read_params(char *db_tmp_path, char *pdb_dir_name, char **argv) {
	// get temp folder path & data dir
	printf("argv1: %s\n", argv[1]);
	setText(&db_tmp_path, argv[1]);
	setText(&pdb_dir_name, argv[2]);
}
///////////////////////////////////////////////////////////////////////////
void ce_load_data(CE *ce, char *db_tmp_path, char *pdb_dir_name, int file_count,
		char **filenames, protein_data *proteins_data) {
	// load data
	ce->fill_dataset_details(pdb_dir_name, filenames);
	ce->parse_dataset_files_make_entries(pdb_dir_name, filenames, file_count,
			db_tmp_path, proteins_data);
}
///////////////////////////////////////////////////////////////////////////
void tmalign_load_data(int line_count, char **filenames, char *dir_name) {
	int i;
	// every 5000 belongs to one protein.
	seq_amino_acids = (char *) malloc(
			3 * MAX_ATOMS * line_count * sizeof(char));
	amino_acid_chain_seq = (char *) malloc(
			1 * 5001 * line_count * sizeof(char));
	atom_index = (int *) malloc(MAX_ATOMS * line_count * sizeof(int));
	// every 15,000 belongs to one protein. order assumed.
	seq_coords = (float *) malloc(3 * MAX_ATOMS * line_count * sizeof(float));
	seq_length = (int *) malloc(MAX_ATOMS * line_count * sizeof(int));

	// read files
	char str[256];
	for (i = 0; i < line_count; i++) {
		strcpy(str, dir_name);
		strcat(str, filenames[i]);
		seq_length[i] = tmalign_read_pdb_structure(str,
				&seq_coords[3 * MAX_ATOMS * i],
				&seq_amino_acids[3 * MAX_ATOMS * i], &zero,
				&atom_index[MAX_ATOMS * i], &amino_acid_chain_seq[5001 * i]);
		//printf("seq: %c\n",amino_acid_chain_seq[5001*i+1]);
	}
}

int tmalign_read_pdb_structure(char *fname, float *prot_backbone, char *ss,
		int *start, int *lenarray, char *seq) {
	char aa[22][4] = { "BCK", "GLY", "ALA", "SER", "CYS", "VAL", "THR", "ILE",
			"PRO", "MET", "ASP", "ASN", "LEU", "LYS", "GLU", "GLN", "ARG",
			"HIS", "PHE", "TYR", "TRP", "CYX" };
	char slc[] = "XGASCVTIPMDNLKEQRHFYWC";

	//printf("Reading File: %s|\n", fname);
	FILE *file = fopen(fname, "r");
	char line[500];
	float x, y, z;
	int t;
	char chain;
	char aaname[3];
	int index = 1, j = 0;

	if (file == 0) {
		perror(fname);
		return -2;
	}
	char buffer1[3], buffer2[4], buffer3;
	while (fgets(line, sizeof line, file) != NULL) {
		strncpy(buffer1, line, 3);
		strncpy(buffer2, line + 12, 4);
		strncpy(&buffer3, line + 16, 1);
		//printf("buffer1: %s, buffer2: %s\n",buffer1,buffer2);
		// if line starts with TER we're done reading
		if (strncmp(buffer1, "TER", 3) == 0 && index - 1 > 0) {
			//printf("Done reading: found TER.\n");
			break;
		}
		if (strncmp(buffer1, "ATO", 3) != 0) {
			//printf("Skipping line: no ATOM line\n");
			continue;
		}
		if (strncmp(buffer2, "CA  ", 4) != 0 && strncmp(buffer2, " CA ", 4) != 0
				&& strncmp(buffer2, "  CA", 4) != 0) {
			//printf("line + 12: %s\n", line + 12);
			//printf("line: %s\n", line);
			//printf("Skipping line: no CA on line in required location\n");
			continue;
		}
		//printf("chain:%c|\n",buffer3);
		sscanf(&line[17], "%s %c%d", aaname, &chain, &t);
		sscanf(&line[30], "%f%f%f", &x, &y, &z);
		if (buffer3 != 'A' && buffer3 != ' ') {
			//printf("Wrong chain\n");
			continue;
		}
		lenarray[index - 1] = t;
		prot_backbone[index * 3 - 3 + *start] = x;
		prot_backbone[index * 3 - 2 + *start] = y;
		prot_backbone[index * 3 - 1 + *start] = z;
		seq[index] = slc[0];
		for (j = -1; j <= 20; ++j) {
			if (strncmp(aa[j + 1], aaname, 3) == 0) {
				seq[index] = slc[j + 1];
				break;
			}
		}
		strcpy(ss + (index - 1) * 3, aaname);
		//printf("aanam: %s, initial4_1.mm1: %d, %f, %f, %f, %c, %s\n",
		//       aaname,lenarray[index-1],x,y,z,
		//      seq[index], ss + (index-1)*3
		//);
		if (index++ == 5000) {
			break;
		}
	}
	//printf("index: %d\n",index);
	fclose(file);

	return --index;
}
////////
int isQuery(char *name, std::set<std::string> lines) {
	string n = name;
	for(auto line : lines)
		if(n.find(line) == 0)	
			return 1;
	return 0;
}

#ifdef PC_TEST
// Target PC
void do_target_pc(int argc, char **argv) {
	struct timeb tp1, tp2, otp1, otp2;
#ifdef SCC_1_CORE_TEST
        rckskel_env_t env1;
        rckskel_env_init(&env1, &argc, &argv);
        int ue_count1 = env1.ue_count, ue_id1 = env1.ue_id;
        printf("starting ue %d\n", ue_id1);

        if (ue_id1 == MASTER_ID) {
#endif
	// read params
	char *db_tmp_path, *pdb_dir_name, *pdb_cm_name;
	printf("reading params\n");
	//read_params(db_tmp_path, pdb_dir_name, argv);
	setText(&db_tmp_path, argv[1]);
	setText(&pdb_dir_name, argv[2]);
	setText(&pdb_cm_name, argv[3]);
	printf("db_tmp_path: %s\n", db_tmp_path);

#ifdef MANY_MANY
	std::ifstream inf(argv[4]);
	std::set<std::string> lines;
	std::string line;
	for (unsigned int i=1; std::getline(inf,line); ++i)
        	lines.insert(line);
	cout << "loaded query names" << endl;
#endif

	// get file count
	int file_count = get_file_count(pdb_dir_name);
	char *filenames[file_count];
	protein_data proteins_data[file_count];
	int cm_file_count = get_file_count(pdb_cm_name);
	char *cm_filenames[cm_file_count];
	fill_dataset_details(pdb_cm_name, cm_filenames);

	CE ce;
	printf("loading CE data\n");
	// load CE data and fill filename as side-effect
	ftime(&tp1);
	ce->fill_dataset_details(pdb_dir_name, filenames);
        ce->parse_dataset_files_make_entries(pdb_dir_name, filenames, file_count,
                        db_tmp_path, proteins_data);

	ftime(&tp2);
	cout << "CE data load time: " << diff_timeb(tp1, tp2) << endl;

	printf("loading TMalign data\n");
	ftime(&tp1);
	tmalign_load_data(file_count, filenames, pdb_dir_name);
	ftime(&tp2);
	cout << "TMalign data load time: " << diff_timeb(tp1, tp2) << endl;
	backbone_ backbone_1;

	printf("loading USM data\n");
	ftime(&tp1);
	USM usm;
	usm.load_data(cm_file_count, pdb_cm_name, cm_filenames);
	ftime(&tp2);
	cout << "USM data load time: " << diff_timeb(tp1, tp2) << endl;

	ftime(&otp1);
#ifndef MANY_MANY
	usm.calculate_pairwise_distances();
#endif
#ifdef MANY_MANY
	usm.calculate_pairwise_distances(lines);
#endif
	ftime(&otp2);
	printf("USM Total time (msec): %ld\n", diff_timeb(otp1, otp2));

	// do TMalign
	ftime(&otp1);
	for (int i = 0; i < file_count; i++) {
#ifdef MANY_MANY
	if(isQuery(proteins_data[i].name, lines) == 0)
		continue;
#endif
                for (int j = 0; j < file_count; j++) {
#ifdef MANY_MANY
			if(isQuery(proteins_data[j].name, lines) == 1)
				continue;
#endif
			// do TMalign
			int cntr;
			char seq1_amino_acids[3 * MAX_ATOMS],
					seq2_amino_acids[3 * MAX_ATOMS];
			char amino_acid_chain_seq1[5001], amino_acid_chain_seq2[5001];
			int atom_indices1[MAX_ATOMS], atom_indices2[MAX_ATOMS];
			backbone_ backbone_1;

			struct timeb tp1, tp2;
			ftime(&tp1);
			// now collect the batch of data
			for (cntr = 0; cntr < 5001; cntr++) {
				amino_acid_chain_seq1[cntr] = amino_acid_chain_seq[(5001 * i)
						+ cntr];
				amino_acid_chain_seq2[cntr] = amino_acid_chain_seq[(5001 * j)
						+ cntr];
			}
			// atom_index
			for (cntr = 0; cntr < MAX_ATOMS; cntr++) {
				atom_indices1[cntr] = atom_index[(MAX_ATOMS * i) + cntr];
				atom_indices2[cntr] = atom_index[(MAX_ATOMS * j) + cntr];
			}
			for (cntr = 0; cntr < 3 * MAX_ATOMS; cntr++) {
				seq1_amino_acids[cntr] = seq_amino_acids[3 * MAX_ATOMS * i
						+ cntr];
				seq2_amino_acids[cntr] = seq_amino_acids[3 * MAX_ATOMS * j
						+ cntr];
			}
			// seq_coords
			for (cntr = 0; cntr < 3 * MAX_ATOMS; cntr++) {
				backbone_1.xa[cntr] = seq_coords[i * 3 * MAX_ATOMS + cntr];
				backbone_1.xa[3 * MAX_ATOMS + cntr] = seq_coords[j * 3
						* MAX_ATOMS + cntr];
			}
			// lengths
			//length_1.nseq1 = seq_length[i];
			//length_1.nseq2 = seq_length[j];

			ftime(&tp1);
			ClientToMasterTransferBlock ldata;
			main_process(seq1_amino_acids, seq2_amino_acids,
					amino_acid_chain_seq1, amino_acid_chain_seq2, &ldata,
					&backbone_1, atom_indices1, atom_indices2, seq_length[i],
					seq_length[j]);
			ftime(&tp2);

			printf("Algo: %d, Protein1: %s, Protein2: %s, len1: %d, len2: %d, "
					"Rmsd: %.2fA, TM1: %.2f, TM2: %.2f, "
					"id: %d, sec_job_start: %ld, sec_result_send: %ld, "
					"processing_time_msec: %ld, data_collect_time_msec: %ld, "
					"idle_time_msec: %ld\n", TMALIGN_ALGO_TYPE,
					proteins_data[i].name, proteins_data[j].name, seq_length[i],
					seq_length[j], ldata.rmsd, ldata.tm1, ldata.tm2, 0,
					tp1.time, tp2.time, diff_timeb(tp1, tp2), 0, 0);
		}
	}
	ftime(&otp2);
	printf("TMalign Total time (msec): %ld\n", diff_timeb(otp1, otp2));

	// do CE
	ftime(&otp1);
	for (int i = 0; i < file_count; i++) {
#ifdef MANY_MANY
        if(isQuery(proteins_data[i].name, lines) == 0)
                continue;
#endif
		for (int j = 0; j < file_count; j++) {
#ifdef MANY_MANY
                        if(isQuery(proteins_data[j].name, lines) == 1)
                                continue;
#endif
			ftime(&tp1);

			double d_[20];
			int nse1 = proteins_data[i].nSe, nse2 = proteins_data[j].nSe;
			int isPrint = 0, lcmp, *align_se1 = new int[nse1 + nse2];
			int *align_se2 = new int[nse1 + nse2];

			int aln_len, gaps;
			float rmsd;
			double z = ce_1(proteins_data[i].name, proteins_data[j].name,
					proteins_data[i].ca, proteins_data[j].ca,
					proteins_data[i].nSe, proteins_data[j].nSe, align_se1,
					align_se2, lcmp, 8, 3.0, 4.0, d_, isPrint, &aln_len, &gaps,
					&rmsd);
			// free stuff
			delete[] align_se1;
			delete[] align_se2;
			ftime(&tp2);

			printf("Algo: %d, Protein1: %s, Protein2: %s, len1: %d, len2: %d, "
					"Rmsd: %.2fA, TM1: %.2f, TM2: %.2f, "
					"id: %d, sec_job_start: %ld, sec_result_send: %ld, "
					"processing_time_msec: %ld, data_collect_time_msec: %ld, "
					"idle_time_msec: %ld\n", CE_ALGO_TYPE,
					proteins_data[i].name, proteins_data[j].name, nse1, nse2,
					rmsd, -1.0, -1.0, 0, tp1.time, tp2.time,
					diff_timeb(tp1, tp2), 0, 0);
		}
	}
	ftime(&otp2);
	printf("CE Total time (msec): %ld\n", diff_timeb(otp1, otp2));
#ifdef SCC_1_CORE_TEST
	}
#endif
}
#endif

#ifndef PC_TEST
// Target SCC N-core
void do_target_scc(int argc, char **argv) {
	struct timeb tp1, tp2, otp1, otp2;
	rckskel_env_t env;
	rckskel_env_init(&env, &argc, &argv);
	int ue_count = env.ue_count, ue_id = env.ue_id;
	printf("starting ue %d\n", ue_id);

	if (ue_id == MASTER_ID) {
		// read params
		char *db_tmp_path, *pdb_dir_name, *pdb_cm_name;
		setText(&db_tmp_path, argv[1]);
		setText(&pdb_dir_name, argv[2]);
		setText(&pdb_cm_name, argv[3]);

		printf("uecount: %d\n", ue_count);
		int ue_ids[ue_count];
		for (int i = 0; i < ue_count; i++) {
			ue_ids[i] = i;
		}
#ifdef MANY_MANY
		std::ifstream inf(argv[4]);
		std::set<std::string> lines;
		std::string line;
		for (unsigned int i=1; std::getline(inf,line); ++i)
        		lines.insert(line);
		cout << "loaded query names" << endl;
#endif

		// Process USM
//		int cm_file_count = get_file_count(pdb_cm_name);
//		char *cm_filenames[cm_file_count];
//		fill_dataset_details(pdb_cm_name, cm_filenames);

//        	ftime(&tp1);
//	        USM usm;
//        	usm.load_data(cm_file_count, pdb_cm_name, cm_filenames);
//	        ftime(&tp2);
//        	cout << "USM data load time: " << diff_timeb(tp1, tp2) << endl;

//	        ftime(&otp1);
#ifndef MANY_MANY
//       	usm.calculate_pairwise_distances();
#endif
#ifdef MANY_MANY
//		usm.calculate_pairwise_distances(lines);
#endif
// 	        ftime(&otp2);
//        	printf("USM Total time (msec): %ld\n", diff_timeb(otp1, otp2));
		//delete &usm;
		//delete[] cm_filenames;

		// process CE and TMalign
#ifndef MANY_MANY
		int file_count = get_file_count(pdb_dir_name);
		printf("file_count: %d\n", file_count);
		char *filenames[file_count];
		fill_dataset_details(pdb_dir_name, filenames);
#endif
#ifdef MANY_MANY
		cout << "many many processing" << endl;
                int t_file_count = get_file_count(pdb_dir_name);
                printf("file_count: %d\n", t_file_count);
                char *t_filenames[t_file_count];
                fill_dataset_details(pdb_dir_name, t_filenames);
		int max_batch_size = 30;

		char *q_filenames[lines.size()], *r_filenames[t_file_count - lines.size()];

		// separate names into query and data
		cout << "separate names into query and data" << endl;
		int q_length = 0, r_length = 0;
		for(int mm_ctr = 0; mm_ctr < t_file_count; mm_ctr++) {
			 if(isQuery(t_filenames[mm_ctr], lines) != 0)
				q_filenames[q_length++] = t_filenames[mm_ctr];
			else
				r_filenames[r_length++] = t_filenames[mm_ctr];
		}

		int batch_size = max_batch_size - q_length;
		cout << "data batch size: " << batch_size << endl;

		//
		for(int mm_b_ctr = 0; mm_b_ctr < r_length; mm_b_ctr += batch_size) {

		int file_count = q_length + ((mm_b_ctr + batch_size > r_length) ? (r_length - mm_b_ctr) : batch_size);

		char *filenames[file_count];
		// create filenames for this batch
		int mm_index = 0;
		for(mm_index = 0; mm_index < q_length; mm_index++)
			filenames[mm_index] = q_filenames[mm_index];
		for(mm_index = 0; mm_index < r_length; mm_index++)
			filenames[q_length + mm_index] = r_filenames[mm_b_ctr + mm_index];
		
#endif
		CE ce;

		ftime(&tp1);
		// load CE data and fill filename as side-effect
		cout << "loading ce data" << endl;

		protein_data proteins_data[file_count];
        	ce.parse_dataset_files_make_entries(pdb_dir_name, filenames, file_count,
                        db_tmp_path, proteins_data);

		ftime(&tp2);
		cout << "CE data load time: " << diff_timeb(tp1, tp2) << endl;

		cout << "loading tmalign data" << endl;
		ftime(&tp1);
		tmalign_load_data(file_count, filenames, pdb_dir_name);
		ftime(&tp2);
		cout << "TM data load time: " << diff_timeb(tp1, tp2) << endl;

		// make task
		printf("making jobs\n");
#ifndef MANY_MANY
		job_pool = (job_descriptor *) malloc(
				2 * file_count * file_count * sizeof(job_descriptor));
		job_indexes = (int *) malloc(2 * file_count * file_count * sizeof(int));
#endif
#ifdef MANY_MANY
		job_pool = (job_descriptor *) malloc(
				2 * q_length * r_length * sizeof(job_descriptor));
		job_indexes = (int *) malloc(2 * q_length * r_length * sizeof(int));
#endif
		int jcount = 0;
		// place jobs in job pool
		for (int i = file_count - 1; i >= 0; i--) {
#ifdef MANY_MANY
			if(isQuery(proteins_data[i].name, lines) == 0)
				continue;
#endif
			for (int j = 0; j < file_count; j++) {
#ifdef MANY_MANY
				if(isQuery(proteins_data[j].name, lines) == 1)
					continue;
#endif
				// ce task
				job_pool[jcount].algo_type = CE_ALGO_TYPE;
				job_pool[jcount].protein1 = &proteins_data[i];
				job_pool[jcount].protein2 = &proteins_data[j];
				job_pool[jcount].protein1_idx = i;
				job_pool[jcount].protein2_idx = j;
				job_indexes[jcount] = jcount++;
				// tmalign task
				job_pool[jcount].algo_type = TMALIGN_ALGO_TYPE;
				job_pool[jcount].protein1 = &proteins_data[i];
				job_pool[jcount].protein2 = &proteins_data[j];
				job_pool[jcount].protein1_idx = i;
				job_pool[jcount].protein2_idx = j;
				job_indexes[jcount] = jcount++;
			}
		}
		cout << "made job lists" << endl;
#ifdef RANDOM_SPLIT
		// create random partitions for CE & TMalign jobs
		std::vector<SortData> ce_jobs = get_partitions(ue_count - 1, CE_ALGO_TYPE, jcount, 0);
		std::vector<SortData> tm_jobs = get_partitions(ue_count - 1, TMALIGN_ALGO_TYPE, jcount, 0);
		//
		ce_jobs.insert(ce_jobs.end(), tm_jobs.begin(), tm_jobs.end());
		std::random_shuffle(ce_jobs.begin(), ce_jobs.end());
		// make job list
		int index = 0;
		for(auto datum : ce_jobs)
			job_indexes[index++] = datum.index;
#endif
#ifdef GREEDY_SPLIT
		// create greedy partitions for CE & TMalign jobs
		std::vector<SortData> sorted_ce_jobs = get_partitions(ue_count - 1, CE_ALGO_TYPE, jcount, 1);
		std::vector<SortData> sorted_tm_jobs = get_partitions(ue_count - 1, TMALIGN_ALGO_TYPE, jcount, 1);
		// make job list
		int index = 0;
		for(auto datum : sorted_ce_jobs)
			job_indexes[index++] = datum.index;
		for(auto datum : sorted_tm_jobs)
			job_indexes[index++] = datum.index;	
#endif
		cout << "sorted jobs" << endl;
		task_t task;
		create_task(&task, FARM_rckskel, '0');
		init_task(&task);
		task.worker = &master_send_job_data;
		task.collector = &master_receive_result;
		task.job_idxs = job_indexes;
		set_task_ues(&task, ue_count - 1, &ue_ids[1]);
		set_unfinished_count(&task, jcount);

		// run task
		check_ues_up(&ue_count);
		ftime(&tp1);
		cout << "start process" << endl;
		rckskel_process_task_tree(&task, ue_count, ue_ids);
		cout << "end process" << endl;
		ftime(&tp2);
		printf("Time at process start = %ld.%d\n", tp1.time, tp1.millitm);
		printf("Time at process end = %ld.%d\n", tp2.time, tp2.millitm);
		cout << "Iteration run time: " << diff_timeb(tp1, tp2) << endl;

#ifdef MANY_MANY
		delete[] job_pool;
		delete[] job_indexes;
		delete[] proteins_data;
}
#endif
		// terminate clients
		for (int i = 0; i < ue_count; i++) {
			terminate(i);
		}

	} else {
		timeb t1;
		// send ready to master
		READY_T ready_t;
		ready_t.ready = 1;
		ftime(&t1);
		RCCE_send((char *) &ready_t, sizeof(READY_T), MASTER_ID);

		// process till die
		while (RCKSKEL_TRUE) {
			// get and process job
			if (client_receive_job(t1) == RCKSKEL_TRUE)
				break;
			ftime(&t1);
			// send result of processing to master
			RCCE_send((char *) &master_in_data,
					sizeof(ClientToMasterTransferBlock), MASTER_ID);
		}
	}

	printf("ue %d exiting.\n", ue_id);

}
#endif
////////
