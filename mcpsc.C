#include "comm_block.h"
#include "ce.h"
#include "usm.hpp"
#include "tmalign_extern.h"
#include "pom/ipdb.h"
#include <sys/timeb.h>
#include <stdio.h>
#include <omp.h>
#define _GNU_SOURCE
#include <stdlib.h>

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

int *job_indexes;
job_descriptor *job_pool;

// local function definitions
void master_send_job_data(int job_idx, int ue_id);
int client_receive_job(timeb t1);
int master_receive_result(int id);
int check_ready(int id);
void terminate(int id);
void tmalign_load_data(int line_count, char **filenames, char *dir_name);
int tmalign_read_pdb_structure(char *fname, float *prot_backbone, char *ss,
		int *start, int *lenarray, char *seq);
int get_file_count(char *dir_name);
void ce_load_data(CE *ce, char *db_tmp_path, char *pdb_dir_name, int file_count,
		char **filenames, protein_data *proteins_data);
void read_params(char *db_tmp_path, char *pdb_dir_name, char **argv);
int get_file_count2(char *dom_files);
void fill_filenames(int file_count1, char **filenames1, char *dom_files);
int get_job_count(char *pair_indices);
void load_job_pairs(int job_pairs[][2], char *pair_indices);

MasterToClientTransferBlock client_in_data;
ClientToMasterTransferBlock master_in_data;

int diff_timeb(timeb t1, timeb t2) {
	return (1000 * t2.time + t2.millitm) - (1000 * t1.time + t1.millitm);
}

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


void dispatch(int ue_count, int job_index) {
	while(1) {
		for(int i = 1; i < ue_count; i++) {
			// if node is accepting send data to it
		}
	}
}

void do_usm(char *pdb_cm_name, char *cm_files, int jcount, int jobs[][2]) {
	struct timeb tp1, tp2, otp1, otp2;
	// do usm and wrap up
	int cm_file_count = get_file_count2(cm_files);
	char *cm_filenames[cm_file_count];
	fill_filenames(cm_file_count, cm_filenames, cm_files);
    ftime(&tp1);
	USM usm;
    usm.load_data(cm_file_count, pdb_cm_name, cm_filenames);
	ftime(&tp2);
    cout << "USM data load time: " << diff_timeb(tp1, tp2) << endl;

	ftime(&otp1);
	// usm.calculate_pairwise_distances();
	usm.calculate_pairwise_distance(jcount, jobs, cm_filenames);
	ftime(&otp2);
	cout << "USM Total time (msec): " << diff_timeb(otp1, otp2) << endl;
}

int main(int argc, char **argv) {
	struct timeb tp1, tp2;
	struct timeb otp1, otp2;
	// read params
	char *db_tmp_path, *pdb_dir_name, *pdb_cm_name, *dom_files, 
		*pair_indices, *cm_files;
	printf("reading params\n");
	//read_params(db_tmp_path, pdb_dir_name, argv);
	setText(&db_tmp_path, argv[1]);
	setText(&pdb_dir_name, argv[2]);
	setText(&pdb_cm_name, argv[3]);
	setText(&dom_files, argv[4]);
	setText(&pair_indices, argv[5]);
	setText(&cm_files, argv[6]);
	printf("db_tmp_path: %s\n", db_tmp_path);

	// get file count, filenames and load pairs for ce, tmalign
	int file_count = get_file_count2(dom_files);
	char *filenames[file_count];
	fill_filenames(file_count, filenames, dom_files);
	int jcount = get_job_count(pair_indices);
	int jobs[jcount][2];
	load_job_pairs(jobs, pair_indices);

	// get file count
	protein_data proteins_data[file_count];

#ifndef ONLY_CE
#ifndef ONLY_TMALIGN
	do_usm(pdb_cm_name, cm_files, jcount, jobs);
#endif
#endif


	int nthreads, tid, counter, chunk = 1;
#ifndef ONLY_CE
#ifndef ONLY_USM
	printf("loading TMalign data\n");
	ftime(&tp1);
	tmalign_load_data(file_count, filenames, pdb_dir_name);
	ftime(&tp2);
	cout << "TMalign data load time: " << diff_timeb(tp1, tp2) << endl;
	//backbone_ backbone_1;
	// do TMalign
	ftime(&otp1);
	printf("jcount: %d\n",jcount);		

#pragma omp parallel shared(jobs,amino_acid_chain_seq,atom_index,seq_amino_acids,seq_coords,jcount,nthreads,chunk) private(counter, tid) 
{
	tid = omp_get_thread_num();
	if (tid == 0) {
		nthreads = omp_get_num_threads();
		printf("Number of threads = %d\n", nthreads);
	}
	printf("Thread %d starting...\n",tid);

	#pragma omp for schedule(dynamic,chunk)
	for (counter = 0; counter < jcount; counter++) {
		int i = jobs[counter][0], j = jobs[counter][1];
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
					filenames[i], filenames[j], seq_length[i],
					seq_length[j], ldata.rmsd, ldata.tm1, ldata.tm2, tid,
					tp1.time, tp2.time, diff_timeb(tp1, tp2), 0, 0);
	}
} // end of parallel block
	ftime(&otp2);
	printf("TMalign Total time (msec): %ld\n", diff_timeb(otp1, otp2));
#endif
#endif

#ifndef ONLY_TMALIGN
#ifndef ONLY_USM
	CE ce;
	printf("loading CE data\n");
	// load CE data and fill filename as side-effect
	ftime(&tp1);
	ce_load_data(&ce, db_tmp_path, pdb_dir_name, file_count, filenames,
			proteins_data);
	ftime(&tp2);
	cout << "CE data load time: " << diff_timeb(tp1, tp2) << endl;

	// do CE
	ftime(&otp1);
	printf("jcount: %d\n",jcount);		

#pragma omp parallel shared(jobs,proteins_data,jcount,nthreads,chunk) private(counter, tid) 
{
	tid = omp_get_thread_num();
	if (tid == 0) {
		nthreads = omp_get_num_threads();
		printf("Number of threads = %d\n", nthreads);
	}
	printf("Thread %d starting...\n",tid);

	#pragma omp for schedule(dynamic,chunk)
	for (counter = 0; counter < jcount; counter++) {
		int i = jobs[counter][0], j = jobs[counter][1];
			// do TMalign
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
					rmsd, -1.0, -1.0, tid, tp1.time, tp2.time,
					diff_timeb(tp1, tp2), 0, 0);
	}
} // end of parallel
	ftime(&otp2);
	printf("CE Total time (msec): %ld\nq", diff_timeb(otp1, otp2));
#endif
#endif
	return 0;
}


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
		char mk_cmd[] = "mkdir ";
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
	//ce->fill_dataset_details(pdb_dir_name, filenames);
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
		//printf("%s\n",filenames[i]);
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
int get_file_count2(char *dom_files) {
	FILE *file = fopen(dom_files, "r");
	int i;
	fscanf (file, "%d", &i);
	fclose(file);
	return i;
}
void fill_filenames(int file_count1, char **filenames1, char *dom_files) {
	FILE *file = fopen(dom_files, "r");
	size_t len = 0;
	ssize_t read;
	char * line = NULL;
	int lcount = 0;
	while ((read = getline(&line, &len, file)) != -1) {
            line[read - 1] = '\0';
			--read;
			if (lcount++ == 0) continue;
			char *domname;
			setText(&domname, line);
			filenames1[lcount - 2] = domname;
	}
	fclose(file);
}
int get_job_count(char *pair_indices) {
	FILE *file = fopen(pair_indices, "r");
	int i;
	fscanf (file, "%d", &i);
	fclose(file);
	return i;
}
void load_job_pairs(int job_pairs[][2], char *pair_indices) {
	FILE *file = fopen(pair_indices, "r");
	int i, count, n1, n2;
	char line[80];
	fgets(line, 80, file);
	sscanf(line, "%d", &count);
	for(i = 0; i < count; i++) {
		fgets(line, 80, file);
		sscanf(line, "%d %d", &n1, &n2);
		job_pairs[i][0] = n1;
		job_pairs[i][1] = n2;
	}
	fclose(file);
}
////////
