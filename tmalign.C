/**************************************************************************
*     This program is to identify the best alignment of two protein 
*     structures that gives the highest TM-score. Input structures must 
*     be in the PDB format. By default, TM-score is normalized by the 
*     second protein. Users can obtain a brief instruction by simply 
*     running the program without arguments. For comments/suggestions,
*     please contact email: zhng@umich.edu.
*     
*     Reference to cite:
*     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2303-9
*     
*     Permission to use, copy, modify, and distribute this program for 
*     any purpose, with or without fee, is hereby granted, provided that
*     the notices on the head, the reference information, and this
*     copyright notice appear in all copies or substantial portions of 
*     the Software. It is provided "as is" without express or implied 
*     warranty.
************************ updating history ********************************
*     2005/06/01: A bug of two-point superposition was fixed.
*     2005/10/19: the program was reformed so that the alignment
*                 results are not dependent on the specific compilers.
*     2006/06/20: select 'A' if there is altLoc when reading PDB file.
*     2007/02/27: rotation matrix from Chain-1 to Chain-2 was added.
*     2007/04/18: added options with TM-score normalized by average
*                 length, shorter length, or longer length of two 
*                 structures.
*     2007/05/23: added additional output file 'TM.sup_all' for showing
*                 all atoms while 'TM.sup' is only for aligned atoms
*     2007/09/19: added a new feature alignment to deal with the problem
*                 of aligning fractional structures (e.g. protein
*                 interfaces).
*     2007/10/16: A bug for irregular bond-length models was fixed.
*     2009/03/14: A new initial alignment was added and previous initial
*                 alignments are further enhanced. This change increased
*                 accuracy by 4% but increases CPU cost by 40%.
*     2009/08/20: A bug for asymmetry alignment result was fixed.
*     2010/08/02: A new RMSD matrix was used to remove obsolete statements.
*                 Staled subroutines were deleted.
*     2011/01/03: The length of pdb file names were extended to 500.
*     2011/01/24: Fixed a bug on output file name created on 2011/01/03.
*     2011/01/30: An open source license is attached to the program.
*     2011/09/03: A new option "-d" is added to allow users to change
*                 TM-score normalization scale. A version number is attached 
*                 to the program from now on.
*     2011/10/11: A new scale (d0) was introduced for alignment search. This
*                 is to mainly improve alignment selection for small proteins 
*                 (e.g. L<50 residues) but also increase alignment coverage
*                 of larger proteins. Second, TM-align output format is changed
*                 and two TM-scores normalized by both chains are reported.
*     2011/10/12: Distance cutoff for gap is increased from 3.85A to 4.25A.
*                 Added 'TMalign -v' to allow user to check version number.
*     2012/01/24: Fix a bug for secondary structure definition
**************************************************************************
* The fortran program was converted to c using f2c. Manual editing was
* done to deprecate dependency on libf2c. Finally minor edits were made
* to make the code callable rather than being main, to return results
* in a data-struct and avoid alignment printing.
*/
////////
#include "tmalign_extern.h"
//
#define FALSE_ 0
#define TRUE_ 1

// macros
#ifndef max 
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min 
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

#ifndef dabs 
#define dabs( a ) ( ((a) < 0) ? (-1*a) : (a) )
#endif

#ifndef abs 
#define abs( a ) ( ((a) < 0) ? (-1*a) : (a) )
#endif

/* Common Block Declarations */
typedef struct {
	int mm1[5000], mm2[5000];
} initial4_;
//#define initial4_1 initial4_
static __thread initial4_ initial4_1; 

typedef struct {
	int nseq1, nseq2;
} length_;
//#define length_1 length_
static __thread length_ length_1;

__thread double w[5000];
//
typedef struct {
	float score[25000000] /* was [5000][5000] */, gap_open__;
	int invmap[5000];
} dpc_;
typedef struct {
	int invmap0[5000];
} alignrst_;
//#define alignrst_1 alignrst_
static __thread alignrst_ alignrst_1;
typedef struct {
	float d0, anseq;
} d0_;
//#define d0_1 d0_
static __thread d0_ d0_1;
typedef struct {
	float d0_min__;
} d0min_;
//#define d0min_1 d0min_
static __thread d0min_ d0min_1;
typedef struct {
	float d00, d002;
} d00_;
//#define d00_1 d00_
static __thread d00_ d00_1;
typedef struct {
	int invmap_i__[5000];
} init_;
//#define init_1 init_
static __thread init_ init_1;
typedef struct {
	float tm, tmmax;
} tm_;
//#define tm_1 tm_
static __thread tm_ tm_1;
typedef struct {
	float d8;
} d8_;
//#define d8_1 d8_
static __thread d8_ d8_1;
typedef struct {
	float zrms;
	int n_al__;
	float rmsd_al__;
} zscore_;
//#define zscore_1 zscore_
static __thread zscore_ zscore_1;
typedef struct {
	int isec[5000], jsec[5000];
} sec_;
//#define sec_1 sec_
static __thread sec_ sec_1;
//struct {
//        int invmap_a__[5000];
//} inv_;
//#define inv_1 inv_
typedef union {
	struct {
		float xt[5000], yt[5000], zt[5000], xb[5000], yb[5000], zb[5000];
	} _1;
	struct {
		float xa[5000], ya[5000], za[5000], xb[5000], yb[5000], zb[5000];
	} _2;
} stru_;
static __thread stru_ stru;
#define stru_1 stru._1
#define stru_2 stru._2
typedef struct {
	int nseqa, nseqb;
} nres_;
//#define nres_1 nres_
static __thread nres_ nres_1;
typedef struct {
	float d__, d0;
} para_;
//#define para_1 para_
static __thread para_ para_1;
typedef struct {
	int n_ali__, ia[5000], ib[5000];
} align_;
//#define align_1 align_
static __thread align_ align_1;
typedef struct {
	int i_ali__[5000], n_cut__;
} nscore_;
//#define nscore_1 nscore_
static __thread nscore_ nscore_1;
typedef struct {
	double score;
} scores_;
//#define scores_1 scores_
static __thread scores_ scores_1;
/* Table of constant values */
int c__1 = 1;
double c_b151 = .3;
double c_b152 = .33333333333333331;
int c__2 = 2;
int c__0 = 0;

//tmalign function defs
double pow_dd(double *a, double *b);
int pow_ii_f2c(int *ap, int *bp);
int u3b_(double *w, double *x, double *y, int *n, int *mode, double *rms,
		double *u, double *t, int *ier);
int dp_(int *nseq1, int *nseq2, dpc_ *dpc_1);
int score_fun__(void);
int tmscore_(float *dx, int *l1, float *x1, float *y1, float *z1, int *l2,
		float *x2, float *y2, float *z2, float *tm, float *rcomm, int *lcomm);
int score_fun8__(void);
int tmscore8_(float *dx, int *l1, float *x1, float *y1, float *z1, int *l2,
		float *x2, float *y2, float *z2, float *tm, float * rcomm, int *lcomm);
int tmscore8_search__(float *dx, int *l1, float *x1, float * y1, float *z1,
		int *l2, float *x2, float *y2, float *z2, float *tm, float *rcomm,
		int *lcomm);
int get_score1__(dpc_ *dpc_1, backbone_ *backbone_1);
int get_score__(dpc_ *dpc_1, backbone_ *backbone_1);
int get_gl__(float *gl, dpc_ *dpc_1, backbone_ *backbone_1);
int make_sec__(float *dis13, float *dis14, float *dis15, float *dis24,
		float * dis25, float *dis35);
double diszy_(int *i__, int *i1, int *i2, backbone_ *backbone_1);
int smooth_(void);
int get_initial5__(dpc_ *dpc_1, backbone_ *backbone_1);
int get_initial4__(dpc_ *dpc_1, backbone_ *backbone_1);
int get_initial3__(dpc_ *dpc_1, backbone_ *backbone_1);
int get_initial2__(dpc_ *dpc_1, backbone_ *backbone_1);
int get_initial__(dpc_ *dpc_1, backbone_ *backbone_1);
int tm_align__(backbone_ *backbone_1);
//

// tmalign specific variables
float rxx, ryy, rzz;
double distance1;
float d0_input__;
int i__, j, nseq, i__1;
float min_seq_len;
int m1[MAX_ATOMS], m2[MAX_ATOMS];
float tmscore_chain1, tmscore_chain2;
float tm8;
float dis2;
float xtm1[MAX_ATOMS], ytm1[MAX_ATOMS], ztm1[MAX_ATOMS], xtm2[MAX_ATOMS],
		ytm2[MAX_ATOMS], ztm2[MAX_ATOMS];
int number_aligned_residues;
int number_equal_residues;
float rmsd;
int number_aligned_residues_dis8;
int lcomm;
float rcomm;
float aligned_seq_identity;
int perform;

// internal function prototypes
int pow_ii_f2c(int *ap, int *bp);
double pow_dd(double *a, double *b);

int time_milli(struct timeb t) {
	return t.time * 1000 + t.millitm;
}

int diff_times(struct timeb t1, struct timeb t2) {
	return time_milli(t2) - time_milli(t1);
}

int main_process(char *ss1, char *ss2, char *seq1, char *seq2,
		ClientToMasterTransferBlock *ldata, backbone_ *backbone_1,
		int *atom_indices1, int *atom_indices2, int len1, int len2) {
	// set w to zero
	int cntr;
	for (cntr = 0; cntr < 5000; cntr++)
		w[cntr] = 1.;
	// set to global
	for (cntr = 0; cntr < MAX_ATOMS; cntr++) {
		initial4_1.mm1[cntr] = atom_indices1[cntr];
		initial4_1.mm2[cntr] = atom_indices2[cntr];
	}
	length_1.nseq1 = len1;
	length_1.nseq2 = len2;
	/* Initialized data */
	/* System generated locals */
	int i__1;
	float r__1, r__2, r__3;
	double d__1;
	/* Subroutine */
	/* Local variables */
	float d0_input__;
	int i__, j;
	float anseq_min__;
	int m1[5000], m2[5000];
	float tm1, tm2;
	float tm8;
	int m_d0__;
	float dis2;
	static __thread float xtm1[5000], ytm1[5000], ztm1[5000], xtm2[5000], ytm2[5000],
			ztm2[5000];
	int n_al__;
	int n_eq__;
	float rmsd;
	int n8_al__;
	int m_ave__, l_fix__, m_fix__, lcomm;
	float rcomm;
	float d0_fix__, seq_id__;

	/* ****** options -----------> */
	/* decided output */
	m_fix__ = -1;
	/* fixed length-scale only for output */
	m_ave__ = -1;
	/* using average length */
	/* diminum d0 for search */
	m_d0__ = -1;
	/* given d0 for output */
	i__ = 0;
	j = 0;

	++i__;
	anseq_min__ = (float) min(length_1.nseq1,length_1.nseq2);
	/* both search and d8_cut use nseq_min */d0_1.anseq = anseq_min__;
	/* length for defining TMscore in search */
	d__1 = (double) anseq_min__;
	d8_1.d8 = pow_dd(&d__1, &c_b151) * 1.5f + 3.5f;
	/* remove pairs with dis>d8 during search */
	if (d0_1.anseq > 19.f) {
		/* L=19, d0=0.168 */
		d__1 = (double) (d0_1.anseq - 15);
		d0_1.d0 = pow_dd(&d__1, &c_b152) * 1.24f - 1.8f;
		/* scale for defining TM-score */
	} else {
		d0_1.d0 = .168f;
	}
	d0min_1.d0_min__ = d0_1.d0 + .8f;
	/* best for search */
	if (d0_1.d0 < d0min_1.d0_min__) {
		d0_1.d0 = d0min_1.d0_min__;
	}
	/*     write(*,*)'d0 in search=',d0 */
	/* min d0 in search=0.968, min d0 in outpu */d00_1.d00 = d0_1.d0;
	/* for quickly calculate TM-score in searc */
	if (d00_1.d00 > 8.f) {
		d00_1.d00 = 8.f;
	}
	if (d00_1.d00 < 4.5f) {
		d00_1.d00 = 4.5f;
	}
	/* Computing 2nd power */
	r__1 = d00_1.d00;
	d00_1.d002 = r__1 * r__1;
	nseq = max(length_1.nseq1,length_1.nseq2);
	/* **** do alignment ************************** */
	//printf("do alignment\n");
	tm_align__(backbone_1);
	/* *********************************************************** */
	/* **   Refine alignment by cutting dis>d8 ------------------------> */
	/* to find invmap(j) */
	n_al__ = 0;
	i__1 = length_1.nseq2;
	for (j = 1; j <= i__1; ++j) {
		if (alignrst_1.invmap0[j - 1] > 0) {
			i__ = alignrst_1.invmap0[j - 1];
			++n_al__;
			xtm1[n_al__ - 1] = backbone_1->xa[i__ * 3 - 3];
			ytm1[n_al__ - 1] = backbone_1->xa[i__ * 3 - 2];
			ztm1[n_al__ - 1] = backbone_1->xa[i__ * 3 - 1];
			//printf("%f,%f,%f\n",backbone_1->xa[i__ * 3 - 3],backbone_1->xa[i__ * 3 - 2],
			//       backbone_1->xa[i__ * 3 - 1]);
			xtm2[n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 3];
			ytm2[n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 2];
			ztm2[n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 1];
			m1[n_al__ - 1] = i__;
			/* for recording residue order */
			m2[n_al__ - 1] = j;
		}
	}
	d0_input__ = d0_1.d0;
	/* scaled by seq_min */
	tmscore8_(&d0_input__, &n_al__, xtm1, ytm1, ztm1, &n_al__, xtm2, ytm2, ztm2,
			&tm_1.tm, &rcomm, &lcomm);
	/* TM-score with dis<d8 only */
	j = 0;
	n_eq__ = 0;
	i__1 = n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		/* Computing 2nd power */
		r__1 = xtm1[i__ - 1] - xtm2[i__ - 1];
		/* Computing 2nd power */
		r__2 = ytm1[i__ - 1] - ytm2[i__ - 1];
		/* Computing 2nd power */
		r__3 = ztm1[i__ - 1] - ztm2[i__ - 1];
		dis2 = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		if (dis2 <= d8_1.d8) {
			++j;
			xtm1[j - 1] = xtm1[i__ - 1];
			ytm1[j - 1] = ytm1[i__ - 1];
			ztm1[j - 1] = ztm1[i__ - 1];
			xtm2[j - 1] = xtm2[i__ - 1];
			ytm2[j - 1] = ytm2[i__ - 1];
			ztm2[j - 1] = ztm2[i__ - 1];
			m1[j - 1] = m1[i__ - 1];
			m2[j - 1] = m2[i__ - 1];
			if (strncmp(ss1 + (m1[i__ - 1] - 1) * 3,
					ss2 + (m2[i__ - 1] - 1) * 3, 3) == 0) {
				++n_eq__;
			}
		}
	}
	n8_al__ = j;
	seq_id__ = (float) n_eq__ / (n8_al__ + 1e-8f);
	/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	/* ^^^^^^^ alignment is done, all cutoffs were based on shorter chain^^^^^^^^ */
	/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	/* *********************************************************** */
	/* **   Output TM-score --------------------------> */
	/* only for showing residue-pair distance */d0min_1.d0_min__ = .5f;
	/*     Based on Chain_1===> */
	/* for TM-score output, consistent stdrd T */d0_1.anseq =
			(float) length_1.nseq1;
	if (d0_1.anseq > 21.f) {
		d__1 = (double) (d0_1.anseq - 15);
		d0_1.d0 = pow_dd(&d__1, &c_b152) * 1.24f - 1.8f;
		/* scale for defining TM-score */
	} else {
		d0_1.d0 = d0min_1.d0_min__;
	}
	if (d0_1.d0 < d0min_1.d0_min__) {
		d0_1.d0 = d0min_1.d0_min__;
	}
	d0_input__ = d0_1.d0;
	tmscore_(&d0_input__, &n8_al__, xtm1, ytm1, ztm1, &n8_al__, xtm2, ytm2,
			ztm2, &tm8, &rcomm, &lcomm);
	/* normal TMscore */
	rmsd = rcomm;
	tm1 = tm8 * n8_al__ / d0_1.anseq;
	/*     Based on Chain_2===> */d0_1.anseq = (float) length_1.nseq2;
	if (d0_1.anseq > 21.f) {
		d__1 = (double) (d0_1.anseq - 15);
		d0_1.d0 = pow_dd(&d__1, &c_b152) * 1.24f - 1.8f;
		/* scale for defining TM-score */
	} else {
		d0_1.d0 = d0min_1.d0_min__;
	}
	if (d0_1.d0 < d0min_1.d0_min__) {
		d0_1.d0 = d0min_1.d0_min__;
	}
	d0_input__ = d0_1.d0;
	tmscore_(&d0_input__, &n8_al__, xtm1, ytm1, ztm1, &n8_al__, xtm2, ytm2,
			ztm2, &tm8, &rcomm, &lcomm);
	/* normal TMscore */
	tm2 = tm8 * n8_al__ / d0_1.anseq;
	/*     Based on Average length===> */
	if (m_ave__ == 1) {
		d0_1.anseq = (length_1.nseq1 + length_1.nseq2) / 2.f;
		if (d0_1.anseq > 21.f) {
			d__1 = (double) (d0_1.anseq - 15);
			d0_1.d0 = pow_dd(&d__1, &c_b152) * 1.24f - 1.8f;
			/* scale for defining TM-sco */
		} else {
			d0_1.d0 = d0min_1.d0_min__;
		}
		if (d0_1.d0 < d0min_1.d0_min__) {
			d0_1.d0 = d0min_1.d0_min__;
		}
		d0_input__ = d0_1.d0;
		tmscore_(&d0_input__, &n8_al__, xtm1, ytm1, ztm1, &n8_al__, xtm2, ytm2,
				ztm2, &tm8, &rcomm, &lcomm);
		/* normal TMscore */
	}
	/*     Based on assigned length===> */
	if (m_fix__ == 1) {
		d0_1.anseq = (float) l_fix__;
		/* input length */
		if (d0_1.anseq > 21.f) {
			d__1 = (double) (d0_1.anseq - 15);
			d0_1.d0 = pow_dd(&d__1, &c_b152) * 1.24f - 1.8f;
			/* scale for defining TM-sco */
		} else {
			d0_1.d0 = d0min_1.d0_min__;
		}
		if (d0_1.d0 < d0min_1.d0_min__) {
			d0_1.d0 = d0min_1.d0_min__;
		}
		d0_input__ = d0_1.d0;
		tmscore_(&d0_input__, &n8_al__, xtm1, ytm1, ztm1, &n8_al__, xtm2, ytm2,
				ztm2, &tm8, &rcomm, &lcomm);
		/* normal TMscore */
	}
	/*     Based on user-specified d0===> */
	if (m_d0__ == 1) {
		d0_1.d0 = d0_fix__;
		d0_input__ = d0_1.d0;
		tmscore_(&d0_input__, &n8_al__, xtm1, ytm1, ztm1, &n8_al__, xtm2, ytm2,
				ztm2, &tm8, &rcomm, &lcomm);
		/* normal TMscore */
	}
	/* ******** for output summary ****************************** */
	/*printf("Name of Chain_1: %s\n",lfname1);
	 printf("Name of Chain_2: %s\n",lfname2);*/

#ifdef VERBOSE_PRINTS
	struct timeb tp1;

	printf("Length of Chain_1: %d\n",length_1.nseq1);
	printf("Length of Chain_2: %d\n\n",length_1.nseq2);
	printf("Aligned length= %d, RMSD= %.2f, Seq_ID=n_identical/n_aligned= %.3f\n\n",
			n8_al__, rmsd, seq_id__);
	printf("TM1: %f, TM2: %f\n", tm1, tm2);
	printf("RMSD1: %f, RMSD2: %f\n", rmsd, rcomm);ftime(&tp1);
	printf("*****intmalign_time3 = %ld.%d\n",tp1.time,tp1.millitm);
#endif

	ldata->tm1 = tm1;
	ldata->tm2 = tm2;
	ldata->rmsd = rmsd;
	ldata->aln_len = n8_al__;
	ldata->seq_id = seq_id__;
	ldata->len1 = length_1.nseq1;
	ldata->len2 = length_1.nseq2;
	ldata->n_eq__ = n_eq__;

	return 0;
}

/* ********************************************************************** */
/* ********************************************************************** */
/*     Structure superposition */
/* ********************************************************************** */
/* ********************************************************************** */
/* ********************************************************************** */
/* Subroutine */int tm_align__(backbone_ *backbone_1) {
	dpc_ *dpc_1 = (dpc_ *) malloc(sizeof(dpc_));
	/* System generated locals */
	int i__1, i__2;
	float r__1;
	/* Local variables */
	int i__, j;
	int id;
	float ddcc, diff;
	float gapp[100];
	int i_gapp__, n_gapp__;
	float tm_old__;
	tm_1.tmmax = 0.f;
	n_gapp__ = 2;
	gapp[0] = -.6f;
	gapp[1] = 0.f;
	ddcc = .4f;
	if (d0_1.anseq <= 40.f) {
		ddcc = .1f;
	}
	/* 11111111111111111111111111111111111111111111111111111111 */
	/*     get initial alignment from gapless threading */
	/* ********************************************************* */
	get_initial__(dpc_1, backbone_1);
	/* gapless threading */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dpc_1->invmap[i__ - 1] = init_1.invmap_i__[i__ - 1];
		/* with highest zcore */
	}
	get_score__(dpc_1, backbone_1);
	/* TM, matrix score(i,j) */
	if (tm_1.tm > tm_1.tmmax) {
		tm_1.tmmax = tm_1.tm;
		i__1 = length_1.nseq2;
		for (j = 1; j <= i__1; ++j) {
			alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
		}
	}
	/* **************************************************************** */
	/*       initerative alignment, for different gap_open: */
	/* **************************************************************** */
	i__1 = n_gapp__;
	for (i_gapp__ = 1; i_gapp__ <= i__1; ++i_gapp__) {
		/* different gap panalties */
		dpc_1->gap_open__ = gapp[i_gapp__ - 1];
		/* gap panalty */
		for (id = 1; id <= 30; ++id) {
			/* maximum interation is 200 */
			dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
			/*     Input: score(i,j), and gap_open */
			/*     Output: invmap(j) */
			/* produce alignment invmap(j) */
			get_score__(dpc_1, backbone_1);
			/*     record the best alignment in whole search ----------> */
			/* calculate TM-score, score(i,j) */
			if (tm_1.tm > tm_1.tmmax) {
				tm_1.tmmax = tm_1.tm;
				i__2 = length_1.nseq2;
				for (j = 1; j <= i__2; ++j) {
					alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
				}
			}
			if (id > 1) {
				diff = (r__1 = tm_1.tm - tm_old__, dabs(r__1));
				if (diff < 1e-6f) {
					goto L111;
				}
			}
			tm_old__ = tm_1.tm;
			/* L11: */
		}
		L111:
		/* L1: */
		;
	}
	/* 222222222222222222222222222222222222222222222222222222222 */
	/*     get initial alignment from secondary structure alignment */
	/* ********************************************************* */
	get_initial2__(dpc_1, backbone_1);
	/* DP for secondary structure */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dpc_1->invmap[i__ - 1] = init_1.invmap_i__[i__ - 1];
		/* with highest zcore */
	}
	get_score__(dpc_1, backbone_1);
	/* TM, score(i,j) */
	if (tm_1.tm > tm_1.tmmax) {
		tm_1.tmmax = tm_1.tm;
		i__1 = length_1.nseq2;
		for (j = 1; j <= i__1; ++j) {
			alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
		}
	}
	if (tm_1.tm <= tm_1.tmmax * .2f) {
		goto L2222;
	}
	/* **************************************************************** */
	/*     initerative alignment, for different gap_open: */
	/* **************************************************************** */
	i__1 = n_gapp__;
	for (i_gapp__ = 1; i_gapp__ <= i__1; ++i_gapp__) {
		/* different gap panalties */
		dpc_1->gap_open__ = gapp[i_gapp__ - 1];
		/* gap panalty */
		for (id = 1; id <= 30; ++id) {
			/* maximum interation is 200 */
			dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
			/*     Input: score(i,j), and gap_open */
			/*     Output: invmap(j) */
			/* produce alignment invmap(j) */
			get_score__(dpc_1, backbone_1);
			/*     write(*,21)gap_open,rmsd_al,n_al,TM */
			/*     record the best alignment in whole search ----------> */
			/* calculate TM-score, score(i,j) */
			if (tm_1.tm > tm_1.tmmax) {
				tm_1.tmmax = tm_1.tm;
				i__2 = length_1.nseq2;
				for (j = 1; j <= i__2; ++j) {
					alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
				}
			}
			if (id > 1) {
				diff = (r__1 = tm_1.tm - tm_old__, dabs(r__1));
				if (diff < 1e-6f) {
					goto L222;
				}
			}
			tm_old__ = tm_1.tm;
			/* L22: */
		}
		L222:
		/* L2: */
		;
	}
	L2222:
	/* 555555555555555555555555555555555555555555555555555555555555555555 */
	/*     get initial alignment of local structure superposition */
	/* ****************************************************************** */
	get_initial5__(dpc_1, backbone_1);
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dpc_1->invmap[i__ - 1] = init_1.invmap_i__[i__ - 1];
		/* with highest zcore */
	}
	get_score__(dpc_1, backbone_1);
	/* TM, matrix score(i,j) */
	if (tm_1.tm > tm_1.tmmax) {
		tm_1.tmmax = tm_1.tm;
		i__1 = length_1.nseq2;
		for (j = 1; j <= i__1; ++j) {
			alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
		}
	}
	if (tm_1.tm <= tm_1.tmmax * ddcc) {
		goto L5555;
	}
	/* **************************************************************** */
	/*     initerative alignment, for different gap_open: */
	/* **************************************************************** */
	i__1 = n_gapp__;
	for (i_gapp__ = 1; i_gapp__ <= i__1; ++i_gapp__) {
		/* different gap panalties */
		dpc_1->gap_open__ = gapp[i_gapp__ - 1];
		/* gap panalty */
		for (id = 1; id <= 2; ++id) {
			/* maximum interation is 200 */
			dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
			/*     Input: score(i,j), and gap_open */
			/*     Output: invmap(j) */
			/* produce alignment invmap(j) */
			get_score__(dpc_1, backbone_1);
			/*     record the best alignment in whole search ----------> */
			/* calculate TM-score, score(i,j) */
			if (tm_1.tm > tm_1.tmmax) {
				tm_1.tmmax = tm_1.tm;
				i__2 = length_1.nseq2;
				for (j = 1; j <= i__2; ++j) {
					alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
				}
			}
			if (id > 1) {
				diff = (r__1 = tm_1.tm - tm_old__, dabs(r__1));
				if (diff < 1e-6f) {
					goto L555;
				}
			}
			tm_old__ = tm_1.tm;
			/* L55: */
		}
		L555:
		/* L5: */
		;
	}
	L5555:
	/* 333333333333333333333333333333333333333333333333333333333333 */
	/*     get initial alignment from invmap0+SS */
	/* ************************************************************ */
	get_initial3__(dpc_1, backbone_1);
	/* invmap0+SS */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dpc_1->invmap[i__ - 1] = init_1.invmap_i__[i__ - 1];
		/* with highest zcore */
	}
	get_score__(dpc_1, backbone_1);
	/* TM, score(i,j) */
	if (tm_1.tm > tm_1.tmmax) {
		tm_1.tmmax = tm_1.tm;
		i__1 = length_1.nseq2;
		for (j = 1; j <= i__1; ++j) {
			alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
		}
	}
	if (tm_1.tm <= tm_1.tmmax * ddcc) {
		goto L3333;
	}
	/* **************************************************************** */
	/*     initerative alignment, for different gap_open: */
	/* **************************************************************** */
	i__1 = n_gapp__;
	for (i_gapp__ = 1; i_gapp__ <= i__1; ++i_gapp__) {
		/* different gap panalties */
		dpc_1->gap_open__ = gapp[i_gapp__ - 1];
		/* gap panalty */
		for (id = 1; id <= 30; ++id) {
			/* maximum interation is 200 */
			dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
			/*     Input: score(i,j), and gap_open */
			/*     Output: invmap(j) */
			/* produce alignment invmap(j) */
			get_score__(dpc_1, backbone_1);
			/*     write(*,21)gap_open,rmsd_al,n_al,TM */
			/*     record the best alignment in whole search ----------> */
			/* calculate TM-score, score(i,j) */
			if (tm_1.tm > tm_1.tmmax) {
				tm_1.tmmax = tm_1.tm;
				i__2 = length_1.nseq2;
				for (j = 1; j <= i__2; ++j) {
					alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
				}
			}
			if (id > 1) {
				diff = (r__1 = tm_1.tm - tm_old__, dabs(r__1));
				if (diff < 1e-6f) {
					goto L333;
				}
			}
			tm_old__ = tm_1.tm;
			/* L33: */
		}
		L333:
		/* L3: */
		;
	}
	L3333:
	/* 444444444444444444444444444444444444444444444444444444444 */
	/*     get initial alignment of pieces from gapless threading */
	/* ********************************************************* */
	get_initial4__(dpc_1, backbone_1);
	/* gapless threading */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dpc_1->invmap[i__ - 1] = init_1.invmap_i__[i__ - 1];
		/* with highest zcore */
	}
	get_score__(dpc_1, backbone_1);
	/* TM, matrix score(i,j) */
	if (tm_1.tm > tm_1.tmmax) {
		tm_1.tmmax = tm_1.tm;
		i__1 = length_1.nseq2;
		for (j = 1; j <= i__1; ++j) {
			alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
		}
	}
	if (tm_1.tm <= tm_1.tmmax * ddcc) {
		goto L4444;
	}
	/* **************************************************************** */
	/*     initerative alignment, for different gap_open: */
	/* **************************************************************** */
	i__1 = n_gapp__;
	for (i_gapp__ = 2; i_gapp__ <= i__1; ++i_gapp__) {
		/* different gap panalties */
		dpc_1->gap_open__ = gapp[i_gapp__ - 1];
		/* gap panalty */
		for (id = 1; id <= 2; ++id) {
			/* maximum interation is 200 */
			dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
			/*     Input: score(i,j), and gap_open */
			/*     Output: invmap(j) */
			/* produce alignment invmap(j) */
			get_score__(dpc_1, backbone_1);
			/*     record the best alignment in whole search ----------> */
			/* calculate TM-score, score(i,j) */
			if (tm_1.tm > tm_1.tmmax) {
				tm_1.tmmax = tm_1.tm;
				i__2 = length_1.nseq2;
				for (j = 1; j <= i__2; ++j) {
					alignrst_1.invmap0[j - 1] = dpc_1->invmap[j - 1];
				}
			}
			/* L44: */
		}
		/* L4: */
	}
	L4444:
	/* ^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^ */
	delete[] dpc_1;
	return 0;
} /* tm_align__ */
/* ************************************************************* */
/*     get initial alignment invmap0(i) from gapless threading */
/* ************************************************************* */
/* Subroutine */int get_initial__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* System generated locals */
	int i__1, i__2;
	/* Local variables */
	int i__, j, l, n1, n2;
	float al, gl;
	int idel;
	float gl_max__;
	int ishift;
	al = (float) min(length_1.nseq1,length_1.nseq2);
	idel = al / 2.f;
	/* minimum size of considered fragment */
	if (idel <= 5) {
		idel = 5;
	}
	n1 = -length_1.nseq2 + idel;
	n2 = length_1.nseq1 - idel;
	gl_max__ = 0.f;
	i__1 = n2;
	for (ishift = n1; ishift <= i__1; ++ishift) {
		l = 0;
		i__2 = length_1.nseq2;
		for (j = 1; j <= i__2; ++j) {
			i__ = j + ishift;
			if (i__ >= 1 && i__ <= length_1.nseq1) {
				++l;
				dpc_1->invmap[j - 1] = i__;
			} else {
				dpc_1->invmap[j - 1] = -1;
			}
		}
		if (l >= idel) {
			get_gl__(&gl, dpc_1, backbone_1);
			if (gl > gl_max__) {
				gl_max__ = gl;
				i__2 = length_1.nseq2;
				for (i__ = 1; i__ <= i__2; ++i__) {
					init_1.invmap_i__[i__ - 1] = dpc_1->invmap[i__ - 1];
				}
			}
		}
	}
	return 0;
} /* get_initial__ */
/* ************************************************************* */
/*     get initial alignment invmap0(i) from secondary structure */
/* ************************************************************* */
/* Subroutine */int get_initial2__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* System generated locals */
	int i__1, i__2;
	/* Local variables */
	int i__, j, j1, j2, j3, j4, j5;
	float dis13, dis14, dis15, dis24, dis25, dis35;
	/* ********* assign secondary structures *************** */
	/*     1->coil, 2->helix, 3->turn, 4->strand */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sec_1.isec[i__ - 1] = 1;
		j1 = i__ - 2;
		j2 = i__ - 1;
		j3 = i__;
		j4 = i__ + 1;
		j5 = i__ + 2;
		if (j1 >= 1 && j5 <= length_1.nseq1) {
			dis13 = diszy_(&c__0, &j1, &j3, backbone_1);
			dis14 = diszy_(&c__0, &j1, &j4, backbone_1);
			dis15 = diszy_(&c__0, &j1, &j5, backbone_1);
			dis24 = diszy_(&c__0, &j2, &j4, backbone_1);
			dis25 = diszy_(&c__0, &j2, &j5, backbone_1);
			dis35 = diszy_(&c__0, &j3, &j5, backbone_1);
			sec_1.isec[i__ - 1] = make_sec__(&dis13, &dis14, &dis15, &dis24,
					&dis25, &dis35);
		}
	}
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sec_1.jsec[i__ - 1] = 1;
		j1 = i__ - 2;
		j2 = i__ - 1;
		j3 = i__;
		j4 = i__ + 1;
		j5 = i__ + 2;
		if (j1 >= 1 && j5 <= length_1.nseq2) {
			dis13 = diszy_(&c__1, &j1, &j3, backbone_1);
			dis14 = diszy_(&c__1, &j1, &j4, backbone_1);
			dis15 = diszy_(&c__1, &j1, &j5, backbone_1);
			dis24 = diszy_(&c__1, &j2, &j4, backbone_1);
			dis25 = diszy_(&c__1, &j2, &j5, backbone_1);
			dis35 = diszy_(&c__1, &j3, &j5, backbone_1);
			sec_1.jsec[i__ - 1] = make_sec__(&dis13, &dis14, &dis15, &dis24,
					&dis25, &dis35);
		}
	}
	smooth_();
	/* ********* score matrix ************************** */
	/* smooth the assignment */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = length_1.nseq2;
		for (j = 1; j <= i__2; ++j) {
			if (sec_1.isec[i__ - 1] == sec_1.jsec[j - 1]) {
				dpc_1->score[i__ + j * 5000 - 5001] = 1.f;
			} else {
				dpc_1->score[i__ + j * 5000 - 5001] = 0.f;
			}
		}
	}
	/* ********* find initial alignment: invmap(j) ************ */
	dpc_1->gap_open__ = -1.f;
	/* should be -1 */
	dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
	/* produce alignment invmap(j) */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		init_1.invmap_i__[i__ - 1] = dpc_1->invmap[i__ - 1];
	}
	/* ^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* get_initial2__ */
/* ************************************************************* */
/*     get initial alignment invmap0(i) from secondary structure */
/*     and previous alignments */
/* ************************************************************* */
/* Subroutine */int get_initial3__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* System generated locals */
	int i__1, i__2;
	/* Local variables */
	int i__, j;
	/* ********* score matrix ************************** */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dpc_1->invmap[i__ - 1] = alignrst_1.invmap0[i__ - 1];
	}
	get_score1__(dpc_1, backbone_1);
	/* get score(i,j) using RMSD martix */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		i__2 = length_1.nseq2;
		for (j = 1; j <= i__2; ++j) {
			if (sec_1.isec[i__ - 1] == sec_1.jsec[j - 1]) {
				dpc_1->score[i__ + j * 5000 - 5001] += .5f;
			} else {
				dpc_1->score[i__ + j * 5000 - 5001] = dpc_1->score[i__
						+ j * 5000 - 5001];
			}
		}
	}
	/* ********* find initial alignment: invmap(j) ************ */
	dpc_1->gap_open__ = -1.f;
	/* should be -1 */
	dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
	/* produce alignment invmap(j) */
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		init_1.invmap_i__[i__ - 1] = dpc_1->invmap[i__ - 1];
	}
	/* ^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* get_initial3__ */
/* ************************************************************* */
/*     get initial alignment invmap0(i) from fragment gapless threading */
/* ************************************************************* */
/* Subroutine */int get_initial4__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* System generated locals */
	int i__1, i__2, i__3;
	/* Local variables */
	int i__, j, k, l, n1, n2;
	float al, gl;
	int mm[10000] /* was [2][5000] */;
	float dcu, dis;
	int ifr[5000], nfr;
	float dcu0;
	int lfr2[10000]
	/* was [2][5000] */, idel, l_fr__, mark, i_fr2__[2], nseq0;
	int *ifr2 = (int *)malloc(sizeof(int)*50000000);
	float r_min__;
	int nseq1___, nseq2___;
	float gl_max__;
	int ishift;
	int contin;
	float fra_min__;
	int lfr_max__;
	float fra_min1__;
	fra_min__ = 4.f;
	/* >=4,minimum fragment for search */
	fra_min1__ = fra_min__ - 1;
	/* cutoff for shift, save time */
	dcu0 = 4.25f;
	/* cc   Find the smallest continuous fragments --------> */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		mm[(i__ << 1) - 2] = initial4_1.mm1[i__ - 1];
	}
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		mm[(i__ << 1) - 1] = initial4_1.mm2[i__ - 1];
	}
	for (k = 1; k <= 2; ++k) {
		dcu = dcu0;
		if (k == 1) {
			nseq0 = length_1.nseq1;
			r_min__ = length_1.nseq1 / 3.f;
			/* minimum fragment, in case too small pro */
		} else {
			nseq0 = length_1.nseq2;
			r_min__ = length_1.nseq2 / 3.f;
			/* minimum fragment, in case too small pro */
		}
		if (r_min__ > fra_min__) {
			r_min__ = fra_min__;
		}
		L20: nfr = 1;
		/* number of fragments */
		j = 1;
		/* number of residue at nf-fragment */
		ifr2[k + (nfr + j * 5000 << 1) - 10003] = 1;
		/* what residue */
		lfr2[k + (nfr << 1) - 3] = j;
		/* length of the fragment */
		i__1 = nseq0;
		for (i__ = 2; i__ <= i__1; ++i__) {
			i__2 = k - 1;
			i__3 = i__ - 1;
			dis = diszy_(&i__2, &i__3, &i__, backbone_1);
			contin = FALSE_;
			if (dcu > dcu0) {
				if (dis < dcu) {
					contin = TRUE_;
				}
			} else if (mm[k + (i__ << 1) - 3]
					== mm[k + (i__ - 1 << 1) - 3] + 1) {
				if (dis < dcu) {
					contin = TRUE_;
				}
			}
			if (contin) {
				++j;
				ifr2[k + (nfr + j * 5000 << 1) - 10003] = i__;
				lfr2[k + (nfr << 1) - 3] = j;
			} else {
				++nfr;
				j = 1;
				ifr2[k + (nfr + j * 5000 << 1) - 10003] = i__;
				lfr2[k + (nfr << 1) - 3] = j;
			}
		}
		lfr_max__ = 0;
		i_fr2__[k - 1] = 1;
		/* ID of the maximum piece */
		i__1 = nfr;
		for (i__ = 1; i__ <= i__1; ++i__) {
			if (lfr_max__ < lfr2[k + (i__ << 1) - 3]) {
				lfr_max__ = lfr2[k + (i__ << 1) - 3];
				i_fr2__[k - 1] = i__;
			}
		}
		if ((float) lfr_max__ < r_min__) {
			dcu += .01f;
			goto L20;
		}
	}
	/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	/* cc   select what piece will be used (this may araise ansysmetry, but */
	/* cc   only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1 */
	/* cc   if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1 */
	mark = 1;
	if (lfr2[(i_fr2__[0] << 1) - 2] < lfr2[(i_fr2__[1] << 1) - 1]) {
		mark = 1;
	} else if (lfr2[(i_fr2__[0] << 1) - 2] > lfr2[(i_fr2__[1] << 1) - 1]) {
		mark = 2;
	} else {
		/* Lfr1=Lfr2 */
		if (length_1.nseq1 <= length_1.nseq2) {
			mark = 1;
		} else {
			mark = 2;
		}
	}
	/* cc */
	l_fr__ = lfr2[mark + (i_fr2__[mark - 1] << 1) - 3];
	i__1 = l_fr__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		ifr[i__ - 1] =
				ifr2[mark + (i_fr2__[mark - 1] + i__ * 5000 << 1) - 10003];
	}
	/* cc */
	if (mark == 1) {
		/* non-redundant to get_initial1 */
		nseq0 = length_1.nseq1;
	} else {
		nseq0 = length_1.nseq2;
	}
	if (l_fr__ == nseq0) {
		n1 = (int) (nseq0 * .1f);
		/* 0 */
		n2 = (int) (nseq0 * .89f);
		/* 2 */
		j = 0;
		i__1 = n2;
		for (i__ = n1; i__ <= i__1; ++i__) {
			++j;
			ifr[j - 1] = ifr[n1 + j - 1];
		}
		l_fr__ = j;
	}
	/* cc   get initial -------------> */
	if (mark == 1) {
		/* nseq1 as the smallest one */
		nseq1___ = l_fr__;
		al = (float) min(nseq1___,length_1.nseq2);
		idel = al / 2.5f;
		/* minimum size of considered fragment */
		if ((float) idel <= fra_min1__) {
			idel = fra_min1__;
		}
		n1 = -length_1.nseq2 + idel;
		/* shift1 */
		n2 = nseq1___ - idel;
		/* shift2 */
		gl_max__ = 0.f;
		i__1 = n2;
		for (ishift = n1; ishift <= i__1; ++ishift) {
			l = 0;
			i__2 = length_1.nseq2;
			for (j = 1; j <= i__2; ++j) {
				i__ = j + ishift;
				if (i__ >= 1 && i__ <= nseq1___) {
					++l;
					dpc_1->invmap[j - 1] = ifr[i__ - 1];
				} else {
					dpc_1->invmap[j - 1] = -1;
				}
			}
			if (l >= idel) {
				get_gl__(&gl, dpc_1, backbone_1);
				if (gl > gl_max__) {
					gl_max__ = gl;
					i__2 = length_1.nseq2;
					for (i__ = 1; i__ <= i__2; ++i__) {
						init_1.invmap_i__[i__ - 1] = dpc_1->invmap[i__ - 1];
					}
				}
			}
		}
	} else {
		nseq2___ = l_fr__;
		al = (float) min(length_1.nseq1,nseq2___);
		idel = al / 2.5f;
		/* minimum size of considered fragment */
		if ((float) idel <= fra_min1__) {
			idel = fra_min1__;
		}
		n1 = -nseq2___ + idel;
		n2 = length_1.nseq1 - idel;
		gl_max__ = 0.f;
		i__1 = n2;
		for (ishift = n1; ishift <= i__1; ++ishift) {
			l = 0;
			i__2 = length_1.nseq2;
			for (j = 1; j <= i__2; ++j) {
				dpc_1->invmap[j - 1] = -1;
			}
			i__2 = nseq2___;
			for (j = 1; j <= i__2; ++j) {
				i__ = j + ishift;
				if (i__ >= 1 && i__ <= length_1.nseq1) {
					++l;
					dpc_1->invmap[ifr[j - 1] - 1] = i__;
				}
			}
			if (l >= idel) {
				get_gl__(&gl, dpc_1, backbone_1);
				if (gl > gl_max__) {
					gl_max__ = gl;
					i__2 = length_1.nseq2;
					for (i__ = 1; i__ <= i__2; ++i__) {
						init_1.invmap_i__[i__ - 1] = dpc_1->invmap[i__ - 1];
					}
				}
			}
		}
	}
	free(ifr2);
	return 0;
} /* get_initial4__ */
/* ************************************************************* */
/*    fifth initial alignement. Using local structure super   * */
/*    position.                                               * */
/* ************************************************************* */
/* Subroutine */int get_initial5__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* Initialized data */
	/* System generated locals */
	int i__1, i__2, i__3, i__4, i__5, i__6;
	float r__1, r__2, r__3;
	/* Local variables */
	int k;
	double t[3], u[9] /* was [3][3] */;
	int i1, j1, m1, m2;
	float d01, d02, dd;
	int al, ii;
	float gl;
	int jj;
	int ns;
	float xx, yy, zz;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	int iii, jjj, ier;
	double rms;
	int n_frag__;
	float glmaxa;
	/* **** setting parameters ************************************ */
	d01 = d0_1.d0 + 1.5f;
	if (d01 < d0min_1.d0_min__) {
		d01 = d0min_1.d0_min__;
	}
	d02 = d01 * d01;
	glmaxa = 0.f;
	al = min(length_1.nseq1,length_1.nseq2);
	if (al > 250) {
		n_frag__ = 50;
	} else if (al > 200) {
		n_frag__ = 40;
	} else if (al > 150) {
		n_frag__ = 30;
	} else {
		n_frag__ = 20;
		/* length of fragment for superposition */
	}
	if (n_frag__ > al / 3) {
		n_frag__ = al / 3;
	}
	ns = 20;
	/* tail length to discard */
	if (ns > al / 3) {
		ns = al / 3;
	}
	m1 = length_1.nseq1 - n_frag__ - ns;
	m2 = length_1.nseq2 - n_frag__ - ns;
	i__1 = m1;
	i__2 = n_frag__;
	for (ii = ns; i__2 < 0 ? ii >= i__1 : ii <= i__1; ii += i__2) {
		i__3 = m2;
		i__4 = n_frag__;
		for (jj = ns; i__4 < 0 ? jj >= i__3 : jj <= i__3; jj += i__4) {
			i__5 = n_frag__;
			for (k = 1; k <= i__5; ++k) {
				iii = ii + k - 1;
				jjj = jj + k - 1;
				r_1__[k * 3 - 3] = backbone_1->xa[iii * 3 - 3];
				r_1__[k * 3 - 2] = backbone_1->xa[iii * 3 - 2];
				r_1__[k * 3 - 1] = backbone_1->xa[iii * 3 - 1];
				r_2__[k * 3 - 3] = backbone_1->xa[(jjj + 5000) * 3 - 3];
				r_2__[k * 3 - 2] = backbone_1->xa[(jjj + 5000) * 3 - 2];
				r_2__[k * 3 - 1] = backbone_1->xa[(jjj + 5000) * 3 - 1];
			}
			/* ********superpose the two structures and rotate it ***************** */
			u3b_(w, r_1__, r_2__, &n_frag__, &c__1, &rms, u, t, &ier);
			/* u rotate r_1 to r_ */
			i__5 = length_1.nseq1;
			for (i1 = 1; i1 <= i__5; ++i1) {
				xx = t[0] + u[0] * backbone_1->xa[i1 * 3 - 3]
						+ u[3] * backbone_1->xa[i1 * 3 - 2]
						+ u[6] * backbone_1->xa[i1 * 3 - 1];
				yy = t[1] + u[1] * backbone_1->xa[i1 * 3 - 3]
						+ u[4] * backbone_1->xa[i1 * 3 - 2]
						+ u[7] * backbone_1->xa[i1 * 3 - 1];
				zz = t[2] + u[2] * backbone_1->xa[i1 * 3 - 3]
						+ u[5] * backbone_1->xa[i1 * 3 - 2]
						+ u[8] * backbone_1->xa[i1 * 3 - 1];
				i__6 = length_1.nseq2;
				for (j1 = 1; j1 <= i__6; ++j1) {
					/* Computing 2nd power */
					r__1 = xx - backbone_1->xa[(j1 + 5000) * 3 - 3];
					/* Computing 2nd power */
					r__2 = yy - backbone_1->xa[(j1 + 5000) * 3 - 2];
					/* Computing 2nd power */
					r__3 = zz - backbone_1->xa[(j1 + 5000) * 3 - 1];
					dd = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
					dpc_1->score[i1 + j1 * 5000 - 5001] = 1 / (dd / d02 + 1);
					/* changing */
				}
			}
			/* ********extract alignement with score(i,j) ***************** */
			dp_(&length_1.nseq1, &length_1.nseq2, dpc_1);
			get_gl__(&gl, dpc_1, backbone_1);
			if (gl > glmaxa) {
				glmaxa = gl;
				i__5 = length_1.nseq2;
				for (j1 = 1; j1 <= i__5; ++j1) {
					init_1.invmap_i__[j1 - 1] = dpc_1->invmap[j1 - 1];
				}
			}
		}
	}
	return 0;
} /* get_initial5__ */
/* ************************************************************* */
/*     smooth the secondary structure assignment */
/* ************************************************************* */
/* Subroutine */int smooth_(void) {
	/* System generated locals */
	int i__1;
	/* Local variables */
	int i__, j;
	/* **   smooth single --------------> */
	/* **   --x-- => ----- */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (sec_1.isec[i__ - 1] == 2 || sec_1.isec[i__ - 1] == 4) {
			j = sec_1.isec[i__ - 1];
			if (sec_1.isec[i__ - 3] != j) {
				if (sec_1.isec[i__ - 2] != j) {
					if (sec_1.isec[i__] != j) {
						if (sec_1.isec[i__ + 1] != j) {
							sec_1.isec[i__ - 1] = 1;
						}
					}
				}
			}
		}
	}
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (sec_1.jsec[i__ - 1] == 2 || sec_1.jsec[i__ - 1] == 4) {
			j = sec_1.jsec[i__ - 1];
			if (sec_1.jsec[i__ - 3] != j) {
				if (sec_1.jsec[i__ - 2] != j) {
					if (sec_1.jsec[i__] != j) {
						if (sec_1.jsec[i__ + 1] != j) {
							sec_1.jsec[i__ - 1] = 1;
						}
					}
				}
			}
		}
	}
	/* **   smooth double --------------> */
	/* **   --xx-- => ------ */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (sec_1.isec[i__ - 1] != 2) {
			if (sec_1.isec[i__] != 2) {
				if (sec_1.isec[i__ + 1] == 2) {
					if (sec_1.isec[i__ + 2] == 2) {
						if (sec_1.isec[i__ + 3] != 2) {
							if (sec_1.isec[i__ + 4] != 2) {
								sec_1.isec[i__ + 1] = 1;
								sec_1.isec[i__ + 2] = 1;
							}
						}
					}
				}
			}
		}
		if (sec_1.isec[i__ - 1] != 4) {
			if (sec_1.isec[i__] != 4) {
				if (sec_1.isec[i__ + 1] == 4) {
					if (sec_1.isec[i__ + 2] == 4) {
						if (sec_1.isec[i__ + 3] != 4) {
							if (sec_1.isec[i__ + 4] != 4) {
								sec_1.isec[i__ + 1] = 1;
								sec_1.isec[i__ + 2] = 1;
							}
						}
					}
				}
			}
		}
	}
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (sec_1.jsec[i__ - 1] != 2) {
			if (sec_1.jsec[i__] != 2) {
				if (sec_1.jsec[i__ + 1] == 2) {
					if (sec_1.jsec[i__ + 2] == 2) {
						if (sec_1.jsec[i__ + 3] != 2) {
							if (sec_1.jsec[i__ + 4] != 2) {
								sec_1.jsec[i__ + 1] = 1;
								sec_1.jsec[i__ + 2] = 1;
							}
						}
					}
				}
			}
		}
		if (sec_1.jsec[i__ - 1] != 4) {
			if (sec_1.jsec[i__] != 4) {
				if (sec_1.jsec[i__ + 1] == 4) {
					if (sec_1.jsec[i__ + 2] == 4) {
						if (sec_1.jsec[i__ + 3] != 4) {
							if (sec_1.jsec[i__ + 4] != 4) {
								sec_1.jsec[i__ + 1] = 1;
								sec_1.jsec[i__ + 2] = 1;
							}
						}
					}
				}
			}
		}
	}
	/* **   connect --------------> */
	/* **   x-x => xxx */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (sec_1.isec[i__ - 1] == 2) {
			if (sec_1.isec[i__] != 2) {
				if (sec_1.isec[i__ + 1] == 2) {
					sec_1.isec[i__] = 2;
				}
			}
		}
		if (sec_1.isec[i__ - 1] == 4) {
			if (sec_1.isec[i__] != 4) {
				if (sec_1.isec[i__ + 1] == 4) {
					sec_1.isec[i__] = 4;
				}
			}
		}
	}
	i__1 = length_1.nseq2;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (sec_1.jsec[i__ - 1] == 2) {
			if (sec_1.jsec[i__] != 2) {
				if (sec_1.jsec[i__ + 1] == 2) {
					sec_1.jsec[i__] = 2;
				}
			}
		}
		if (sec_1.jsec[i__ - 1] == 4) {
			if (sec_1.jsec[i__] != 4) {
				if (sec_1.jsec[i__ + 1] == 4) {
					sec_1.jsec[i__] = 4;
				}
			}
		}
	}
	return 0;
} /* smooth_ */
/* ************************************************************ */
/*     assign secondary structure: */
/* ************************************************************ */
double diszy_(int *i__, int *i1, int *i2, backbone_ *backbone_1) {
	/* System generated locals */
	float ret_val, r__1, r__2, r__3;
	/* Builtin functions */
	double sqrt(double);
	/* Computing 2nd power */
	r__1 = backbone_1->xa[(*i1 + *i__ * 5000) * 3 - 3]
			- backbone_1->xa[(*i2 + *i__ * 5000) * 3 - 3];
	/* Computing 2nd power */
	r__2 = backbone_1->xa[(*i1 + *i__ * 5000) * 3 - 2]
			- backbone_1->xa[(*i2 + *i__ * 5000) * 3 - 2];
	/* Computing 2nd power */
	r__3 = backbone_1->xa[(*i1 + *i__ * 5000) * 3 - 1]
			- backbone_1->xa[(*i2 + *i__ * 5000) * 3 - 1];
	ret_val = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
	return ret_val;
} /* diszy_ */
/* ************************************************************ */
/*     assign secondary structure: */
/* ************************************************************ */
int make_sec__(float *dis13, float *dis14, float *dis15, float *dis24,
		float * dis25, float *dis35) {
	/* System generated locals */
	int ret_val;
	float r__1;
	/* Local variables */
	float delta;
	ret_val = 1;
	delta = 2.1f;
	if ((r__1 = *dis15 - 6.37f, dabs(r__1)) < delta) {
		if ((r__1 = *dis14 - 5.18f, dabs(r__1)) < delta) {
			if ((r__1 = *dis25 - 5.18f, dabs(r__1)) < delta) {
				if ((r__1 = *dis13 - 5.45f, dabs(r__1)) < delta) {
					if ((r__1 = *dis24 - 5.45f, dabs(r__1)) < delta) {
						if ((r__1 = *dis35 - 5.45f, dabs(r__1)) < delta) {
							ret_val = 2;
							/* helix */
							return ret_val;
						}
					}
				}
			}
		}
	}
	delta = 1.42f;
	if ((r__1 = *dis15 - 13, dabs(r__1)) < delta) {
		if ((r__1 = *dis14 - 10.4f, dabs(r__1)) < delta) {
			if ((r__1 = *dis25 - 10.4f, dabs(r__1)) < delta) {
				if ((r__1 = *dis13 - 6.1f, dabs(r__1)) < delta) {
					if ((r__1 = *dis24 - 6.1f, dabs(r__1)) < delta) {
						if ((r__1 = *dis35 - 6.1f, dabs(r__1)) < delta) {
							ret_val = 4;
							/* strand */
							return ret_val;
						}
					}
				}
			}
		}
	}
	if (*dis15 < 8.f) {
		ret_val = 3;
	}
	return ret_val;
} /* make_sec__ */
/* *************************************************************** */
/*     quickly calculate TM-score with given invmap(i) in 3 iterations */
/* *************************************************************** */
/* Subroutine */int get_gl__(float *gl, dpc_ *dpc_1, backbone_ *backbone_1) {
	/* Initialized data */
	/* System generated locals */
	int i__1;
	float r__1, r__2, r__3;
	/* Local variables */
	int i__, j, l;
	double t[3], u[9] /* was [3][3] */;
	float g2, g3, xx, yy, zz;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	float xo1[5000], yo1[5000], zo1[5000], xo2[5000], yo2[5000],
			zo2[5000];
	int ier;
	double rms;
	float d002t, dis2[5000];
	/* cc   RMSD: */
	/* cc */
	/*     calculate RMSD between aligned structures and rotate the structures --> */
	/* armsd is float */zscore_1.n_al__ = 0;
	i__1 = length_1.nseq2;
	for (j = 1; j <= i__1; ++j) {
		i__ = dpc_1->invmap[j - 1];
		/* j aligned to i */
		if (i__ > 0) {
			++zscore_1.n_al__;
			r_1__[zscore_1.n_al__ * 3 - 3] = backbone_1->xa[i__ * 3 - 3];
			r_1__[zscore_1.n_al__ * 3 - 2] = backbone_1->xa[i__ * 3 - 2];
			r_1__[zscore_1.n_al__ * 3 - 1] = backbone_1->xa[i__ * 3 - 1];
			r_2__[zscore_1.n_al__ * 3 - 3] = backbone_1->xa[(j + 5000) * 3 - 3];
			r_2__[zscore_1.n_al__ * 3 - 2] = backbone_1->xa[(j + 5000) * 3 - 2];
			r_2__[zscore_1.n_al__ * 3 - 1] = backbone_1->xa[(j + 5000) * 3 - 1];
			xo1[zscore_1.n_al__ - 1] = backbone_1->xa[i__ * 3 - 3];
			yo1[zscore_1.n_al__ - 1] = backbone_1->xa[i__ * 3 - 2];
			zo1[zscore_1.n_al__ - 1] = backbone_1->xa[i__ * 3 - 1];
			xo2[zscore_1.n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 3];
			yo2[zscore_1.n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 2];
			zo2[zscore_1.n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 1];
		}
	}
	u3b_(w, r_1__, r_2__, &zscore_1.n_al__, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	*gl = 0.f;
	i__1 = zscore_1.n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xx = t[0] + u[0] * xo1[i__ - 1] + u[3] * yo1[i__ - 1]
				+ u[6] * zo1[i__ - 1];
		yy = t[1] + u[1] * xo1[i__ - 1] + u[4] * yo1[i__ - 1]
				+ u[7] * zo1[i__ - 1];
		zz = t[2] + u[2] * xo1[i__ - 1] + u[5] * yo1[i__ - 1]
				+ u[8] * zo1[i__ - 1];
		/* Computing 2nd power */
		r__1 = xx - xo2[i__ - 1];
		/* Computing 2nd power */
		r__2 = yy - yo2[i__ - 1];
		/* Computing 2nd power */
		r__3 = zz - zo2[i__ - 1];
		dis2[i__ - 1] = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
		/* Computing 2nd power */
		r__1 = d0_1.d0;
		*gl += 1 / (dis2[i__ - 1] / (r__1 * r__1) + 1);
	}
	/* cc   for next iteration-------------> */
	d002t = d00_1.d002;
	L21: j = 0;
	i__1 = zscore_1.n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (dis2[i__ - 1] <= d002t) {
			++j;
			r_1__[j * 3 - 3] = xo1[i__ - 1];
			r_1__[j * 3 - 2] = yo1[i__ - 1];
			r_1__[j * 3 - 1] = zo1[i__ - 1];
			r_2__[j * 3 - 3] = xo2[i__ - 1];
			r_2__[j * 3 - 2] = yo2[i__ - 1];
			r_2__[j * 3 - 1] = zo2[i__ - 1];
		}
	}
	if (j < 3 && zscore_1.n_al__ > 3) {
		d002t += .5f;
		goto L21;
	}
	l = j;
	u3b_(w, r_1__, r_2__, &l, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	g2 = 0.f;
	i__1 = zscore_1.n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xx = t[0] + u[0] * xo1[i__ - 1] + u[3] * yo1[i__ - 1]
				+ u[6] * zo1[i__ - 1];
		yy = t[1] + u[1] * xo1[i__ - 1] + u[4] * yo1[i__ - 1]
				+ u[7] * zo1[i__ - 1];
		zz = t[2] + u[2] * xo1[i__ - 1] + u[5] * yo1[i__ - 1]
				+ u[8] * zo1[i__ - 1];
		/* Computing 2nd power */
		r__1 = xx - xo2[i__ - 1];
		/* Computing 2nd power */
		r__2 = yy - yo2[i__ - 1];
		/* Computing 2nd power */
		r__3 = zz - zo2[i__ - 1];
		dis2[i__ - 1] = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
		/* Computing 2nd power */
		r__1 = d0_1.d0;
		g2 += 1 / (dis2[i__ - 1] / (r__1 * r__1) + 1);
	}
	/* cc   for next iteration-------------> */
	d002t = d00_1.d002 + 1;
	L22: j = 0;
	i__1 = zscore_1.n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		if (dis2[i__ - 1] <= d002t) {
			++j;
			r_1__[j * 3 - 3] = xo1[i__ - 1];
			r_1__[j * 3 - 2] = yo1[i__ - 1];
			r_1__[j * 3 - 1] = zo1[i__ - 1];
			r_2__[j * 3 - 3] = xo2[i__ - 1];
			r_2__[j * 3 - 2] = yo2[i__ - 1];
			r_2__[j * 3 - 1] = zo2[i__ - 1];
		}
	}
	if (j < 3 && zscore_1.n_al__ > 3) {
		d002t += .5f;
		goto L22;
	}
	l = j;
	u3b_(w, r_1__, r_2__, &l, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	g3 = 0.f;
	i__1 = zscore_1.n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xx = t[0] + u[0] * xo1[i__ - 1] + u[3] * yo1[i__ - 1]
				+ u[6] * zo1[i__ - 1];
		yy = t[1] + u[1] * xo1[i__ - 1] + u[4] * yo1[i__ - 1]
				+ u[7] * zo1[i__ - 1];
		zz = t[2] + u[2] * xo1[i__ - 1] + u[5] * yo1[i__ - 1]
				+ u[8] * zo1[i__ - 1];
		/* Computing 2nd power */
		r__1 = xx - xo2[i__ - 1];
		/* Computing 2nd power */
		r__2 = yy - yo2[i__ - 1];
		/* Computing 2nd power */
		r__3 = zz - zo2[i__ - 1];
		dis2[i__ - 1] = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
		/* Computing 2nd power */
		r__1 = d0_1.d0;
		g3 += 1 / (dis2[i__ - 1] / (r__1 * r__1) + 1);
	}
	if (g2 > *gl) {
		*gl = g2;
	}
	if (g3 > *gl) {
		*gl = g3;
	}
	/* ^^^^^^^^^^^^^^^^ GL done ^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* get_gl__ */
/* *************************************************************** */
/*     with invmap(i) calculate TM-score and martix score(i,j) for rotation */
/* *************************************************************** */
/* Subroutine */int get_score__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* Initialized data */
	/* System generated locals */
	int i__1, i__2;
	float r__1, r__2, r__3;
	/* Local variables */
	float d0_input__;
	int i__, j;
	double t[3], u[9] /* was [3][3] */;
	float dd, xx, yy, zz;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	int ier;
	double rms;
	float xtm1[5000], ytm1[5000], ztm1[5000], xtm2[5000], ytm2[5000],
			ztm2[5000];
	int lcomm;
	float rcomm;
	/* cc   RMSD: */
	/* cc */
	/*     calculate RMSD between aligned structures and rotate the structures --> */
	/* armsd is float */zscore_1.n_al__ = 0;
	i__1 = length_1.nseq2;
	for (j = 1; j <= i__1; ++j) {
		i__ = dpc_1->invmap[j - 1];
		/* j aligned to i */
		if (i__ > 0) {
			++zscore_1.n_al__;
			/* cc   for TM-score: */
			xtm1[zscore_1.n_al__ - 1] = backbone_1->xa[i__ * 3 - 3];
			/* for TM-score */
			ytm1[zscore_1.n_al__ - 1] = backbone_1->xa[i__ * 3 - 2];
			ztm1[zscore_1.n_al__ - 1] = backbone_1->xa[i__ * 3 - 1];
			xtm2[zscore_1.n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 3];
			ytm2[zscore_1.n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 2];
			ztm2[zscore_1.n_al__ - 1] = backbone_1->xa[(j + 5000) * 3 - 1];
			/* cc   for rotation matrix: */
			r_1__[zscore_1.n_al__ * 3 - 3] = backbone_1->xa[i__ * 3 - 3];
			r_1__[zscore_1.n_al__ * 3 - 2] = backbone_1->xa[i__ * 3 - 2];
			r_1__[zscore_1.n_al__ * 3 - 1] = backbone_1->xa[i__ * 3 - 1];
		}
	}
	/* **   calculate TM-score for the given alignment-----------> */
	d0_input__ = d0_1.d0;
	tmscore8_search__(&d0_input__, &zscore_1.n_al__, xtm1, ytm1, ztm1,
			&zscore_1.n_al__, xtm2, ytm2, ztm2, &tm_1.tm, &rcomm, &lcomm);
	/* simplified search engine */tm_1.tm = tm_1.tm * zscore_1.n_al__
			/ d0_1.anseq;
	/* **   calculate score matrix score(i,j)------------------> */
	/* TM-score */
	i__1 = zscore_1.n_al__;
	for (i__ = 1; i__ <= i__1; ++i__) {
		r_2__[i__ * 3 - 3] = xtm1[i__ - 1];
		r_2__[i__ * 3 - 2] = ytm1[i__ - 1];
		r_2__[i__ * 3 - 1] = ztm1[i__ - 1];
	}
	u3b_(w, r_1__, r_2__, &zscore_1.n_al__, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xx = t[0] + u[0] * backbone_1->xa[i__ * 3 - 3]
				+ u[3] * backbone_1->xa[i__ * 3 - 2]
				+ u[6] * backbone_1->xa[i__ * 3 - 1];
		yy = t[1] + u[1] * backbone_1->xa[i__ * 3 - 3]
				+ u[4] * backbone_1->xa[i__ * 3 - 2]
				+ u[7] * backbone_1->xa[i__ * 3 - 1];
		zz = t[2] + u[2] * backbone_1->xa[i__ * 3 - 3]
				+ u[5] * backbone_1->xa[i__ * 3 - 2]
				+ u[8] * backbone_1->xa[i__ * 3 - 1];
		i__2 = length_1.nseq2;
		for (j = 1; j <= i__2; ++j) {
			/* Computing 2nd power */
			r__1 = xx - backbone_1->xa[(j + 5000) * 3 - 3];
			/* Computing 2nd power */
			r__2 = yy - backbone_1->xa[(j + 5000) * 3 - 2];
			/* Computing 2nd power */
			r__3 = zz - backbone_1->xa[(j + 5000) * 3 - 1];
			dd = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			/* Computing 2nd power */
			r__1 = d0_1.d0;
			dpc_1->score[i__ + j * 5000 - 5001] = 1 / (dd / (r__1 * r__1) + 1);
		}
	}
	/* ^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* get_score__ */
/* *************************************************************** */
/*     with invmap(i) calculate score(i,j) using RMSD rotation */
/* *************************************************************** */
/* Subroutine */int get_score1__(dpc_ *dpc_1, backbone_ *backbone_1) {
	/* Initialized data */
	/* System generated locals */
	int i__1, i__2;
	float r__1, r__2, r__3;
	/* Local variables */
	int i__, j;
	double t[3], u[9] /* was [3][3] */;
	float d01, d02, dd, xx, yy, zz;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	int ier;
	double rms;
	/* cc   RMSD: */
	/* cc */
	/*     calculate RMSD between aligned structures and rotate the structures --> */
	/* armsd is float */zscore_1.n_al__ = 0;
	i__1 = length_1.nseq2;
	for (j = 1; j <= i__1; ++j) {
		i__ = dpc_1->invmap[j - 1];
		/* j aligned to i */
		if (i__ > 0) {
			++zscore_1.n_al__;
			/* cc   for rotation matrix: */
			r_1__[zscore_1.n_al__ * 3 - 3] = backbone_1->xa[i__ * 3 - 3];
			r_1__[zscore_1.n_al__ * 3 - 2] = backbone_1->xa[i__ * 3 - 2];
			r_1__[zscore_1.n_al__ * 3 - 1] = backbone_1->xa[i__ * 3 - 1];
			r_2__[zscore_1.n_al__ * 3 - 3] = backbone_1->xa[(j + 5000) * 3 - 3];
			r_2__[zscore_1.n_al__ * 3 - 2] = backbone_1->xa[(j + 5000) * 3 - 2];
			r_2__[zscore_1.n_al__ * 3 - 1] = backbone_1->xa[(j + 5000) * 3 - 1];
		}
	}
	/* **   calculate score matrix score(i,j)------------------> */
	u3b_(w, r_1__, r_2__, &zscore_1.n_al__, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	d01 = d0_1.d0 + 1.5f;
	if (d01 < d0min_1.d0_min__) {
		d01 = d0min_1.d0_min__;
	}
	d02 = d01 * d01;
	i__1 = length_1.nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xx = t[0] + u[0] * backbone_1->xa[i__ * 3 - 3]
				+ u[3] * backbone_1->xa[i__ * 3 - 2]
				+ u[6] * backbone_1->xa[i__ * 3 - 1];
		yy = t[1] + u[1] * backbone_1->xa[i__ * 3 - 3]
				+ u[4] * backbone_1->xa[i__ * 3 - 2]
				+ u[7] * backbone_1->xa[i__ * 3 - 1];
		zz = t[2] + u[2] * backbone_1->xa[i__ * 3 - 3]
				+ u[5] * backbone_1->xa[i__ * 3 - 2]
				+ u[8] * backbone_1->xa[i__ * 3 - 1];
		i__2 = length_1.nseq2;
		for (j = 1; j <= i__2; ++j) {
			/* Computing 2nd power */
			r__1 = xx - backbone_1->xa[(j + 5000) * 3 - 3];
			/* Computing 2nd power */
			r__2 = yy - backbone_1->xa[(j + 5000) * 3 - 2];
			/* Computing 2nd power */
			r__3 = zz - backbone_1->xa[(j + 5000) * 3 - 1];
			dd = r__1 * r__1 + r__2 * r__2 + r__3 * r__3;
			dpc_1->score[i__ + j * 5000 - 5001] = 1 / (dd / d02 + 1);
		}
	}
	/* ^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* get_score1__ */
/* ************************************************************************ */
/* ************************************************************************ */
/*     This is a subroutine to compare two structures and find the */
/*     superposition that has the maximum TM-score. */
/*     L1--Length of the first structure */
/*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure */
/*     L2--Length of the second structure */
/*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure */
/*     TM--TM-score of the comparison */
/*     Rcomm--RMSD of two structures in the common aligned residues */
/*     Lcomm--Length of the common aligned regions */
/*     Note: */
/*     1, Always put native as the second structure, by which TM-score */
/*        is normalized. */
/*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after */
/*        TM-score superposition. */
/* ************************************************************************ */
/* ************************************************************************ */
/* ** dis<8, simplified search engine */
/* Subroutine */int tmscore8_search__(float *dx, int *l1, float *x1, float * y1,
		float *z1, int *l2, float *x2, float *y2, float *z2, float *tm,
		float *rcomm, int *lcomm) {
	/* Initialized data */
	/* System generated locals */
	int i__1, i__2, i__3, i__4;
	double sqrt(double);
	/* Local variables */
	float d0_search__;
	int i__, j, k, m;
	double t[3], u[9] /* was [3][3] */;
	int l_ini_min__;
	double score_max__;
	int ka, il, ll;
	float xa[5000], ya[5000], za[5000];
	int it, ka0;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	int n_init_max__, il0[5000];
	int ier, neq;
	double rms;
	int n_it__, k_ali__[5000], l_ini__[100];
	float armsd;
	int k_ali0__[5000], il_max__, i_init__, l_init__, n_init__,
			i_shift__, n_shift__;
	/* [1,n_ali],align residues for the */
	/* cc   RMSD: */
	/* cc */
	/* ******** convert input data *************************** */
	/*     because L1=L2 in this special case----------> */
	/* armsd is float */
	/* Parameter adjustments */
	--z2;
	--y2;
	--x2;
	--z1;
	--y1;
	--x1;
	/* Function Body */nres_1.nseqa = *l1;
	nres_1.nseqb = *l2;
	i__1 = nres_1.nseqa;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xa[i__ - 1] = x1[i__];
		ya[i__ - 1] = y1[i__];
		za[i__ - 1] = z1[i__];
		stru_1.xb[i__ - 1] = x2[i__];
		stru_1.yb[i__ - 1] = y2[i__];
		stru_1.zb[i__ - 1] = z2[i__];
		align_1.ia[i__ - 1] = i__;
		align_1.ib[i__ - 1] = i__;
	}
	align_1.n_ali__ = *l1;
	/* number of aligned residues */
	*lcomm = *l1;
	/* **********+///// */
	/*     parameters: */
	/* **************** */
	/* **   d0-------------> */para_1.d0 = *dx;
	if (para_1.d0 < d0min_1.d0_min__) {
		para_1.d0 = d0min_1.d0_min__;
	}
	/* **   d0_search -----> */
	d0_search__ = para_1.d0;
	if (d0_search__ > 8.f) {
		d0_search__ = 8.f;
	}
	if (d0_search__ < 4.5f) {
		d0_search__ = 4.5f;
	}
	/* **   iterative parameters -----> */
	n_it__ = 20;
	/* maximum number of iterations */
	n_init_max__ = 6;
	/* maximum number of L_init */
	n_init__ = 0;
	l_ini_min__ = 4;
	if (align_1.n_ali__ < 4) {
		l_ini_min__ = align_1.n_ali__;
	}
	i__1 = n_init_max__ - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		++n_init__;
		i__2 = n_init__ - 1;
		l_ini__[n_init__ - 1] = align_1.n_ali__ / pow_ii_f2c(&c__2, &i__2);
		if (l_ini__[n_init__ - 1] <= l_ini_min__) {
			l_ini__[n_init__ - 1] = l_ini_min__;
			goto L402;
		}
	}
	++n_init__;
	l_ini__[n_init__ - 1] = l_ini_min__;
	L402:
	/* ***************************************************************** */
	/*     find the maximum score starting from local structures superposition */
	/* ***************************************************************** */
	score_max__ = -1.;
	/* TM-score */
	i__1 = n_init__;
	for (i_init__ = 1; i_init__ <= i__1; ++i_init__) {
		l_init__ = l_ini__[i_init__ - 1];
		il_max__ = align_1.n_ali__ - l_init__ + 1;
		k = 0;
		i__2 = il_max__;
		for (i__ = 1; i__ <= i__2; i__ += 40) {
			/* this is the simplification! */
			++k;
			il0[k - 1] = i__;
		}
		if (il0[k - 1] < il_max__) {
			++k;
			il0[k - 1] = il_max__;
		}
		n_shift__ = k;
		i__2 = n_shift__;
		for (i_shift__ = 1; i_shift__ <= i__2; ++i_shift__) {
			il = il0[i_shift__ - 1];
			ll = 0;
			ka = 0;
			i__3 = l_init__;
			for (i__ = 1; i__ <= i__3; ++i__) {
				k = il + i__ - 1;
				/* [1,n_ali] common aligned */
				r_1__[i__ * 3 - 3] = xa[align_1.ia[k - 1] - 1];
				r_1__[i__ * 3 - 2] = ya[align_1.ia[k - 1] - 1];
				r_1__[i__ * 3 - 1] = za[align_1.ia[k - 1] - 1];
				r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[k - 1] - 1];
				r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[k - 1] - 1];
				r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[k - 1] - 1];
				++ll;
				++ka;
				k_ali__[ka - 1] = k;
			}
			u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
			/* u rotate r_1 to r_2 */
			if (i_init__ == 1) {
				/* global superposition */
				armsd = sqrt(rms / ll);
				*rcomm = armsd;
			}
			i__3 = nres_1.nseqa;
			for (j = 1; j <= i__3; ++j) {
				stru_1.xt[j - 1] = t[0] + u[0] * xa[j - 1] + u[3] * ya[j - 1]
						+ u[6] * za[j - 1];
				stru_1.yt[j - 1] = t[1] + u[1] * xa[j - 1] + u[4] * ya[j - 1]
						+ u[7] * za[j - 1];
				stru_1.zt[j - 1] = t[2] + u[2] * xa[j - 1] + u[5] * ya[j - 1]
						+ u[8] * za[j - 1];
			}
			para_1.d__ = d0_search__ - 1;
			score_fun8__();
			/* init, get scores, n_cut+i_ali(i) for i */
			if (score_max__ < scores_1.score) {
				score_max__ = scores_1.score;
				ka0 = ka;
				i__3 = ka0;
				for (i__ = 1; i__ <= i__3; ++i__) {
					k_ali0__[i__ - 1] = k_ali__[i__ - 1];
				}
			}
			/* **   iteration for extending ----------------------------------> */para_1.d__ =
					d0_search__ + 1;
			i__3 = n_it__;
			for (it = 1; it <= i__3; ++it) {
				ll = 0;
				ka = 0;
				i__4 = nscore_1.n_cut__;
				for (i__ = 1; i__ <= i__4; ++i__) {
					m = nscore_1.i_ali__[i__ - 1];
					/* [1,n_ali] */
					r_1__[i__ * 3 - 3] = xa[align_1.ia[m - 1] - 1];
					r_1__[i__ * 3 - 2] = ya[align_1.ia[m - 1] - 1];
					r_1__[i__ * 3 - 1] = za[align_1.ia[m - 1] - 1];
					r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[m - 1] - 1];
					r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[m - 1] - 1];
					r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[m - 1] - 1];
					++ka;
					k_ali__[ka - 1] = m;
					++ll;
				}
				u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
				/* u rotate r_1 to r_2 */
				i__4 = nres_1.nseqa;
				for (j = 1; j <= i__4; ++j) {
					stru_1.xt[j - 1] = t[0] + u[0] * xa[j - 1]
							+ u[3] * ya[j - 1] + u[6] * za[j - 1];
					stru_1.yt[j - 1] = t[1] + u[1] * xa[j - 1]
							+ u[4] * ya[j - 1] + u[7] * za[j - 1];
					stru_1.zt[j - 1] = t[2] + u[2] * xa[j - 1]
							+ u[5] * ya[j - 1] + u[8] * za[j - 1];
				}
				score_fun8__();
				/* get scores, n_cut+i_ali(i) for iterati */
				if (score_max__ < scores_1.score) {
					score_max__ = scores_1.score;
					ka0 = ka;
					i__4 = ka;
					for (i__ = 1; i__ <= i__4; ++i__) {
						k_ali0__[i__ - 1] = k_ali__[i__ - 1];
					}
				}
				if (it == n_it__) {
					goto L302;
				}
				if (nscore_1.n_cut__ == ka) {
					neq = 0;
					i__4 = nscore_1.n_cut__;
					for (i__ = 1; i__ <= i__4; ++i__) {
						if (nscore_1.i_ali__[i__ - 1] == k_ali__[i__ - 1]) {
							++neq;
						}
					}
					if (nscore_1.n_cut__ == neq) {
						goto L302;
					}
				}
				/* L301: */
			}
			/* for iteration */
			L302:
			/* L300: */
			;
		}
		/* for shift */
		/* L333: */
	}
	/* ******* return the final rotation **************** */
	/* for initial length, L_ali/M */
	ll = 0;
	i__1 = ka0;
	for (i__ = 1; i__ <= i__1; ++i__) {
		m = k_ali0__[i__ - 1];
		/* record of the best alignment */
		r_1__[i__ * 3 - 3] = xa[align_1.ia[m - 1] - 1];
		r_1__[i__ * 3 - 2] = ya[align_1.ia[m - 1] - 1];
		r_1__[i__ * 3 - 1] = za[align_1.ia[m - 1] - 1];
		r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[m - 1] - 1];
		r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[m - 1] - 1];
		r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[m - 1] - 1];
		++ll;
	}
	u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	i__1 = nres_1.nseqa;
	for (j = 1; j <= i__1; ++j) {
		x1[j] = t[0] + u[0] * xa[j - 1] + u[3] * ya[j - 1] + u[6] * za[j - 1];
		y1[j] = t[1] + u[1] * xa[j - 1] + u[4] * ya[j - 1] + u[7] * za[j - 1];
		z1[j] = t[2] + u[2] * xa[j - 1] + u[5] * ya[j - 1] + u[8] * za[j - 1];
	}
	*tm = score_max__;
	/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* tmscore8_search__ */
/* ************************************************************************ */
/* ************************************************************************ */
/*     This is a subroutine to compare two structures and find the */
/*     superposition that has the maximum TM-score. */
/*     L1--Length of the first structure */
/*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure */
/*     L2--Length of the second structure */
/*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure */
/*     TM--TM-score of the comparison */
/*     Rcomm--RMSD of two structures in the common aligned residues */
/*     Lcomm--Length of the common aligned regions */
/*     Note: */
/*     1, Always put native as the second structure, by which TM-score */
/*        is normalized. */
/*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after */
/*        TM-score superposition. */
/* ************************************************************************ */
/* ************************************************************************ */
/* **   dis<8, but same search engine */
/* Subroutine */int tmscore8_(float *dx, int *l1, float *x1, float *y1,
		float *z1, int *l2, float *x2, float *y2, float *z2, float *tm,
		float * rcomm, int *lcomm) {
	/* Initialized data */
	/* System generated locals */
	int i__1, i__2, i__3, i__4;
	/* Builtin functions */
	int pow_ii(int *, int *);
	double sqrt(double);
	/* Local variables */
	float d0_search__;
	int i__, j, k, m;
	double t[3], u[9] /* was [3][3] */;
	int l_ini_min__;
	double score_max__;
	int ka, il, ll;
	float xa[5000], ya[5000], za[5000];
	int it, ka0;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	int n_init_max__;
	int ier, neq;
	double rms;
	int n_it__, k_ali__[5000], l_ini__[100];
	float armsd;
	int k_ali0__[5000], il_max__, i_init__, l_init__, n_init__;
	/* [1,n_ali],align residues for the */
	/* cc   RMSD: */
	/* cc */
	/* ******** convert input data *************************** */
	/*     because L1=L2 in this special case----------> */
	/* armsd is float */
	/* Parameter adjustments */
	--z2;
	--y2;
	--x2;
	--z1;
	--y1;
	--x1;
	/* Function Body */nres_1.nseqa = *l1;
	nres_1.nseqb = *l2;
	i__1 = nres_1.nseqa;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xa[i__ - 1] = x1[i__];
		ya[i__ - 1] = y1[i__];
		za[i__ - 1] = z1[i__];
		stru_1.xb[i__ - 1] = x2[i__];
		stru_1.yb[i__ - 1] = y2[i__];
		stru_1.zb[i__ - 1] = z2[i__];
		align_1.ia[i__ - 1] = i__;
		align_1.ib[i__ - 1] = i__;
	}
	align_1.n_ali__ = *l1;
	/* number of aligned residues */
	*lcomm = *l1;
	/* **********+///// */
	/*     parameters: */
	/* **************** */
	/* **   d0-------------> */para_1.d0 = *dx;
	if (para_1.d0 < d0min_1.d0_min__) {
		para_1.d0 = d0min_1.d0_min__;
	}
	/* **   d0_search -----> */
	d0_search__ = para_1.d0;
	if (d0_search__ > 8.f) {
		d0_search__ = 8.f;
	}
	if (d0_search__ < 4.5f) {
		d0_search__ = 4.5f;
	}
	/* **   iterative parameters -----> */
	n_it__ = 20;
	/* maximum number of iterations */
	n_init_max__ = 6;
	/* maximum number of L_init */
	n_init__ = 0;
	l_ini_min__ = 4;
	if (align_1.n_ali__ < 4) {
		l_ini_min__ = align_1.n_ali__;
	}
	i__1 = n_init_max__ - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		++n_init__;
		i__2 = n_init__ - 1;
		l_ini__[n_init__ - 1] = align_1.n_ali__ / pow_ii_f2c(&c__2, &i__2);
		if (l_ini__[n_init__ - 1] <= l_ini_min__) {
			l_ini__[n_init__ - 1] = l_ini_min__;
			goto L402;
		}
	}
	++n_init__;
	l_ini__[n_init__ - 1] = l_ini_min__;
	L402:
	/* ***************************************************************** */
	/*     find the maximum score starting from local structures superposition */
	/* ***************************************************************** */
	score_max__ = -1.;
	/* TM-score */
	i__1 = n_init__;
	for (i_init__ = 1; i_init__ <= i__1; ++i_init__) {
		l_init__ = l_ini__[i_init__ - 1];
		il_max__ = align_1.n_ali__ - l_init__ + 1;
		i__2 = il_max__;
		for (il = 1; il <= i__2; ++il) {
			/* on aligned residues, [1,nseqA] */
			ll = 0;
			ka = 0;
			i__3 = l_init__;
			for (i__ = 1; i__ <= i__3; ++i__) {
				k = il + i__ - 1;
				/* [1,n_ali] common aligned */
				r_1__[i__ * 3 - 3] = xa[align_1.ia[k - 1] - 1];
				r_1__[i__ * 3 - 2] = ya[align_1.ia[k - 1] - 1];
				r_1__[i__ * 3 - 1] = za[align_1.ia[k - 1] - 1];
				r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[k - 1] - 1];
				r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[k - 1] - 1];
				r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[k - 1] - 1];
				++ll;
				++ka;
				k_ali__[ka - 1] = k;
			}
			u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
			/* u rotate r_1 to r_2 */
			if (i_init__ == 1) {
				/* global superposition */
				armsd = sqrt(rms / ll);
				*rcomm = armsd;
			}
			i__3 = nres_1.nseqa;
			for (j = 1; j <= i__3; ++j) {
				stru_1.xt[j - 1] = t[0] + u[0] * xa[j - 1] + u[3] * ya[j - 1]
						+ u[6] * za[j - 1];
				stru_1.yt[j - 1] = t[1] + u[1] * xa[j - 1] + u[4] * ya[j - 1]
						+ u[7] * za[j - 1];
				stru_1.zt[j - 1] = t[2] + u[2] * xa[j - 1] + u[5] * ya[j - 1]
						+ u[8] * za[j - 1];
			}
			para_1.d__ = d0_search__ - 1;
			score_fun8__();
			/* init, get scores, n_cut+i_ali(i) for i */
			if (score_max__ < scores_1.score) {
				score_max__ = scores_1.score;
				ka0 = ka;
				i__3 = ka0;
				for (i__ = 1; i__ <= i__3; ++i__) {
					k_ali0__[i__ - 1] = k_ali__[i__ - 1];
				}
			}
			/* **   iteration for extending ----------------------------------> */para_1.d__ =
					d0_search__ + 1;
			i__3 = n_it__;
			for (it = 1; it <= i__3; ++it) {
				ll = 0;
				ka = 0;
				i__4 = nscore_1.n_cut__;
				for (i__ = 1; i__ <= i__4; ++i__) {
					m = nscore_1.i_ali__[i__ - 1];
					/* [1,n_ali] */
					r_1__[i__ * 3 - 3] = xa[align_1.ia[m - 1] - 1];
					r_1__[i__ * 3 - 2] = ya[align_1.ia[m - 1] - 1];
					r_1__[i__ * 3 - 1] = za[align_1.ia[m - 1] - 1];
					r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[m - 1] - 1];
					r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[m - 1] - 1];
					r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[m - 1] - 1];
					++ka;
					k_ali__[ka - 1] = m;
					++ll;
				}
				u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
				/* u rotate r_1 to r_2 */
				i__4 = nres_1.nseqa;
				for (j = 1; j <= i__4; ++j) {
					stru_1.xt[j - 1] = t[0] + u[0] * xa[j - 1]
							+ u[3] * ya[j - 1] + u[6] * za[j - 1];
					stru_1.yt[j - 1] = t[1] + u[1] * xa[j - 1]
							+ u[4] * ya[j - 1] + u[7] * za[j - 1];
					stru_1.zt[j - 1] = t[2] + u[2] * xa[j - 1]
							+ u[5] * ya[j - 1] + u[8] * za[j - 1];
				}
				score_fun8__();
				/* get scores, n_cut+i_ali(i) for iterati */
				if (score_max__ < scores_1.score) {
					score_max__ = scores_1.score;
					ka0 = ka;
					i__4 = ka;
					for (i__ = 1; i__ <= i__4; ++i__) {
						k_ali0__[i__ - 1] = k_ali__[i__ - 1];
					}
				}
				if (it == n_it__) {
					goto L302;
				}
				if (nscore_1.n_cut__ == ka) {
					neq = 0;
					i__4 = nscore_1.n_cut__;
					for (i__ = 1; i__ <= i__4; ++i__) {
						if (nscore_1.i_ali__[i__ - 1] == k_ali__[i__ - 1]) {
							++neq;
						}
					}
					if (nscore_1.n_cut__ == neq) {
						goto L302;
					}
				}
				/* L301: */
			}
			/* for iteration */
			L302:
			/* L300: */
			;
		}
		/* for shift */
		/* L333: */
	}
	/* ******* return the final rotation **************** */
	/* for initial length, L_ali/M */
	ll = 0;
	i__1 = ka0;
	for (i__ = 1; i__ <= i__1; ++i__) {
		m = k_ali0__[i__ - 1];
		/* record of the best alignment */
		r_1__[i__ * 3 - 3] = xa[align_1.ia[m - 1] - 1];
		r_1__[i__ * 3 - 2] = ya[align_1.ia[m - 1] - 1];
		r_1__[i__ * 3 - 1] = za[align_1.ia[m - 1] - 1];
		r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[m - 1] - 1];
		r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[m - 1] - 1];
		r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[m - 1] - 1];
		++ll;
	}
	u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	i__1 = nres_1.nseqa;
	for (j = 1; j <= i__1; ++j) {
		x1[j] = t[0] + u[0] * xa[j - 1] + u[3] * ya[j - 1] + u[6] * za[j - 1];
		y1[j] = t[1] + u[1] * xa[j - 1] + u[4] * ya[j - 1] + u[7] * za[j - 1];
		z1[j] = t[2] + u[2] * xa[j - 1] + u[5] * ya[j - 1] + u[8] * za[j - 1];
	}
	*tm = score_max__;
	/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* tmscore8_ */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     1, collect those residues with dis<d; */
/*     2, calculate score_GDT, score_maxsub, score_TM */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */int score_fun8__(void) {
	/* System generated locals */
	int i__1;
	float r__1, r__2, r__3;
	/* Builtin functions */
	double sqrt(double);
	/* Local variables */
	int i__, j, k;
	float score_sum__, dis, d_tmp__;
	/* [1,n_ali],align residues for the */
	d_tmp__ = para_1.d__;
	L21:
	nscore_1.n_cut__ = 0;
	/* number of residue-pairs dis<d, for iter */
	score_sum__ = 0.f;
	/* TMscore */
	i__1 = align_1.n_ali__;
	for (k = 1; k <= i__1; ++k) {
		i__ = align_1.ia[k - 1];
		/* [1,nseqA] reoder number of structureA */
		j = align_1.ib[k - 1];
		/* [1,nseqB] */
		/* Computing 2nd power */
		r__1 = stru_2.xa[i__ - 1] - stru_2.xb[j - 1];
		/* Computing 2nd power */
		r__2 = stru_2.ya[i__ - 1] - stru_2.yb[j - 1];
		/* Computing 2nd power */
		r__3 = stru_2.za[i__ - 1] - stru_2.zb[j - 1];
		dis = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		if (dis < d_tmp__) {
			++nscore_1.n_cut__;
			nscore_1.i_ali__[nscore_1.n_cut__ - 1] = k;
			/* [1,n_ali], mark the residue-pairs in di */
		}
		if (dis <= d8_1.d8) {
			/* Computing 2nd power */
			r__1 = dis / para_1.d0;
			score_sum__ += 1 / (r__1 * r__1 + 1);
		}
	}
	if (nscore_1.n_cut__ < 3 && align_1.n_ali__ > 3) {
		d_tmp__ += .5f;
		goto L21;
	}
	scores_1.score = score_sum__ / (float) nres_1.nseqb;
	/* TM-score */
	return 0;
} /* score_fun8__ */
/* ************************************************************************ */
/* ************************************************************************ */
/*     This is a subroutine to compare two structures and find the */
/*     superposition that has the maximum TM-score. */
/*     L1--Length of the first structure */
/*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure */
/*     L2--Length of the second structure */
/*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure */
/*     TM--TM-score of the comparison */
/*     Rcomm--RMSD of two structures in the common aligned residues */
/*     Lcomm--Length of the common aligned regions */
/*     Note: */
/*     1, Always put native as the second structure, by which TM-score */
/*        is normalized. */
/*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after */
/*        TM-score superposition. */
/* ************************************************************************ */
/* ************************************************************************ */
/* **  normal TM-score: */
/* Subroutine */int tmscore_(float *dx, int *l1, float *x1, float *y1,
		float *z1, int *l2, float *x2, float *y2, float *z2, float *tm,
		float *rcomm, int *lcomm) {
	/* Initialized data */
	/* System generated locals */
	int i__1, i__2, i__3, i__4;
	double sqrt(double);
	/* Local variables */
	float d0_search__;
	int i__, j, k, m;
	double t[3], u[9] /* was [3][3] */;
	int l_ini_min__;
	double score_max__;
	int ka, il, ll;
	float xa[5000], ya[5000], za[5000];
	int it, ka0;
	double r_1__[15000] /* was [3][5000] */, r_2__[15000]
	/* was [3][5000] */;
	int n_init_max__;
	int ier, neq;
	double rms;
	int n_it__, k_ali__[5000], l_ini__[100];
	float armsd;
	int k_ali0__[5000], il_max__, i_init__, l_init__, n_init__;
	/* [1,n_ali],align residues for the */
	/* cc   RMSD: */
	/* cc */
	/* ******** convert input data *************************** */
	/*     because L1=L2 in this special case----------> */
	/* armsd is float */
	/* Parameter adjustments */
	--z2;
	--y2;
	--x2;
	--z1;
	--y1;
	--x1;
	/* Function Body */nres_1.nseqa = *l1;
	nres_1.nseqb = *l2;
	i__1 = nres_1.nseqa;
	for (i__ = 1; i__ <= i__1; ++i__) {
		xa[i__ - 1] = x1[i__];
		ya[i__ - 1] = y1[i__];
		za[i__ - 1] = z1[i__];
		stru_1.xb[i__ - 1] = x2[i__];
		stru_1.yb[i__ - 1] = y2[i__];
		stru_1.zb[i__ - 1] = z2[i__];
		align_1.ia[i__ - 1] = i__;
		align_1.ib[i__ - 1] = i__;
	}
	align_1.n_ali__ = *l1;
	/* number of aligned residues */
	*lcomm = *l1;
	/* **********+///// */
	/*     parameters: */
	/* **************** */
	/* **   d0-------------> */
	/*      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8 */para_1.d0 = *dx;
	if (para_1.d0 < d0min_1.d0_min__) {
		para_1.d0 = d0min_1.d0_min__;
	}
	/* **   d0_search -----> */
	d0_search__ = para_1.d0;
	if (d0_search__ > 8.f) {
		d0_search__ = 8.f;
	}
	if (d0_search__ < 4.5f) {
		d0_search__ = 4.5f;
	}
	/* **   iterative parameters -----> */
	n_it__ = 20;
	/* maximum number of iterations */
	n_init_max__ = 6;
	/* maximum number of L_init */
	n_init__ = 0;
	l_ini_min__ = 4;
	if (align_1.n_ali__ < 4) {
		l_ini_min__ = align_1.n_ali__;
	}
	i__1 = n_init_max__ - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		++n_init__;
		i__2 = n_init__ - 1;
		l_ini__[n_init__ - 1] = align_1.n_ali__ / pow_ii_f2c(&c__2, &i__2);
		if (l_ini__[n_init__ - 1] <= l_ini_min__) {
			l_ini__[n_init__ - 1] = l_ini_min__;
			goto L402;
		}
	}
	++n_init__;
	l_ini__[n_init__ - 1] = l_ini_min__;
	L402:
	/* ***************************************************************** */
	/*     find the maximum score starting from local structures superposition */
	/* ***************************************************************** */
	score_max__ = -1.;
	/* TM-score */
	i__1 = n_init__;
	for (i_init__ = 1; i_init__ <= i__1; ++i_init__) {
		l_init__ = l_ini__[i_init__ - 1];
		il_max__ = align_1.n_ali__ - l_init__ + 1;
		i__2 = il_max__;
		for (il = 1; il <= i__2; ++il) {
			/* on aligned residues, [1,nseqA] */
			ll = 0;
			ka = 0;
			i__3 = l_init__;
			for (i__ = 1; i__ <= i__3; ++i__) {
				k = il + i__ - 1;
				/* [1,n_ali] common aligned */
				r_1__[i__ * 3 - 3] = xa[align_1.ia[k - 1] - 1];
				r_1__[i__ * 3 - 2] = ya[align_1.ia[k - 1] - 1];
				r_1__[i__ * 3 - 1] = za[align_1.ia[k - 1] - 1];
				r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[k - 1] - 1];
				r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[k - 1] - 1];
				r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[k - 1] - 1];
				++ll;
				++ka;
				k_ali__[ka - 1] = k;
			}
			u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
			/* u rotate r_1 to r_2 */
			if (i_init__ == 1) {
				/* global superposition */
				armsd = sqrt(rms / ll);
				*rcomm = armsd;
			}
			i__3 = nres_1.nseqa;
			for (j = 1; j <= i__3; ++j) {
				stru_1.xt[j - 1] = t[0] + u[0] * xa[j - 1] + u[3] * ya[j - 1]
						+ u[6] * za[j - 1];
				stru_1.yt[j - 1] = t[1] + u[1] * xa[j - 1] + u[4] * ya[j - 1]
						+ u[7] * za[j - 1];
				stru_1.zt[j - 1] = t[2] + u[2] * xa[j - 1] + u[5] * ya[j - 1]
						+ u[8] * za[j - 1];
			}
			para_1.d__ = d0_search__ - 1;
			score_fun__();
			/* init, get scores, n_cut+i_ali(i) for it */
			if (score_max__ < scores_1.score) {
				score_max__ = scores_1.score;
				ka0 = ka;
				i__3 = ka0;
				for (i__ = 1; i__ <= i__3; ++i__) {
					k_ali0__[i__ - 1] = k_ali__[i__ - 1];
				}
			}
			/* **   iteration for extending ----------------------------------> */para_1.d__ =
					d0_search__ + 1;
			i__3 = n_it__;
			for (it = 1; it <= i__3; ++it) {
				ll = 0;
				ka = 0;
				i__4 = nscore_1.n_cut__;
				for (i__ = 1; i__ <= i__4; ++i__) {
					m = nscore_1.i_ali__[i__ - 1];
					/* [1,n_ali] */
					r_1__[i__ * 3 - 3] = xa[align_1.ia[m - 1] - 1];
					r_1__[i__ * 3 - 2] = ya[align_1.ia[m - 1] - 1];
					r_1__[i__ * 3 - 1] = za[align_1.ia[m - 1] - 1];
					r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[m - 1] - 1];
					r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[m - 1] - 1];
					r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[m - 1] - 1];
					++ka;
					k_ali__[ka - 1] = m;
					++ll;
				}
				u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
				/* u rotate r_1 to r_2 */
				i__4 = nres_1.nseqa;
				for (j = 1; j <= i__4; ++j) {
					stru_1.xt[j - 1] = t[0] + u[0] * xa[j - 1]
							+ u[3] * ya[j - 1] + u[6] * za[j - 1];
					stru_1.yt[j - 1] = t[1] + u[1] * xa[j - 1]
							+ u[4] * ya[j - 1] + u[7] * za[j - 1];
					stru_1.zt[j - 1] = t[2] + u[2] * xa[j - 1]
							+ u[5] * ya[j - 1] + u[8] * za[j - 1];
				}
				score_fun__();
				/* get scores, n_cut+i_ali(i) for iteratio */
				if (score_max__ < scores_1.score) {
					score_max__ = scores_1.score;
					ka0 = ka;
					i__4 = ka;
					for (i__ = 1; i__ <= i__4; ++i__) {
						k_ali0__[i__ - 1] = k_ali__[i__ - 1];
					}
				}
				if (it == n_it__) {
					goto L302;
				}
				if (nscore_1.n_cut__ == ka) {
					neq = 0;
					i__4 = nscore_1.n_cut__;
					for (i__ = 1; i__ <= i__4; ++i__) {
						if (nscore_1.i_ali__[i__ - 1] == k_ali__[i__ - 1]) {
							++neq;
						}
					}
					if (nscore_1.n_cut__ == neq) {
						goto L302;
					}
				}
				/* L301: */
			}
			/* for iteration */
			L302:
			/* L300: */
			;
		}
		/* for shift */
		/* L333: */
	}
	/* ******* return the final rotation **************** */
	/* for initial length, L_ali/M */
	ll = 0;
	i__1 = ka0;
	for (i__ = 1; i__ <= i__1; ++i__) {
		m = k_ali0__[i__ - 1];
		/* record of the best alignment */
		r_1__[i__ * 3 - 3] = xa[align_1.ia[m - 1] - 1];
		r_1__[i__ * 3 - 2] = ya[align_1.ia[m - 1] - 1];
		r_1__[i__ * 3 - 1] = za[align_1.ia[m - 1] - 1];
		r_2__[i__ * 3 - 3] = stru_1.xb[align_1.ib[m - 1] - 1];
		r_2__[i__ * 3 - 2] = stru_1.yb[align_1.ib[m - 1] - 1];
		r_2__[i__ * 3 - 1] = stru_1.zb[align_1.ib[m - 1] - 1];
		++ll;
	}
	u3b_(w, r_1__, r_2__, &ll, &c__1, &rms, u, t, &ier);
	/* u rotate r_1 to r_2 */
	i__1 = nres_1.nseqa;
	for (j = 1; j <= i__1; ++j) {
		x1[j] = t[0] + u[0] * xa[j - 1] + u[3] * ya[j - 1] + u[6] * za[j - 1];
		y1[j] = t[1] + u[1] * xa[j - 1] + u[4] * ya[j - 1] + u[7] * za[j - 1];
		z1[j] = t[2] + u[2] * xa[j - 1] + u[5] * ya[j - 1] + u[8] * za[j - 1];
	}
	*tm = score_max__;
	/* ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* tmscore_ */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/*     1, collect those residues with dis<d; */
/*     2, calculate score_GDT, score_maxsub, score_TM */
/* cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */int score_fun__(void) {
	/* System generated locals */
	int i__1;
	float r__1, r__2, r__3;
	/* Builtin functions */
	double sqrt(double);
	/* Local variables */
	int i__, j, k;
	float score_sum__, dis, d_tmp__;
	/* [1,n_ali],align residues for the */
	d_tmp__ = para_1.d__;
	L21:
	nscore_1.n_cut__ = 0;
	/* number of residue-pairs dis<d, for iter */
	score_sum__ = 0.f;
	/* TMscore */
	i__1 = align_1.n_ali__;
	for (k = 1; k <= i__1; ++k) {
		i__ = align_1.ia[k - 1];
		/* [1,nseqA] reoder number of structureA */
		j = align_1.ib[k - 1];
		/* [1,nseqB] */
		/* Computing 2nd power */
		r__1 = stru_2.xa[i__ - 1] - stru_2.xb[j - 1];
		/* Computing 2nd power */
		r__2 = stru_2.ya[i__ - 1] - stru_2.yb[j - 1];
		/* Computing 2nd power */
		r__3 = stru_2.za[i__ - 1] - stru_2.zb[j - 1];
		dis = sqrt(r__1 * r__1 + r__2 * r__2 + r__3 * r__3);
		if (dis < d_tmp__) {
			++nscore_1.n_cut__;
			nscore_1.i_ali__[nscore_1.n_cut__ - 1] = k;
			/* [1,n_ali], mark the residue-pairs in di */
		}
		/* Computing 2nd power */
		r__1 = dis / para_1.d0;
		score_sum__ += 1 / (r__1 * r__1 + 1);
	}
	if (nscore_1.n_cut__ < 3 && align_1.n_ali__ > 3) {
		d_tmp__ += .5f;
		goto L21;
	}
	scores_1.score = score_sum__ / (float) nres_1.nseqb;
	/* TM-score */
	return 0;
} /* score_fun__ */
/* ******************************************************************* */
/*     Dynamic programming for alignment. */
/*     Input: score(i,j), and gap_open */
/*     Output: invmap(j) */
/*     Please note this subroutine is not a correct implementation of */
/*     the N-W dynamic programming because the score tracks back only */
/*     one layer of the matrix. This code was exploited in TM-align */
/*     because it is about 1.5 times faster than a complete N-W code */
/*     and does not influence much the final structure alignment result. */
/* ******************************************************************* */
/* Subroutine */int dp_(int *nseq1, int *nseq2, dpc_ *dpc_1) {
	/* System generated locals */
	int i__1, i__2;
	/* Local variables */
	float d__, h__;
	int i__, j;
	float v;
	//int dir[25010001] /* was [5001][5001] */;
	//float val[25010001] /* was [5001][5001] */;
	int *dir = (int *)malloc(sizeof(int)*25010001);
	float *val = (float *)malloc(sizeof(float)*25010001);
	/* **   initialize the matrix: */
	val[0] = 0.f;
	i__1 = *nseq1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		dir[i__] = FALSE_;
		val[i__] = 0.f;
	}
	i__1 = *nseq2;
	for (j = 1; j <= i__1; ++j) {
		dir[j * 5001] = FALSE_;
		val[j * 5001] = 0.f;
		dpc_1->invmap[j - 1] = -1;
	}
	/* **   decide matrix and path: */
	i__1 = *nseq2;
	for (j = 1; j <= i__1; ++j) {
		i__2 = *nseq1;
		for (i__ = 1; i__ <= i__2; ++i__) {
			d__ = val[i__ - 1 + (j - 1) * 5001]
					+ dpc_1->score[i__ + j * 5000 - 5001];
			h__ = val[i__ - 1 + j * 5001];
			if (dir[i__ - 1 + j * 5001]) {
				h__ += dpc_1->gap_open__;
			}
			v = val[i__ + (j - 1) * 5001];
			if (dir[i__ + (j - 1) * 5001]) {
				v += dpc_1->gap_open__;
			}
			if (d__ >= h__ && d__ >= v) {
				dir[i__ + j * 5001] = TRUE_;
				val[i__ + j * 5001] = d__;
			} else {
				dir[i__ + j * 5001] = FALSE_;
				if (v >= h__) {
					val[i__ + j * 5001] = v;
				} else {
					val[i__ + j * 5001] = h__;
				}
			}
		}
	}
	/* **   extract the alignment: */
	i__ = *nseq1;
	j = *nseq2;
	while (i__ > 0 && j > 0) {
		if (dir[i__ + j * 5001]) {
			dpc_1->invmap[j - 1] = i__;
			--i__;
			--j;
		} else {
			h__ = val[i__ - 1 + j * 5001];
			if (dir[i__ - 1 + j * 5001]) {
				h__ += dpc_1->gap_open__;
			}
			v = val[i__ + (j - 1) * 5001];
			if (dir[i__ + (j - 1) * 5001]) {
				v += dpc_1->gap_open__;
			}
			if (v >= h__) {
				--j;
			} else {
				--i__;
			}
		}
	}
	free(dir);free(val);
	/* ^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^ */
	return 0;
} /* dp_ */
/* ccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc */
/*  w    - w(m) is weight for atom pair  c m                 (given) */
/*  x    - x(i,m) are coordinates of atom c m in set x       (given) */
/*  y    - y(i,m) are coordinates of atom c m in set y       (given) */
/*  n    - n is number of atom pairs                         (given) */
/*  mode  - 0:calculate rms only                             (given) */
/*          1:calculate rms,u,t                              (takes longer) */
/*  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result) */
/*  u    - u(i,j) is   rotation  matrix for best superposition  (result) */
/*  t    - t(i)   is translation vector for best superposition  (result) */
/*  ier  - 0: a unique optimal superposition has been determined(result) */
/*       -1: superposition is not unique but optimal */
/*       -2: no result obtained because of negative weights w */
/*           or all weights equal to zero. */
/* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
/* Subroutine */int u3b_(double *w, double *x, double *y, int *n, int *mode,
		double *rms, double *u, double *t, int *ier) {
	/* Initialized data */
	double sqrt3 = 1.73205080756888;
	double tol = .01;
	double zero = 0.;
	int ip[9] = { 1, 2, 4, 2, 3, 5, 4, 5, 6 };
	int ip2312[4] = { 2, 3, 1, 2 };
	/* System generated locals */
	int i__1;
	double d__1, d__2;
	/* Builtin functions */
	double sqrt(double), atan2(double, double), cos(double), sin(double);
	/* Local variables */
	double a[9] /* was [3][3] */, b[9] /* was [3][3] */, d__, e[3], g,
			h__;
	int i__, j, k, l, m;
	double p, r__[9] /* was [3][3] */, e0;
	int m1;
	double wc, xc[3], yc[3], rr[6], ss[6], cof, det, cth, sth, spur,
			sigma, sqrth;
	/* Parameter adjustments */
	--t;
	u -= 4;
	y -= 4;
	x -= 4;
	--w;
	/* Function Body */
	wc = zero;
	*rms = zero;
	e0 = zero;
	for (i__ = 1; i__ <= 3; ++i__) {
		xc[i__ - 1] = zero;
		yc[i__ - 1] = zero;
		t[i__] = zero;
		for (j = 1; j <= 3; ++j) {
			r__[i__ + j * 3 - 4] = zero;
			u[i__ + j * 3] = zero;
			a[i__ + j * 3 - 4] = zero;
			if (i__ == j) {
				u[i__ + j * 3] = 1.f;
				a[i__ + j * 3 - 4] = 1.f;
			}
		}
	}
	*ier = -1;
	if (*n < 1) {
		return 0;
	}
	*ier = -2;
	i__1 = *n;
	for (m = 1; m <= i__1; ++m) {
		if (w[m] < 0.f) {
			return 0;
		}
		wc += w[m];
		for (i__ = 1; i__ <= 3; ++i__) {
			xc[i__ - 1] += w[m] * x[i__ + m * 3];
			yc[i__ - 1] += w[m] * y[i__ + m * 3];
		}
	}
	if (wc <= zero) {
		return 0;
	}
	for (i__ = 1; i__ <= 3; ++i__) {
		xc[i__ - 1] /= wc;
		yc[i__ - 1] /= wc;
	}
	i__1 = *n;
	for (m = 1; m <= i__1; ++m) {
		for (i__ = 1; i__ <= 3; ++i__) {
			/* Computing 2nd power */
			d__1 = x[i__ + m * 3] - xc[i__ - 1];
			/* Computing 2nd power */
			d__2 = y[i__ + m * 3] - yc[i__ - 1];
			e0 += w[m] * (d__1 * d__1 + d__2 * d__2);
			d__ = w[m] * (y[i__ + m * 3] - yc[i__ - 1]);
			for (j = 1; j <= 3; ++j) {
				r__[i__ + j * 3 - 4] += d__ * (x[j + m * 3] - xc[j - 1]);
			}
		}
	}
	det = r__[0] * (r__[4] * r__[8] - r__[7] * r__[5])
			- r__[3] * (r__[1] * r__[8] - r__[7] * r__[2])
			+ r__[6] * (r__[1] * r__[5] - r__[4] * r__[2]);
	sigma = det;
	m = 0;
	for (j = 1; j <= 3; ++j) {
		i__1 = j;
		for (i__ = 1; i__ <= i__1; ++i__) {
			++m;
			rr[m - 1] = r__[i__ * 3 - 3] * r__[j * 3 - 3]
					+ r__[i__ * 3 - 2] * r__[j * 3 - 2]
					+ r__[i__ * 3 - 1] * r__[j * 3 - 1];
		}
	}
	spur = (rr[0] + rr[2] + rr[5]) / 3.f;
	cof = (rr[2] * rr[5] - rr[4] * rr[4] + rr[0] * rr[5] - rr[3] * rr[3]
			+ rr[0] * rr[2] - rr[1] * rr[1]) / 3.f;
	det *= det;
	for (i__ = 1; i__ <= 3; ++i__) {
		e[i__ - 1] = spur;
	}
	if (spur <= zero) {
		goto L40;
	}
	d__ = spur * spur;
	h__ = d__ - cof;
	g = (spur * cof - det) / 2.f - spur * h__;
	if (h__ <= zero) {
		if (*mode == 0) {
			goto L50;
		} else {
			goto L30;
		}
	}
	sqrth = sqrt(h__);
	d__ = h__ * h__ * h__ - g * g;
	if (d__ < zero) {
		d__ = zero;
	}
	d__ = atan2(sqrt(d__), -g) / 3.f;
	cth = sqrth * cos(d__);
	sth = sqrth * sqrt3 * sin(d__);
	e[0] = spur + cth + cth;
	e[1] = spur - cth + sth;
	e[2] = spur - cth - sth;
	if (*mode == 0) {
		goto L50;
	}
	for (l = 1; l <= 3; l += 2) {
		d__ = e[l - 1];
		ss[0] = (d__ - rr[2]) * (d__ - rr[5]) - rr[4] * rr[4];
		ss[1] = (d__ - rr[5]) * rr[1] + rr[3] * rr[4];
		ss[2] = (d__ - rr[0]) * (d__ - rr[5]) - rr[3] * rr[3];
		ss[3] = (d__ - rr[2]) * rr[3] + rr[1] * rr[4];
		ss[4] = (d__ - rr[0]) * rr[4] + rr[1] * rr[3];
		ss[5] = (d__ - rr[0]) * (d__ - rr[2]) - rr[1] * rr[1];
		if (abs(ss[0]) >= abs(ss[2])) {
			j = 1;
			if (abs(ss[0]) < abs(ss[5])) {
				j = 3;
			}
		} else if (abs(ss[2]) >= abs(ss[5])) {
			j = 2;
		} else {
			j = 3;
		}
		d__ = zero;
		j = (j - 1) * 3;
		for (i__ = 1; i__ <= 3; ++i__) {
			k = ip[i__ + j - 1];
			a[i__ + l * 3 - 4] = ss[k - 1];
			d__ += ss[k - 1] * ss[k - 1];
		}
		if (d__ > zero) {
			d__ = 1.f / sqrt(d__);
		}
		for (i__ = 1; i__ <= 3; ++i__) {
			a[i__ + l * 3 - 4] *= d__;
		}
	}
	d__ = a[0] * a[6] + a[1] * a[7] + a[2] * a[8];
	if (e[0] - e[1] > e[1] - e[2]) {
		m1 = 3;
		m = 1;
	} else {
		m1 = 1;
		m = 3;
	}
	p = zero;
	for (i__ = 1; i__ <= 3; ++i__) {
		a[i__ + m1 * 3 - 4] -= d__ * a[i__ + m * 3 - 4];
		/* Computing 2nd power */
		d__1 = a[i__ + m1 * 3 - 4];
		p += d__1 * d__1;
	}
	if (p <= tol) {
		p = 1.f;
		for (i__ = 1; i__ <= 3; ++i__) {
			if (p < (d__1 = a[i__ + m * 3 - 4], abs(d__1))) {
				goto L21;
			}
			p = (d__1 = a[i__ + m * 3 - 4], abs(d__1));
			j = i__;
			L21: ;
		}
		k = ip2312[j - 1];
		l = ip2312[j];
		/* Computing 2nd power */
		d__1 = a[k + m * 3 - 4];
		/* Computing 2nd power */
		d__2 = a[l + m * 3 - 4];
		p = sqrt(d__1 * d__1 + d__2 * d__2);
		if (p <= tol) {
			goto L40;
		}
		a[j + m1 * 3 - 4] = zero;
		a[k + m1 * 3 - 4] = -a[l + m * 3 - 4] / p;
		a[l + m1 * 3 - 4] = a[k + m * 3 - 4] / p;
	} else {
		p = 1.f / sqrt(p);
		for (i__ = 1; i__ <= 3; ++i__) {
			a[i__ + m1 * 3 - 4] *= p;
		}
	}
	a[3] = a[7] * a[2] - a[1] * a[8];
	a[4] = a[8] * a[0] - a[2] * a[6];
	a[5] = a[6] * a[1] - a[0] * a[7];
	L30: for (l = 1; l <= 2; ++l) {
		d__ = zero;
		for (i__ = 1; i__ <= 3; ++i__) {
			b[i__ + l * 3 - 4] = r__[i__ - 1] * a[l * 3 - 3]
					+ r__[i__ + 2] * a[l * 3 - 2] + r__[i__ + 5] * a[l * 3 - 1];
			/* Computing 2nd power */
			d__1 = b[i__ + l * 3 - 4];
			d__ += d__1 * d__1;
		}
		if (d__ > zero) {
			d__ = 1.f / sqrt(d__);
		}
		for (i__ = 1; i__ <= 3; ++i__) {
			b[i__ + l * 3 - 4] *= d__;
		}
	}
	d__ = b[0] * b[3] + b[1] * b[4] + b[2] * b[5];
	p = zero;
	for (i__ = 1; i__ <= 3; ++i__) {
		b[i__ + 2] -= d__ * b[i__ - 1];
		/* Computing 2nd power */
		d__1 = b[i__ + 2];
		p += d__1 * d__1;
	}
	if (p <= tol) {
		p = 1.f;
		for (i__ = 1; i__ <= 3; ++i__) {
			if (p < (d__1 = b[i__ - 1], abs(d__1))) {
				goto L22;
			}
			p = (d__1 = b[i__ - 1], abs(d__1));
			j = i__;
			L22: ;
		}
		k = ip2312[j - 1];
		l = ip2312[j];
		/* Computing 2nd power */
		d__1 = b[k - 1];
		/* Computing 2nd power */
		d__2 = b[l - 1];
		p = sqrt(d__1 * d__1 + d__2 * d__2);
		if (p <= tol) {
			goto L40;
		}
		b[j + 2] = zero;
		b[k + 2] = -b[l - 1] / p;
		b[l + 2] = b[k - 1] / p;
	} else {
		p = 1.f / sqrt(p);
		for (i__ = 1; i__ <= 3; ++i__) {
			b[i__ + 2] *= p;
		}
	}
	b[6] = b[1] * b[5] - b[4] * b[2];
	b[7] = b[2] * b[3] - b[5] * b[0];
	b[8] = b[0] * b[4] - b[3] * b[1];
	for (i__ = 1; i__ <= 3; ++i__) {
		for (j = 1; j <= 3; ++j) {
			u[i__ + j * 3] = b[i__ - 1] * a[j - 1] + b[i__ + 2] * a[j + 2]
					+ b[i__ + 5] * a[j + 5];
		}
	}
	L40: for (i__ = 1; i__ <= 3; ++i__) {
		t[i__] = yc[i__ - 1] - u[i__ + 3] * xc[0] - u[i__ + 6] * xc[1]
				- u[i__ + 9] * xc[2];
	}
	L50: for (i__ = 1; i__ <= 3; ++i__) {
		if (e[i__ - 1] < zero) {
			e[i__ - 1] = zero;
		}
		e[i__ - 1] = sqrt(e[i__ - 1]);
	}
	*ier = 0;
	if (e[1] <= e[0] * 1e-5) {
		*ier = -1;
	}
	d__ = e[2];
	if (sigma < 0.f) {
		d__ = -d__;
		if (e[1] - e[2] <= e[0] * 1e-5) {
			*ier = -1;
		}
	}
	d__ = d__ + e[1] + e[0];
	*rms = e0 - d__ - d__;
	if (*rms < 0.f) {
		*rms = 0.f;
	}
	return 0;
} /* u3b_ */

int pow_ii_f2c(int *ap, int *bp) {
	int pow, x, n;
	unsigned long u;
	x = *ap;
	n = *bp;

	if (n <= 0) {
		if (n == 0 || x == 1)
			return 1;
		if (x != -1)
			return x == 0 ? 1 / x : 0;
		n = -n;
	}
	u = n;
	for (pow = 1;;) {
		if (u & 01)
			pow *= x;
		if (u >>= 1)
			x *= x;
		else
			break;
	}
	return (pow);
}

double pow_dd(double *a, double *b) {
	return pow(*a, *b);
}
