//
#include "comm_block.h"
#include <sys/timeb.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef _TMALIGN_EXTERN_
#define _TMALIGN_EXTERN_

#define MAX_ATOMS 5000

typedef struct {
	float xa[30000] /* was [3][5000][2] */;
} backbone_;
// Function Prototypes
// entry point for external codes
extern int main_process(char *ss1, char *ss2, char *seq1, char *seq2,
		ClientToMasterTransferBlock *master_in_data, backbone_ *backbone_1,
		int *atom_indices1, int *atom_indices2, int len1, int len2);
#endif
//
