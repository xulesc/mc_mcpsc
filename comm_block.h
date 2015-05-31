//
#ifndef _COMM_BLOCK_
#define _COMM_BLOCK_
typedef struct {
        short die;
        int job_index;
        int algo_type;
} MasterToClientTransferBlock;

typedef struct {
        long seconds_at_job_start, seconds_at_result_send;
        
        int ready, job_index, aln_len, len1, len2, algo_type;
        float rmsd;
        // milliseconds
        long processing_time, data_collect_time, idle_time;
        //ce
        float z;
        int gaps;
        //tmalign
        float tm1, tm2;
        float seq_id;
        int n_eq__;
} ClientToMasterTransferBlock;
#endif

//