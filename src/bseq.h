#ifndef MM_BSEQ_H
#define MM_BSEQ_H

#include <stdint.h>
#include <stdio.h>

struct bseq_file_s;
typedef struct bseq_file_s bseq_file_t;

typedef struct {
	// int rid;//, n_num, m_m, *pos_n;//rid: index of reads; n_num: number of N char; m_m: array size of pos_n; pos_n: store positon of N
	char *name, *seq, *qual;
} bseq1_t;

bseq_file_t *bseq_open(const char *fn);
void bseq_close(bseq_file_t *fp);
// bseq1_t *bseq_read(bseq_file_t *fp, size_t *n_, int seq_len);
bseq1_t *bseq_read(bseq_file_t *fp, size_t *n_);
// void bseq_read_second(bseq1_t *&seqs, bseq_file_t *fp, int *n_, int seq_len);
int bseq_read_len(bseq_file_t *fp);
int bseq_eof(bseq_file_t *fp);

#endif
