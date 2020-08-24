#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "bseq.h"
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char seq_nt4_table[256];

struct bseq_file_s {
	int is_eof;
	gzFile fp;
	kseq_t *ks;
};

bseq_file_t *bseq_open(const char *fn)
{
	bseq_file_t *fp;
	gzFile f;
	f = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	if (f == 0) return 0;
	fp = (bseq_file_t*)calloc(1, sizeof(bseq_file_t));
	fp->fp = f;
	fp->ks = kseq_init(fp->fp);
	return fp;
}

void bseq_close(bseq_file_t *fp)
{
	kseq_destroy(fp->ks);
	gzclose(fp->fp);
	free(fp);
}

bseq1_t *bseq_read(bseq_file_t *fp, size_t *n_)
{
	size_t m, n;
	bseq1_t *seqs;
	kseq_t *ks = fp->ks;
	m = n = 0; seqs = 0;
	while (kseq_read(ks) >= 0) {
		bseq1_t *s;
		// assert(ks->seq.l <= INT32_MAX);
		if (n >= m) {
			m = m? m<<1 : (1 << 20);
			seqs = (bseq1_t*)realloc(seqs, m * sizeof(bseq1_t));
		}
		s = &seqs[n];
		s->seq = strdup(ks->seq.s);
		s->name = strdup(ks->name.s);
		
		s->qual = NULL;
		if (ks->qual.s) s->qual = strdup(ks->qual.s);
		// fprintf(stderr, "%s\n", s->seq);
		// if (qualfp) fprintf(qualfp, "%s\n", ks->qual.s);
		// fprintf(stderr, "%s\n", ks->name.s);
		// if (idsfp) fprintf(idsfp, "%s\n", ks->name.s);

		// if (seq_len != (int)ks->seq.l){
		// 	fprintf(stderr, "Length of reads are different. The program can not compress it.\n");
		// 	exit(1);
		// }
		
		n++;
		// if (n >= 10000000) break;//for debug
		// if (n >= 100000) break;//for debug
		// if (n > 10) exit(0);
	}
	if (n == 0) fp->is_eof = 1;
	*n_ = n;
	return seqs;
}

// void bseq_read_second(bseq1_t *&seqs, bseq_file_t *fp, int *n_, int seq_len)
// {
// 	int m, n;
// 	kseq_t *ks = fp->ks;
// 	n = *n_;
// 	m = 256; 
// 	while (m <= n) {
// 		m <<= 1;
// 	}
// 	// fprintf(stderr, "m: %d, n: %d;\n", m, n);

// 	while (kseq_read(ks) >= 0) {
// 		bseq1_t *s;
// 		// assert(ks->seq.l <= INT32_MAX);
// 		if (n >= m) {
// 			m = m? m<<1 : 256;
// 			seqs = (bseq1_t*)realloc(seqs, m * sizeof(bseq1_t));
// 		}
// 		s = &seqs[n];
// 		// s->name = strdup(ks->name.s);
// 		s->seq = strdup(ks->seq.s);
// 		if (seq_len != (int)ks->seq.l){
// 			fprintf(stderr, "Length of reads are different. The program can not compress it.\n");
// 			exit(1);
// 		}
		
// 		n++;
// 		// if (n >= 10000000) break;//for debug
// 		// if (n >= 10000) break;//for debug
// 	}
// 	if (n == 0) fp->is_eof = 1;
// 	*n_ = n;
// 	// return seqs;
// }

int bseq_read_len(bseq_file_t *fp)
{
	size_t n = 0;
	kseq_t *ks = fp->ks;
	while (kseq_read(ks) >= 0) {
		n++;
		break;
	}
	if (n == 0) fp->is_eof = 1;
	return (int)ks->seq.l;
}

int bseq_eof(bseq_file_t *fp)
{
	return fp->is_eof;
}
