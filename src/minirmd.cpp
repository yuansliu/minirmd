# include <stdio.h>
# include <iostream>
# include <fstream>
# include <vector>
# include <queue>
# include <assert.h>
# include <chrono>
# include <thread>
# include <mutex>
# include <stdlib.h>
# include <string.h>
# include <time.h>
# include <algorithm>
# include <unistd.h>
# include "minirmd.h"

using namespace std;

int nthreads = 24;
int L;
size_t max_rid;
bseq1_t *seq, *seq1, *seq2;
bool *isremove;
int64_t *indexremove;
int difthr = 0;
bool isRC = false;
bool isPE = false;
bool iskf;
int *avgQual;
string rf1, rf2, rsf, kf, logpath;
string header_prefix;
const size_t BUFFER_SIZE = 1ULL << 25;

// bool debug = false;
// size_t did, rcdid;

char complement[256];

uint64_t *ridvec;
size_t ridvecid; // init = 0 before run sortBuckets();
segmemt_v sgv;
mutex segmtx;

mm128_t *minis;
mm192_t *lminis;


inline bool getReads(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	seq1 = bseq_read(fp, &max_rid);
	bseq_close(fp);
	return true;
}

inline bool getReads2(char const *fn) {
	bseq_file_t *fp;
	fp = bseq_open(fn);
	if (fp == 0) return false;
	size_t aa = max_rid;
	seq2 = bseq_read(fp, &max_rid);
	bseq_close(fp);
	if (aa != max_rid) {
		fprintf(stderr, "The two files contain different number of reads. EXIT\n");
		exit(0);
	}
	return true;
}

inline void init() {
	memset(complement, 125, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';
}

inline void displayHelp(const char *prog)
{
	printf("minirmd v1, by Yuansheng Liu, October 2020.\n");
	printf("Usage: %s -i <file> -f <file> -o <output> [option parameters]\n", prog);
	printf("\t options:\n");
	// printf("-----------\n");
	printf("\t\t -i reads file\n");
	printf("\t\t -f reads file, if paired end\n");
	printf("\t\t -o the output file\n");
	printf("\t\t -l the log file location\n");
	printf("\t\t -d number of allowed mismatch\n");
	printf("\t\t -k the file to store values of k\n");
	printf("\t\t -r remove duplicates on reverse-complement strand\n");
	printf("\t\t -t the number of threads\n");
	printf("\t\t -h print help message\n");

	printf("Example:\n\t\t");
	printf("./minirmd -i test.fastq -o test_rm_1.fastq -d 1\n\t\t");
	printf("./minirmd -i test_1.fastq -f test_2.fastq -o test_rm_2.fastq -d 2\n\n");
}

inline void getPars(int argc, char* argv[]) {
	bool is1 = false, is2 = false; //four
	int oc;
	iskf = false;
	while ((oc = getopt(argc, argv, "i:f:o:l:t:d:k:hr")) >= 0) {
		switch (oc) {
			case 'i':
				rf1 = optarg;
				is1 = true;
				break;
			case 'f':
				rf2 = optarg;
				isPE = true;
				break;
			case 'o':
				rsf = optarg;
				is2 = true;
				break;
			case 'k':
				kf = optarg;
				iskf = true;
				break;
			case 'r':
				isRC = true;
				break;
			case 'l':
				logpath = optarg;
				break;
			case 't':
				nthreads = atoi(optarg);
				break;
			case 'd':
				difthr = atoi(optarg);
				break;
			case 'h':
				displayHelp(argv[0]);
				exit(0);
			case '?':
				std::cerr << "Error parameters.\n Please run 'xxx -h'\n";
				exit(1);
				break;
		}
	}

	if (!is1 || !is2) {
		// fprintf(stderr, "Required parameters are not provided!!\n\n");
		displayHelp(argv[0]);
		exit(1);
	}
	
	std::ifstream f;

	if (iskf) {
		f.open(kf);
		if (f.fail()) {
			displayHelp(argv[0]);
			// fprintf(stderr, "Kmer file '%s' does not exist.\n", kf.c_str());
			exit(1);
		}
		f.close();
	}

	f.open(rf1);
	if (f.fail()) {
		displayHelp(argv[0]);
		// fprintf(stderr, "Reads file '%s' does not exist.\n", rf1.c_str());
		exit(1);
	}
	f.close();

	if (isPE) {
		f.open(rf2);
		if (f.fail()) {
			displayHelp(argv[0]);
			// fprintf(stderr, "Reads file '%s' does not exist.\n", rf2.c_str());
			exit(1);
		}
		f.close();
	}

    if (logpath.length() == 0)
	{
		logpath = "output.log";
	}

	header_prefix = ">";

	if (rsf.find(".fq") != -1 || rsf.find("fastaq") != -1 || rsf.find("fastq") != -1)
	{
		header_prefix = "@";
	}

	/*std::ofstream fo;
	fo.open(rsf);
	if (fo.fail()) {
		fprintf(stderr, "The output file '%s' can not be created.\n", rsf.c_str());
		exit(1);
	}
	fo.close();*/
}

inline void reverseComplement(char* start, int len) {
	char* left = start; // sequence starts
	char* right = start + len - 1;
	while (right > left) {
		char tmp = complement[(uint8_t)*left];
		*left = complement[(uint8_t)*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = complement[(uint8_t)*left];
}

inline void reverseComplementPE(char* start, int len) {
	reverseComplement(start, len>>1);
	reverseComplement(start + (len>>1), len>>1);
}

static inline uint64_t hash64(uint64_t key, uint64_t mask) {
	key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
	key = key ^ key >> 24;
	key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
	key = key ^ key >> 14;
	key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
	key = key ^ key >> 28;
	key = (key + (key << 31)) & mask;
	return key;
}

size_t rid_pthread;

// --- for minimizer
int bsize = 21;
mm128_v *B;
mm192_v *BL;

size_t *Bsize, *BLsize;

// int max_kmer = 29;
int *kmervec, *lkmervec, kmervecsize, lkmervecsize;
// mutex *bmtx;

void calcMinimizersFun() { // before calling this method, MUST set rid_pthread = 0
	// char *str = (char*)alloca((L + 3) * sizeof(char));
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));
	int i, l, k, c;
	mm128_t min, info;
	uint64_t shift1, mask, kmer;
	int Bmask = (1<<bsize) - 1, bucketidx;
	
	while (1) {
		size_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		// cout << "rid: " << rid << endl;
		// minis[rid] = new mm128_t[kmervecsize];

		// strcpy(str, seq[rid].seq);

		// if (isRC) {
		// 	strcpy(rcstr, seq[rid].seq);
		// 	reverseComplement(rcstr, L);
		// }
		for (int id = 0; id < kmervecsize; ++id) {
			// mm_sketch(seq[rid].seq, L, kmervec[i], rid, &minimizer);
			if (isPE) {
				if (id&1) {
					seq = seq2;
					// fprintf(stderr, "id: %d; seq2\n", id);
				} else {
					// fprintf(stderr, "id: %d; seq1\n", id);
					seq = seq1;
				}
			}

			k = kmervec[id];
			shift1 = 2 * (k - 1);
			mask = (1ULL<<2*k) - 1;
			kmer = 0;
			min = { UINT64_MAX, UINT64_MAX };

			for (i = 0, l = 0; i < L; i++) {
				c = seq_nt4_table[(uint8_t)seq[rid].seq[i]];
				info = { UINT64_MAX, UINT64_MAX };
				kmer = (kmer << 2 | c) & mask;           // forward k-mer
				if (++l >= k) {
					info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i << 1 | 0;
				}
				if (info.x < min.x) {
					min = info;
				}
			}
			
			if (isRC) {
				// mm_sketch_rc(rcstr, L, kmervec[i], rid, &minimizer);
				if (id == 0 || isPE) {
					strcpy(rcstr, seq[rid].seq);
					reverseComplement(rcstr, L);
				}
				kmer = 0;
				for (i = 0, l = 0; i < L; i++) {
					c = seq_nt4_table[(uint8_t)rcstr[i]];
					info = { UINT64_MAX, UINT64_MAX };
					kmer = (kmer << 2 | c) & mask;           // forward k-mer
					if (++l >= k) {
						info.x = hash64(kmer, mask), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
					}
					if (info.x < min.x) {
						min = info;
					}
				}
			}

			minis[rid*kmervecsize + id] = min;

			bucketidx = (min.x & Bmask) + id*(1<<bsize);
			__sync_fetch_and_add(&Bsize[bucketidx], 1);
		}

		if (lkmervecsize == 0 && seq1[rid].qual != NULL) {
			int sum = 0;
			char *q = seq1[rid].qual;
			while (*q != '\0') {
				sum += *q;
				++q;
			}
			if (isPE && seq2[rid].qual != NULL) {
				q = seq2[rid].qual;
				while (*q != '\0') {
					sum += *q;
					++q;
				}
			}
			avgQual[rid] = sum;
		}
	}
}

inline void calcMinimizers() {
	minis = new mm128_t[max_rid*kmervecsize];
	// for (size_t rid = 0; rid < max_rid; ++rid) {
	// 	minis[rid] = new mm128_t[kmervecsize];
	// }
	// cout << seq[0].seq << endl;
	rid_pthread = 0;
	std::vector<thread> threadVec;
	// cout << "nthreads: " << nthreads << endl;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(calcMinimizersFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
}

void calcLongMinimizersFun() { // before calling this method, MUST set rid_pthread = 0
	// char *str = (char*)alloca((L + 3) * sizeof(char));
	char *rcstr = (char*)alloca((L + 3) * sizeof(char));
	int i, l, k, c;
	uint64_t maskX, maskY;
	mm128_t kmer = {0, 0};
	mm192_t min(UINT64_MAX, UINT64_MAX, UINT64_MAX), info;
	int Bmask = (1<<bsize) - 1, bucketidx;
	
	while (1) {
		size_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;
		
		for (int id = 0; id < lkmervecsize; ++id) {
			// mm_sketch(seq[rid].seq, L, kmervec[i], rid, &minimizer);
			if (isPE) {
				if (id&1) {
					seq = seq2;
					// fprintf(stderr, "id: %d; seq2\n", id);
				} else {
					// fprintf(stderr, "id: %d; seq1\n", id);
					seq = seq1;
				}
			}

			k = lkmervec[id];
			maskX = (1ULL<<(2*31)) - 1;
			maskY = (1ULL<<(2*(k-31))) - 1;
			kmer = {0, 0};
			min = { UINT64_MAX, UINT64_MAX, UINT64_MAX };

			for (i = 0, l = 0; i < L; i++) {
				c = seq_nt4_table[(uint8_t)seq[rid].seq[i]];
				info = { UINT64_MAX, UINT64_MAX, UINT64_MAX };
				kmer.y = ((kmer.y << 2) | (kmer.x >> 60)) & maskY;
				kmer.x = (kmer.x << 2 | c) & maskX;           // forward k-mer

				if (++l >= k) {
					info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 0;
				}
				if (info.x < min.x) {
					min = info;
				}
			}
			
			if (isRC) {
				// mm_sketch_rc(rcstr, L, kmervec[i], rid, &minimizer);
				if (id == 0 || isPE) {
					strcpy(rcstr, seq[rid].seq);
					reverseComplement(rcstr, L);
				}
				kmer = {0, 0};
				for (i = 0, l = 0; i < L; i++) {
					c = seq_nt4_table[(uint8_t)rcstr[i]];
					info = { UINT64_MAX, UINT64_MAX, UINT64_MAX };
					kmer.y = ((kmer.y << 2) | (kmer.x >> 60)) & maskY;
					kmer.x = (kmer.x << 2 | c) & maskX;           // forward k-mer

					if (++l >= k) {
						info.x.x = hash64(kmer.x, maskX), info.x.y = hash64(kmer.y, maskY), info.y = (uint64_t)rid<<32 | (uint32_t)i<<1 | 1;
					}
					if (info.x < min.x) {
						min = info;
					}
				}
			}
			lminis[rid*lkmervecsize + id] = min;

			bucketidx = (((min.x.x & Bmask) + (min.x.y & Bmask)) & Bmask) + id*(1<<bsize);
			__sync_fetch_and_add(&BLsize[bucketidx], 1);
		}

		if (seq1[rid].qual != NULL) {
			int sum = 0;
			char *q = seq1[rid].qual;
			while (*q != '\0') {
				sum += *q;
				++q;
			}
			if (isPE && seq2[rid].qual != NULL) {
				q = seq2[rid].qual;
				while (*q != '\0') {
					sum += *q;
					++q;
				}
			}
			avgQual[rid] = sum;
		}
	}
}

inline void calcLongMinimizers() {
	lminis = new mm192_t[max_rid*lkmervecsize];
	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(calcLongMinimizersFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
}

int ivcmp(const void *a_, const void *b_) {
	uint64_t a = *(uint64_t*)a_, b = *(uint64_t*)b_;
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));

	strcmp(stra, seq[(uint32_t)(a >> 32)].seq);
	if (a & 1) reverseComplement(stra, L);

	strcmp(strb, seq[(uint32_t)(b >> 32)].seq);
	if (b & 1) reverseComplement(strb, L);
	return strcmp(stra, strb);
}

void sortBucketsFun() {
	size_t j, start_a, n;
	// mmrec_t ttmmrec;
	segmemt_t sg;
	uint32_t pos_a;
	uint64_t aa;

	uint64_v *iv = new uint64_v[L];

	while (1) {
		size_t bid = __sync_fetch_and_add(&rid_pthread, 1);
		if (bid >= (1ul << bsize)) break;

		mm128_v *b = &B[bid];
		// cout << "bb: " << bb << endl;
		if (b->n > 0) {
			// cout << "b->n: " << b->n << endl;
			radix_sort_128x(b->a, b->a + b->n);

			for (j = 1, n = 1, start_a = 0; j <= b->n; ++j) {
				if (j == b->n || b->a[j].x != b->a[j-1].x) {
					// fprintf(stderr, "j: %lu;\nstart_a: %lu\nn: %lu", j, start_a, n);
					assert(j - start_a == n);

					if (n >= 2) { // tree
						for (int k = 0; k < L; ++k) {
							kv_init(iv[k]);
						}
						for (size_t k = 0; k < n; ++k) {
							aa = b->a[start_a + k].y;
							pos_a = (uint32_t)aa>>1;
							kv_push(uint64_t, iv[pos_a], aa);
						}

						for (int k = 0; k < L; ++k) {
							if (iv[k].n > 1) {

								if (iv[k].n > 20000 && iv[k].n < 100000) { //
									qsort(iv[k].a, iv[k].n, sizeof(uint64_t), ivcmp);
								}

								sg.a = __sync_fetch_and_add(&ridvecid, iv[k].n);
								// sg.a = ridvecid;
								for (int i1 = 0; i1 < iv[k].n; ++i1) {
									ridvec[sg.a + i1] = iv[k].a[i1];
								}
								// sg.b = ridvecid;
								sg.b = iv[k].n;
								// [sg.a, sg.b)
								segmtx.lock();
								kv_push(segmemt_t, sgv, sg);
								segmtx.unlock();
							}
						}

						for (int k = 0; k < L; ++k) {
							kv_destroy(iv[k]);
						}
					} 
					start_a = j, n = 1;
				} else {++n;}
			}
		}
	}
	delete[] iv;
}

bool segcmp(const segmemt_t &a, const segmemt_t &b) {
	return a.b > b.b;
}

int sgvnl; //the number of large clusters

void bucketsFun(int kmervecidx) {
	int Bmask = (1<<bsize) - 1, bucketidx;
	mm128_t minimizer;

	while (1) {
		size_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		if (!isremove[rid]) {
			minimizer = minis[rid*kmervecsize + kmervecidx];
			bucketidx = minimizer.x & Bmask;

			size_t n = __sync_fetch_and_add(&B[bucketidx].n, 1);
			B[bucketidx].a[n] = minimizer;
		}

	}
}

inline void sortBuckets(int kmervecidx) {
	int mask = (1<<bsize) - 1, bucketidx;
	mm128_t minimizer;

	for (int i = 0; i <= mask; ++i) {
		kv_resize(mm128_t, B[i], Bsize[kmervecidx*(1<<bsize) + i]);
	}

	if (max_rid > (1UL<<23)) {
	// if (true) {
		rid_pthread = 0;
		std::vector<thread> threadVec;
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(bucketsFun, kmervecidx));
		}
		std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
			thr.join();
		});
		threadVec.clear();
	} else {
		for (size_t rid = 0; rid < max_rid; ++rid) {
			if (!isremove[rid]) {
				// minimizer = minis[rid][kmervecidx];
				minimizer = minis[rid*kmervecsize + kmervecidx];
				bucketidx = minimizer.x & mask;
				mm128_v *p = &B[bucketidx];
				kv_push(mm128_t, *p, minimizer);
			}
		}
	}

	kv_init(sgv);
	kv_resize(segmemt_t, sgv, 1<<20);
	// kv_init(lsgv);
	// kv_resize(segmemt_t, lsgv, 1<<10);

	ridvecid = 0;

	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(sortBucketsFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();

	sort (sgv.a, sgv.a + sgv.n, segcmp);

	sgvnl = 0;
	for (size_t i = 0; i < sgv.n; ++i) {
		// if (sgv.a[i].b < 20000) break;
		if (sgv.a[i].b < 10000) break;
		++ sgvnl;
	}

	/*cout << "sgv.n: " << sgv.n << endl;
	for (size_t i = 0; i < 10; ++i) {
		cout << sgv.a[i].b << endl;
	}*/
}

bool mm192cmp(const mm192_t &a, const mm192_t &b) {
	return a.x.y < b.x.y || (a.x.y == b.x.y && a.x.x < b.x.x);
}

void sortLongBucketsFun() {
	size_t j, start_a, n;
	// mmrec_t ttmmrec;
	segmemt_t sg;
	uint32_t pos_a;
	uint64_t aa;

	uint64_v *iv = new uint64_v[L];

	while (1) {
		size_t bid = __sync_fetch_and_add(&rid_pthread, 1);
		if (bid >= (1ul << bsize)) break;

		mm192_v *b = &BL[bid];
		// cout << "bb: " << bb << endl;
		if (b->n > 0) {
			// cout << "b->n: " << b->n << endl;
			// radix_sort_128x(b->a, b->a + b->n);
			sort(b->a, b->a + b->n, mm192cmp);

			for (j = 1, n = 1, start_a = 0; j <= b->n; ++j) {
				if (j == b->n || b->a[j].x != b->a[j-1].x) {
					// fprintf(stderr, "j: %lu;\nstart_a: %lu\nn: %lu", j, start_a, n);
					assert(j - start_a == n);

					if (n >= 2) { // tree
						for (int k = 0; k < L; ++k) {
							kv_init(iv[k]);
						}
						for (size_t k = 0; k < n; ++k) {
							aa = b->a[start_a + k].y;
							pos_a = (uint32_t)aa>>1;
							kv_push(uint64_t, iv[pos_a], aa);
						}

						for (int k = 0; k < L; ++k) {
							if (iv[k].n > 1) {

								if (iv[k].n > 20000 && iv[k].n < 100000) { //
									qsort(iv[k].a, iv[k].n, sizeof(uint64_t), ivcmp);
								}

								sg.a = __sync_fetch_and_add(&ridvecid, iv[k].n);
								// sg.a = ridvecid;
								for (int i1 = 0; i1 < iv[k].n; ++i1) {
									ridvec[sg.a + i1] = iv[k].a[i1];
								}
								// sg.b = ridvecid;
								sg.b = iv[k].n;
								// [sg.a, sg.b)
								segmtx.lock();
								kv_push(segmemt_t, sgv, sg);
								segmtx.unlock();
							}
						}

						for (int k = 0; k < L; ++k) {
							kv_destroy(iv[k]);
						}
					} 
					start_a = j, n = 1;
				} else {++n;}
			}
		}
	}
	delete[] iv;
}

void longBucketsFun(int kmervecidx) {
	int Bmask = (1<<bsize) - 1, bucketidx;
	mm192_t minimizer;

	while (1) {
		size_t rid = __sync_fetch_and_add(&rid_pthread, 1);
		if (rid >= max_rid) break;

		if (!isremove[rid]) {
			minimizer = lminis[rid*lkmervecsize + kmervecidx];
			bucketidx = ((minimizer.x.x & Bmask) + (minimizer.x.y & Bmask)) & Bmask;

			size_t n = __sync_fetch_and_add(&BL[bucketidx].n, 1);
			BL[bucketidx].a[n] = minimizer;
		}

	}
}

inline void sortLongBuckets(int kmervecidx) {
	int mask = (1<<bsize) - 1, bucketidx;
	mm192_t minimizer;

	for (int i = 0; i <= mask; ++i) {
		kv_resize(mm192_t, BL[i], BLsize[kmervecidx*(1<<bsize) + i]);
	}

	if (max_rid > (1UL<<23)) {
	// if (true) {
		rid_pthread = 0;
		std::vector<thread> threadVec;
		for (int i = 0; i < nthreads; ++i) {
			threadVec.push_back(std::thread(longBucketsFun, kmervecidx));
		}
		std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
			thr.join();
		});
		threadVec.clear();
	} else {
		for (size_t rid = 0; rid < max_rid; ++rid) {
			if (!isremove[rid]) {
				// minimizer = minis[rid][kmervecidx];
				minimizer = lminis[rid*lkmervecsize + kmervecidx];
				bucketidx = ((minimizer.x.x & mask) + (minimizer.x.y & mask)) & mask;
				mm192_v *p = &BL[bucketidx];
				kv_push(mm192_t, *p, minimizer);
			}
		}
	}

	kv_init(sgv);
	kv_resize(segmemt_t, sgv, 1<<20);
	// kv_init(lsgv);
	// kv_resize(segmemt_t, lsgv, 1<<10);

	ridvecid = 0;

	rid_pthread = 0;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		threadVec.push_back(std::thread(sortLongBucketsFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();

	sort (sgv.a, sgv.a + sgv.n, segcmp);
	
	sgvnl = 0;
	for (size_t i = 0; i < sgv.n; ++i) {
		// if (sgv.a[i].b < 20000) break;
		if (sgv.a[i].b < 10000) break;
		++ sgvnl;
	}

	/*cout << "sgv.n: " << sgv.n << endl;
	for (size_t i = 0; i < 5; ++i) {
		cout << sgv.a[i].b << endl;
	}*/
}

void processClusterSEFun() {
	size_t beg, end;

	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));
	// cout << "text xxx\n";
	uint32_t aa, bb;
	int difnum;

	while (1) {
		size_t i = __sync_fetch_and_add(&rid_pthread, 1);
		if (i >= sgv.n) break;

		/*if (sgv.a[i].b > 10000) {
			segmtx.lock();
			kv_push(segmemt_t, lsgv, sgv.a[i]);
			segmtx.unlock();
			continue;
		}*/
		beg = sgv.a[i].a;
		end = sgv.a[i].a + sgv.a[i].b;

		if (!isRC) {
			for (size_t i1 = beg; i1 < end; ++i1) {
				aa = ridvec[i1]>>32;
				if (!isremove[aa]) {
					strcpy(stra, seq[aa].seq);
					for (size_t i2 = i1 + 1; i2 < end; ++i2) {
						bb = ridvec[i2]>>32;
						if (!isremove[bb]) {
							strcpy(strb, seq[bb].seq);
							difnum = 0;
							for (int ss = 0; ss < L; ++ss) {
								if (stra[ss] != strb[ss]) ++difnum;
							}
							if (difnum <= difthr) {
								if (avgQual[aa] >= avgQual[bb]) {
									isremove[bb] = true;
									indexremove[bb] = aa;
								} else {
									isremove[aa] = true;
									indexremove[bb] = aa;
									// break;
								}
							}
						}
					}
				}
			}
		} else { // is RC
			for (size_t i1 = beg; i1 < end; ++i1) {
				aa = ridvec[i1]>>32;
				if (!isremove[aa]) {
					strcpy(stra, seq[aa].seq);
					if (ridvec[i1]&1) {
						reverseComplement(stra, L);
					}

					for (size_t i2 = i1 + 1; i2 < end; ++i2) {
						bb = ridvec[i2]>>32;
						if (!isremove[bb]) {
							strcpy(strb, seq[bb].seq);
							if (ridvec[i2]&1) {
								reverseComplement(strb, L);
							}
							difnum = 0;
							for (int ss = 0; ss < L; ++ss) {
								if (stra[ss] != strb[ss]) ++difnum;
							}
							if (difnum <= difthr) {
								if (avgQual[aa] >= avgQual[bb]) {
									isremove[bb] = true;
									indexremove[bb] = aa;
								} else {
									isremove[aa] = true;
									indexremove[bb] = aa;
									// break;
								}
							}
						}
					}
				}
			}
		}
	}
}

void processClusterPEFun() {
	size_t beg, end;

	char *stra = (char*)alloca(((L<<1) + 2) * sizeof(char));
	char *strb = (char*)alloca(((L<<1) + 2) * sizeof(char));
	// cout << "processClusterPEFun()\n";
	uint32_t aa, bb;
	int difnum;

	while (1) {
		size_t i = __sync_fetch_and_add(&rid_pthread, 1);
		if (i >= sgv.n) break;

		beg = sgv.a[i].a;
		end = sgv.a[i].a + sgv.a[i].b;

		if (!isRC) {
			for (size_t i1 = beg; i1 < end; ++i1) {
				aa = ridvec[i1]>>32;
				if (!isremove[aa]) {
					strcpy(stra, seq1[aa].seq);
					strcpy(stra + L, seq2[aa].seq);

					for (size_t i2 = i1 + 1; i2 < end; ++i2) {
						bb = ridvec[i2]>>32;
						if (!isremove[bb]) {
							strcpy(strb, seq1[bb].seq);
							strcpy(strb + L, seq2[bb].seq);

							difnum = 0;
							for (int ss = 0; ss < (L<<1); ++ss) {
								if (stra[ss] != strb[ss]) ++difnum;
							}
							if (difnum <= difthr) {
								if (avgQual[aa] >= avgQual[bb]) {
									isremove[bb] = true;
									indexremove[bb] = aa;
								} else {
									isremove[aa] = true;
									indexremove[bb] = aa;
									// break; 
								}
							}
						}
					}
				}
			}
		} else { // is RC
			for (size_t i1 = beg; i1 < end; ++i1) {
				aa = ridvec[i1]>>32;
				if (!isremove[aa]) {
					strcpy(stra, seq1[aa].seq);
					strcpy(stra + L, seq2[aa].seq);
					if (ridvec[i1]&1) {
						// reverseComplement(stra, 2*L);
						reverseComplementPE(stra, L<<1);
					}

					for (size_t i2 = i1 + 1; i2 < end; ++i2) {
						bb = ridvec[i2]>>32;
						if (!isremove[bb]) {
							strcpy(strb, seq1[bb].seq);
							strcpy(strb + L, seq2[bb].seq);
							if (ridvec[i2]&1) {
								// reverseComplement(strb, 2*L);
								reverseComplementPE(strb, L<<1);
							}
							difnum = 0;
							for (int ss = 0; ss < (L<<1); ++ss) {
								if (stra[ss] != strb[ss]) ++difnum;
							}
							if (difnum <= difthr) {
								if (avgQual[aa] >= avgQual[bb]) {
									isremove[bb] = true;
									indexremove[bb] = aa;
								} else {
									isremove[aa] = true;
									indexremove[bb] = aa;
									// break; 
								}
							}
						}
					}
				}
			}
		}
	}
}

inline void processCluster() {
	// rid_pthread = 0;
	rid_pthread = sgvnl;
	std::vector<thread> threadVec;
	for (int i = 0; i < nthreads; ++i) {
		if (isPE) threadVec.push_back(std::thread(processClusterPEFun));
		else threadVec.push_back(std::thread(processClusterSEFun));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});
	threadVec.clear();
}

size_t lend, cgnum;

void processLargeClusterSEFun() {
	char *stra = (char*)alloca((L + 1) * sizeof(char));
	char *strb = (char*)alloca((L + 1) * sizeof(char));
	// cout << "text xxx\n";
	uint32_t aa, bb;
	int difnum;

	while (1) {
		size_t i1 = __sync_fetch_and_add(&rid_pthread, 1);
		if (i1 >= lend) break;

		if (!isRC) {
			aa = ridvec[i1]>>32;
			if (!isremove[aa]) {
				strcpy(stra, seq[aa].seq);

				for (size_t i2 = i1 + 1; i2 < lend; ++i2) {
					bb = ridvec[i2]>>32;
					if (!isremove[bb]) {
						strcpy(strb, seq[bb].seq);
						difnum = 0;
						for (int ss = 0; ss < L; ++ss) {
							if (stra[ss] != strb[ss]) ++difnum;
						}
						if (difnum <= difthr) {
							if (avgQual[aa] >= avgQual[bb]) {
								isremove[bb] = true;
								indexremove[bb] = aa;
							} else {
								isremove[aa] = true;
 								indexremove[bb] = aa;                               
								// break;
							}
						}
					}
					if (cgnum > 40000 && i2 - i1 > 5000) {
						break;
					}
				}
			}
		} else { // is RC
			aa = ridvec[i1]>>32;
			if (!isremove[aa]) {
				strcpy(stra, seq[aa].seq);
				if (ridvec[i1]&1) {
					reverseComplement(stra, L);
				}

				for (size_t i2 = i1 + 1; i2 < lend; ++i2) {
					bb = ridvec[i2]>>32;
					if (!isremove[bb]) {
						strcpy(strb, seq[bb].seq);
						if (ridvec[i2]&1) {
							reverseComplement(strb, L);
						}
						difnum = 0;
						for (int ss = 0; ss < L; ++ss) {
							if (stra[ss] != strb[ss]) ++difnum;
						}
						if (difnum <= difthr) {
							if (avgQual[aa] >= avgQual[bb]) {
								isremove[bb] = true;
								indexremove[bb] = aa;
							} else {
								isremove[aa] = true;
								indexremove[bb] = aa;
								// break;
							}
						}
					}
					if (cgnum > 40000 && i2 - i1 > 5000) {
						break;
					}
				}
			}
		}
	}
}

void processLargeClusterPEFun() {
	size_t beg, end;

	char *stra = (char*)alloca(((L<<1) + 2) * sizeof(char));
	char *strb = (char*)alloca(((L<<1) + 2) * sizeof(char));
	// cout << "processClusterPEFun()\n";
	uint32_t aa, bb;
	int difnum;

	while (1) {
		size_t i1 = __sync_fetch_and_add(&rid_pthread, 1);
		if (i1 >= lend) break;

		if (!isRC) {
			aa = ridvec[i1]>>32;
			if (!isremove[aa]) {
				strcpy(stra, seq1[aa].seq);
				strcpy(stra + L, seq2[aa].seq);

				for (size_t i2 = i1 + 1; i2 < end; ++i2) {
					bb = ridvec[i2]>>32;
					if (!isremove[bb]) {
						strcpy(strb, seq1[bb].seq);
						strcpy(strb + L, seq2[bb].seq);
						difnum = 0;
						for (int ss = 0; ss < (L<<1); ++ss) {
							if (stra[ss] != strb[ss]) ++difnum;
						}
						if (difnum <= difthr) {
							if (avgQual[aa] >= avgQual[bb]) {
								isremove[bb] = true;
								indexremove[bb] = aa;
							} else {
								isremove[aa] = true;
								indexremove[bb] = aa;
								// break; 
							}
						}
					}
					if (cgnum > 40000 && i2 - i1 > 5000) {
						break;
					}
				}
			}
		} else { // is RC
			aa = ridvec[i1]>>32;
			if (!isremove[aa]) {
				strcpy(stra, seq1[aa].seq);
				strcpy(stra + L, seq2[aa].seq);
				if (ridvec[i1]&1) {
					reverseComplementPE(stra, L<<1);
				}

				for (size_t i2 = i1 + 1; i2 < end; ++i2) {
					bb = ridvec[i2]>>32;
					if (!isremove[bb]) {
						strcpy(strb, seq1[bb].seq);
						strcpy(strb + L, seq2[bb].seq);
						if (ridvec[i2]&1) {
							reverseComplementPE(stra, L<<1);
						}
						difnum = 0;
						for (int ss = 0; ss < L; ++ss) {
							if (stra[ss] != strb[ss]) ++difnum;
						}
						if (difnum <= difthr) {
							if (avgQual[aa] >= avgQual[bb]) {
								isremove[bb] = true;
								indexremove[bb] = aa;
							} else {
								isremove[aa] = true;
								indexremove[bb] = aa;
								// break; 
							}
						}
					}
					if (cgnum > 40000 && i2 - i1 > 5000) {
						break;
					}
				}
			}
		}
	}
}

inline void processLargeCluster() {
	std::vector<thread> threadVec;

	// for (size_t i = 0; i < lsgv.n; ++i) {
	// 	rid_pthread = lsgv.a[i].a;
	// 	lend = lsgv.a[i].a + lsgv.a[i].b;
	for (size_t i = 0; i < sgvnl; ++i) {
		rid_pthread = sgv.a[i].a;
		lend = sgv.a[i].a + sgv.a[i].b;
		cgnum = sgv.a[i].b;

		for (int j = 0; j < nthreads; ++j) {
			if (isPE) threadVec.push_back(std::thread(processLargeClusterPEFun));
			else threadVec.push_back(std::thread(processLargeClusterSEFun));
		}
		std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
			thr.join();
		});
		threadVec.clear();
	}
}

inline void removeDuplicate() {
	// CStopWatch rdt;
	// rdt.start();

	// stopwatch.resume();
	seq = seq1;

	// cout << seq1[0].seq << endl;
	// cout << seq2[0].seq << endl;
	// if (isPE && (kmer^1)) {
	// 	seq = seq2;
	// }
	calcMinimizers();
	// cout << "after calcMinimizers()\n";
	// cout << "Time of calcMinimizers() = " << stopwatch.stop() << std::endl;
	// stopwatch.resume();
	// exit(0);

	for (int kmervecidx = 0; kmervecidx < kmervecsize; ++kmervecidx) {
		if (kmervecidx > 0 && kmervec[kmervecidx - 1] == kmervec[kmervecidx])
			continue;

		for (int i = 0; i < (1 << bsize); ++i) {
			kv_init(B[i]);
		}
		if (isPE) {
			if (kmervecidx&1) {
				seq = seq2;
			} else {
				seq = seq1;
			}
		}
		sortBuckets(kmervecidx);
		// cout << "Time of sortBuckets() = " << stopwatch.stop() << std::endl;
		// stopwatch.resume();
		
		processCluster();

		processLargeCluster();
		kv_destroy(sgv);
		// kv_destroy(lsgv);

		for (int i = 0; i < (1 << bsize); ++i) {
			kv_destroy(B[i]);
		}
		// cout << "Time of processBuckets() = " << stopwatch.stop() << std::endl;
		// stopwatch.resume();
	}

	// cout << "Time of removeDuplicate() = " << rdt.stop() << std::endl;
	// cout << "--------------\n";
}

inline void removeLongDuplicate() {
	// CStopWatch rdt;
	// rdt.start();

	// stopwatch.resume();
	seq = seq1;

	calcLongMinimizers();
	// cout << "Time of calcLongMinimizers() = " << stopwatch.stop() << std::endl;
	// stopwatch.resume();
	// exit(0);

	for (int kmervecidx = 0; kmervecidx < lkmervecsize; ++ kmervecidx) {
		if (kmervecidx > 0 && lkmervec[kmervecidx - 1] == lkmervec[kmervecidx])
			continue;

		// cout << "kmer: " << lkmervec[kmervecidx] << endl;

		for (int i = 0; i < (1 << bsize); ++i) {
			kv_init(BL[i]);
		}
		
		sortLongBuckets(kmervecidx);
		// cout << "Time of sortLongBuckets() = " << stopwatch.stop() << std::endl;
		// stopwatch.resume();
		
		// cout << "before processCluster\n";
		processCluster();

		// cout << "before processLargeCluster\n";
		processLargeCluster();
		// cout << "before kv_destroy(sgv)\n";
		kv_destroy(sgv);
		// kv_destroy(lsgv);

		// cout << "before kv_destroy(BL[i]);\n";
		for (int i = 0; i < (1 << bsize); ++i) {
			kv_destroy(BL[i]);
		}
		// cout << "Time of processBuckets() = " << stopwatch.stop() << std::endl;
		// stopwatch.resume();
	}
	// cout << "over removeLongDuplicate()\n";
	// cout << "Time of removeLongDuplicate() = " << rdt.stop() << std::endl;
	// cout << "--------------\n";
}

int maxStrlen()
{
	int supposed_max = 0;
	int curr_len = 0;

	for (int rid = 0; rid < max_rid; rid++)
	{
		curr_len = strlen(seq1[rid].seq);
		if (curr_len > supposed_max)
		{
			supposed_max = curr_len;
		}
	}

	return supposed_max;
}


int main(int argc, char* argv[]) {
	// stopwatch.start();
	init();
	// L = getReadsLength(argv[1]);
	getPars(argc, argv);

	bool rb = getReads(rf1.c_str());

	if (!rb) {
		fprintf(stderr, "Read the file %s fail.\n", rf1.c_str());
		return 0;
	}

	if (isPE) {
		rb = getReads2(rf2.c_str());
		if (!rb) {
			fprintf(stderr, "Read the file %s fail.\n", rf2.c_str());
			return 0;
		}
	}
	L = maxStrlen();
	if (iskf) {
		kmervecsize = 0;
		kmervec = new int[32];
		lkmervecsize = 0;
		lkmervec = new int[32];

		int tkmer;
		FILE *fp = fopen(kf.c_str(), "r");
		while (fscanf(fp, "%d", &tkmer) != EOF) {
			if (tkmer > 62 && tkmer < 13) {
				fprintf(stderr, "kmer %d is ignored\n", tkmer);
			} else {
				if (tkmer < 32) {
					kmervec[kmervecsize] = tkmer;
					++ kmervecsize;
				} else {
					lkmervec[lkmervecsize] = tkmer;
					++ lkmervecsize;
				}
			}
		}
		fclose(fp);

		sort (kmervec, kmervec + kmervecsize, greater<int>() );
		sort (lkmervec, lkmervec + lkmervecsize, greater<int>() );
	} else {
		if (L < 105) {		
			if (difthr == 0) {
				kmervecsize = 1;
			} else 
			if (difthr == 1) {
				kmervecsize = 9;
			} else 
			if (difthr == 2) {
				kmervecsize = 10;
			} else {
				kmervecsize = 12;
			}
			lkmervecsize = 0;
			if (isRC) kmervecsize = 16;
			if (difthr == 0) kmervecsize = 1;

			if (kmervecsize > 0)  {
				kmervec = new int[kmervecsize];
				if (kmervecsize == 1) {
					kmervec[0] = 29;
					if (isPE) {
						lkmervecsize = 1;
						lkmervec = new int[lkmervecsize];
						lkmervec[0] = 41;
						kmervecsize = 0;
					}
				} else {
					kmervec[0] = 28;
					if (isPE) {
						lkmervecsize = 1;
						lkmervec = new int[lkmervecsize];
						lkmervec[0] = 41;
					}
				}

				if (isRC) kmervec[0] = 31;

				for (int i = 1; i < kmervecsize; ++i) {
					kmervec[i] = kmervec[i-1] - 1;
				}
				// if (lkmervecsize == 0) max_kmer = kmervec[0];
			}

			if (max_rid > (1ULL<<27)) { //134217728 214364631
				// 51 49 48 47 46 32 
				// 31 30 29 28
				lkmervecsize = 6;
				if (isRC) {
					lkmervecsize = 8;
					if (difthr >= 3) lkmervecsize = 9;
				}
				lkmervec = new int[lkmervecsize];
				lkmervec[0] = 51;
				lkmervec[1] = 49;
				lkmervec[2] = 48;
				lkmervec[3] = 47;
				lkmervec[4] = 46;
				if(isRC) {
					lkmervec[5] = 45;			
					lkmervec[6] = 44;
					if (difthr >= 3) lkmervec[7] = 33;		
				} 
				lkmervec[lkmervecsize - 1] = 32;

				kmervecsize = 4;

				kmervec = new int[kmervecsize];
				kmervec[0] = 31;
				kmervec[1] = 30;
				kmervec[2] = 29;
				kmervec[3] = 28;

				if (difthr == 0) {
					lkmervecsize = 1;
					kmervecsize = 0;
				}
			}
		} else 
		if (L >= 105) {
			if (difthr == 0) {
				lkmervecsize = 1;
			} else 
			if (difthr == 1) {
				lkmervecsize = 9;
			} else 
			if (difthr == 2) {
				lkmervecsize = 10;
			} else {
				lkmervecsize = 12;
			}
			if (isRC && difthr > 0) {
				lkmervecsize += 2;
				if (difthr >= 3) ++ lkmervecsize;
			}

			lkmervec = new int[lkmervecsize];
			// 
			kmervecsize = 4;
			if (isRC) ++ kmervecsize;

			kmervec = new int[kmervecsize];
			kmervec[0] = 31;
			kmervec[1] = 30;
			kmervec[2] = 29;
			kmervec[3] = 28;
			if (isRC) kmervec[4] = 27;

			if (difthr == 0) {
				kmervecsize = 1;
			}

			if (L < 120) {
				lkmervec[0] = 41;
			} else 
			if (L >= 120 && L < 140) {
				lkmervec[0] = 45;
			} else 
			if (L >= 140 && L < 153) {
				lkmervec[0] = 51;
			} else {
				lkmervec[0] = 61;
			}
			if (isRC && difthr > 0) ++ lkmervec[0];

			for (int i = 1; i < lkmervecsize; ++i) {
				lkmervec[i] = lkmervec[i-1] - 1;
			}
		}

	}
	// cout << "isRC: " << isRC << endl;
	// lkmervecsize = 0; // for test
	fprintf(stderr, "large kmer:\n");
	cout << "lkmervecsize: " << lkmervecsize << endl;
	for (int i = 0; i < lkmervecsize; ++i) {
		fprintf(stderr, "%d: %d\n", i, lkmervec[i]);
	}

	fprintf(stderr, "---\n");
	fprintf(stderr, "small kmer:\n");
	cout << "kmervecsize: " << kmervecsize << endl;
	for (int i = 0; i < kmervecsize; ++i) {
		fprintf(stderr, "%d: %d\n", i, kmervec[i]);
	}
	// exit (0);
	// fprintf(stderr, "%lu: %s\n", did, seq[did].seq);

	cout << "L: " << L << endl;
	cout << "max_rid: " << max_rid << endl;
	// cout << "Time of read file = " << stopwatch.stop() << std::endl;
	// stopwatch.resume();

	// char str[1<<10];
	// cout << "pls input a string: ";
	// scanf("%s", str);
	isremove = new bool[max_rid];
	memset(isremove, false, sizeof(bool)*max_rid);
	indexremove = new int64_t[max_rid];
	memset(indexremove, -1, sizeof(int64_t) * max_rid);
	avgQual = new int[max_rid];
	memset(avgQual, 0, sizeof(int)*max_rid);

	ridvec = new uint64_t[max_rid];

	if (lkmervecsize > 0) {
		BL = (mm192_v*)calloc(1 << bsize, sizeof(mm192_v));
		BLsize = new size_t[(1 << bsize) * lkmervecsize];
		memset(BLsize, 0, sizeof(size_t)*(1 << bsize) * lkmervecsize);

		removeLongDuplicate();
		free(BL);
		delete[] lminis;	
		delete[] BLsize;	
	}

	if (kmervecsize > 0)  {
		B = (mm128_v*)calloc(1 << bsize, sizeof(mm128_v));
		Bsize = new size_t[(1 << bsize) * kmervecsize];
		memset(Bsize, 0, sizeof(size_t)*(1 << bsize) * kmervecsize);

		removeDuplicate();
		free(B);
		delete[] minis;
		delete[] Bsize;
	}

	delete[] ridvec;
	delete[] lkmervec;
	delete[] kmervec;    
	
	// stopwatch.resume();

	int64_t curr_parent = -1;
	std::string *parents_buffer = new std::string[max_rid];

	if (!isPE) {
		seq = seq1;
		// FILE *fpo = fopen("rmdmini.txt", "w");
		// FILE *fp = fopen("result.fastq", "w");
		uint32_t removed = 0;
		// FILE *fp = fopen(rsf.c_str(), "w");
		ofstream fsout;
		fsout.open(rsf);

		string buf;
		buf.reserve(BUFFER_SIZE);

		for (size_t rid = 0; rid < max_rid; ++rid) {
			// kseq_destroy
			if (!isremove[rid]) {
				// using a buffer to speed
				// fprintf(fp, "%s\n%s\n", seq[rid].name, seq[rid].seq);
				buf.append(header_prefix + string(seq[rid].name) + "\n" + string(seq[rid].seq) + "\n");
				if (seq[rid].qual != NULL) {
					// fprintf(fp, "+\n%s\n", seq[rid].qual);
					buf.append("+\n" + string(seq[rid].qual) + "\n");
				}
				if (buf.size() > (size_t(BUFFER_SIZE*0.95))) {
					fsout << buf; 
					buf.clear();
				}
			} // else {++removed;}
			// else {
			// 	fprintf(fpo, "%s\n", seq[rid].seq);
			// }

            curr_parent = indexremove[rid];
			if (curr_parent != -1)
			{
				if (parents_buffer[curr_parent].length() == 0)
				{
					parents_buffer[curr_parent].append(header_prefix + string(seq[curr_parent].name) + ".\n" + string(seq[curr_parent].seq) + "\n");
				}

				parents_buffer[curr_parent].append(header_prefix + string(seq[rid].name) + "\n" + string(seq[rid].seq) + "\n");
			}
		}
		if (buf.size() > 0) {
			fsout << buf; 
			buf.clear();
		}
		fsout.close();
		// fclose(fp);
		// cout << "removed: " << removed << endl;
		// fclose(fpo);
	} else {
		// FILE *fp1 = fopen("result1.fastq", "w");
		// FILE *fp2 = fopen("result2.fastq", "w");
		// FILE *fp1 = fopen((rsf + "_1").c_str(), "w");
		// FILE *fp2 = fopen((rsf + "_2").c_str(), "w");

		ofstream fsout0;
		fsout0.open(rsf + "_1");		
		ofstream fsout1;
		fsout1.open(rsf + "_2");

		string buf0, buf1;
		buf0.reserve(BUFFER_SIZE >> 1);
		buf1.reserve(BUFFER_SIZE >> 1);

		for (size_t rid = 0; rid < max_rid; ++rid) {
			// kseq_destroy
			if (!isremove[rid]) {
				buf0.append(header_prefix + string(seq1[rid].name) + "\n" + string(seq1[rid].seq) + "\n");
				buf1.append(header_prefix + string(seq2[rid].name) + "\n" + string(seq2[rid].seq) + "\n");
				// fprintf(fp1, "%s\n%s\n", seq1[rid].name, seq1[rid].seq);
				// fprintf(fp2, "%s\n%s\n", seq2[rid].name, seq2[rid].seq);
				if (seq1[rid].qual != NULL) {
					buf0.append("+\n" + string(seq1[rid].qual) + "\n");
					// fprintf(fp1, "+\n%s\n", seq1[rid].qual);
				}
				if (seq2[rid].qual != NULL) {
					buf1.append("+\n" + string(seq2[rid].qual) + "\n");
					// fprintf(fp2, "+\n%s\n", seq2[rid].qual);
				}

				if (buf0.size() > (size_t((BUFFER_SIZE >> 1)*0.95))) {
					fsout0 << buf0; 
					buf0.clear();

					fsout1 << buf1; 
					buf1.clear();
				}
			}

            curr_parent = indexremove[rid];
			if (curr_parent != -1)
			{
				if (parents_buffer[curr_parent].length() == 0)
				{
					parents_buffer[curr_parent].append(header_prefix + string(seq[curr_parent].name) + ".\n" + string(seq[curr_parent].seq) + "\n");
				}

				parents_buffer[curr_parent].append(header_prefix + string(seq[rid].name) + "\n" + string(seq[rid].seq) + "\n");
			}
		}
		if (buf0.size() > 0) {
			fsout0 << buf0; 
			buf0.clear();

			fsout1 << buf1; 
			buf1.clear();
		}

		

		fsout0.close();
		fsout1.close();
		// fclose(fp1);
		// fclose(fp2);
	}
	// cout << "Time of saving file = " << stopwatch.stop() << std::endl;
	// delete[] seq;

	FILE *pFileDups;
	pFileDups = fopen(logpath.c_str(), "w");		

    for (int rid = 0; rid < max_rid; rid++)
	{
		if (parents_buffer[rid].length() > 0)
		{
			fprintf(pFileDups, "%s", parents_buffer[rid].c_str());
		}
	}

	fclose(pFileDups);

	return 0;
}

////////////

unsigned char seq_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	// A  B  C   D  E  F  G   H  I  J  K   L  M  N  O
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	// A  B  C   D  E  F  G   H  I  J  K   L  M  N  O
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

#include "ksort.h"
#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 
KSORT_INIT_GENERIC(uint32_t)
