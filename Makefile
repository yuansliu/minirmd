
CC = g++

# CPPFLAGS = -g -O3 -Wall -std=c++14
CPPFLAGS = -O3 -Wno-used-function -std=c++14 

# TBBLIBS = -lpthread -lz -ltbb
LIBS = -lpthread -lz #-lm #-ltbb

# BSCSRCS = libbsc/bsc.cpp libbsc/libbsc/adler32/adler32.cpp libbsc/libbsc/bwt/divsufsort/divsufsort.c libbsc/libbsc/bwt/bwt.cpp libbsc/libbsc/coder/coder.cpp libbsc/libbsc/coder/qlfc/qlfc.cpp libbsc/libbsc/coder/qlfc/qlfc_model.cpp libbsc/libbsc/filters/detectors.cpp libbsc/libbsc/filters/preprocessing.cpp libbsc/libbsc/libbsc/libbsc.cpp libbsc/libbsc/lzp/lzp.cpp libbsc/libbsc/platform/platform.cpp 

# SRCS = read.cpp hash.cpp compress.cpp fqreader.cpp decompress.cpp main.cpp
SRCS = src/minirmd.cpp src/bseq.c # src/sketch.c
# SRCS = src/minirmd.c src/bseq.c src/sketch.c

# GROUNDSRCS = gsrc/minirmd.c gsrc/bseq.c

OBJS = $(SRCS: .c = .o)

EXEC = minirmd

$(EXEC) : $(OBJS)
	$(CC) $(CPPFLAGS) $^ $(LIBS) -o $@

# ground: $(GROUNDSRCS)
# 	$(CC) $(CPPFLAGS) $(TBBLIBS) $^ -o $@

%.o : %.c
	$(CC) -c $(CPPFLAGS) $<

clean:
	rm -f *.o $(EXEC) *.out
