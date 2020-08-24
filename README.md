# minirmd

minirmd is a tool to remove duplicate reads. The program is written in C++11 (tested with g++ >= 5.4.0) and works on Linux.


## Download & Compile

	git clone https://github.com/yuansliu/bfmem.git
	cd bfmem
	make

## Usage

	./bfmem -r H.all.fa -q M.all.fa -o hm-100.txt [options]

Parameters

	-r  	reference genome, a multi-FASTA file
	-q  	query genome, a multi-FASTA file
	-o  	output file

Options

	-l  	minimal length of matches; default is 100
	-k  	length of k-mer (some default values are set)
	-t  	number of threads (default t=10)
	-s  	strands; 
			default is foward; 
			'r' is reverse-complement; 
			'b' is both foward and reverse-complement.
	-h  	print help message

### Format conversion & Result comparison 

Convert the output of bfMEM to the format of copMEM (i.e., MEMs with the same order).

	./formatconvert f1 f2

*f1* is the result of bfMEM; *f2* stores the format converted result.

	./compareMEM f1 f2

Compare the MEMs in the two files *f1* and *f2*. The program compares two files row by row.

## Citation
Yuansheng Liu, Leo Yu Zhang, Jinyan Li. Fast detection of maximal exact matches via fixed sampling of query k-mers and Bloom filtering of index k-mers. Bioinformatics, 2019.

## Contacts
If any bugs during you run our code, please email to <yyuanshengliu@gmail.com>
