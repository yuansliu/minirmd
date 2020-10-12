# minirmd

minirmd is a tool to remove duplicate reads. The program is written in C++11 (tested with g++ >= 5.4.0) and works on Linux.


## Download & Compile

	git clone https://github.com/yuansliu/minirmd.git
	cd minirmd
	make


## Install from bioconda
	conda install -c bioconda minirmd

## Test
	sh runtest.sh
	
## Usage

	./minirmd -i test.fastq -o test_rm_1.fastq -d 1
	./minirmd -i test_1.fastq -f test_2.fastq -o test_rm_2.fastq -d 2

Parameters

	-i  	reads file, a FASTQ file
	-f  	reads file, a FASTQ file; if paired-end
	-o  	the output file

Options

	-d  	number of allowed mismatch; default is 0
	-k  	the file to sotre values of k
	-r  	remove duplicates on reverse-complement strand 
	-t  	number of threads (default 24)
	-h  	print help message

## Note
For the following fastq format file, the second "header123" will be ignored by the tool.
	
	@header123
	CGTATGCTACGTAC
	+header123
	KKKKKIIKKKKIII

## Contacts
If any bugs, please email to <yyuanshengliu@gmail.com>
