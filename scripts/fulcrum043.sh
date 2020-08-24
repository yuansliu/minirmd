#!/bin/bash
set -e

if [[ $1 == "-i" ]]; then
	while getopts ":i:c:b:n:" opt; do
		case "$opt" in
			i) filename=$OPTARG;;
			c) c=$OPTARG;;
			b) b=$OPTARG;;
			n) n=$OPTARG;;
		esac
	done
	# echo $filename // for test
	inf=$(basename "$filename" .fastq)_rm_"$c"_b"$b"_in
	outf=$(basename "$filename" .fastq)_rm_"$c"_b"$b"_out
	mkdir -p $inf
	mkdir -p $outf
	python linenumfast05.py -b $b -t s -i $filename -o $inf/
	python parallel05.py -t s -n $n -q 0 -i $inf/ -o $outf/ -c $c
	rm -rf $inf 
	# $outf
elif [[ $1 == "-1" ]]; then
	while getopts ":1:2:c:b:n:" opt; do
		case "$opt" in
			1) f1=$OPTARG;;
			2) f2=$OPTARG;;
			c) c=$OPTARG;;
			b) b=$OPTARG;;
			n) n=$OPTARG;;
		esac
	done
	pef=`echo $(basename "$f1" .fastq) | cut -d \_ -f 1`_pe.fulcrum
	inf=`echo $(basename "$f1" .fastq) | cut -d \_ -f 1`_pe_rm_"$c"_b"$b"_in
	outf=`echo $(basename "$f1" .fastq) | cut -d \_ -f 1`_pe_rm_"$c"_b"$b"_out
	mkdir -p $inf
	mkdir -p $outf
	python merger.py $f1 $f2 $pef
	python linenumfast05.py -b $b -t p -i $pef -o $inf/
	python parallel05.py -t p -n $n -q 0 -i $inf/ -o $outf/ -c $c
	rm -rf $inf 
	# $outf
fi

