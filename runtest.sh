echo "Reads in the file testData/test.fastq:\n"
awk 'BEGIN {n=0} {if(n++ % 4 == 1) print $1}' testData/test.fastq
./minirmd -i testData/test.fastq -o testData/test_rm_0.fastq -d 0 1>/dev/null 2>/dev/null
echo '\nReads after removing duplicate:\n'
awk 'BEGIN {n=0} {if(n++ % 4 == 1) print $1}' testData/test_rm_0.fastq
echo "\nThe result is stored in the file testData/test_rm_0.fastq\n"