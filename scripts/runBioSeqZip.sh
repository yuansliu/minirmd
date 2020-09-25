nohup bash -c '/usr/bin/time -v ./bioseqzip_collapse -i SRR921889.fastq -f fastq -t 24' </dev/null 1> SRR921889.out 2> SRR921889.err
nohup bash -c '/usr/bin/time -v ./bioseqzip_collapse -i SRR948355_1.fastq -p SRR948355_2.fastq -f fastq -t 24' </dev/null 1> SRR948355.out 2> SRR948355.err
nohup bash -c '/usr/bin/time -v ./bioseqzip_collapse -i SRR377645.fastq -f fastq -t 24' </dev/null 1> SRR377645.out 2> SRR377645.err

nohup bash -c '/usr/bin/time -v ./bioseqzip_collapse -i simulated.fastq -f fastq -t 24' </dev/null 1> simulated.out 2> simulated.err
nohup bash -c '/usr/bin/time -v ./bioseqzip_collapse -i denisova_1.fastq -p denisova_2.fastq -f fastq -t 24' </dev/null 1> denisova.out 2> denisova.err
nohup bash -c '/usr/bin/time -v ./bioseqzip_collapse -i SRR12175235_1.fastq -p SRR12175235_2.fastq -f fastq -t 24' </dev/null 1> SRR12175235.out 2> SRR12175235.err
