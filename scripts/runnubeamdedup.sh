nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i SRR921889.fastq' </dev/null 1> SRR921889.out 2> SRR921889.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i SRR921889.fastq -o SRR921889.uniq.s.fastq -s 0' </dev/null 1> SRR921889_s.out 2> SRR921889_s.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i1 SRR948355_1.fastq -i2 SRR948355_2.fastq' </dev/null 1> SRR948355.out 2> SRR948355.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i1 SRR948355_1.fastq -i2 SRR948355_2.fastq -o SRR948355.uniq.s.fastq -s 0' </dev/null 1> SRR948355_s.out 2> SRR948355_s.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i SRR377645.fastq' </dev/null 1> SRR377645.out 2> SRR377645.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i SRR377645.fastq -o SRR377645.uniq.s.fastq -s 0' </dev/null 1> SRR377645_s.out 2> SRR377645_s.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i simulated.fastq' </dev/null 1> simulated.out 2> simulated.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i simulated.fastq -o simulated.uniq.s.fastq -s 0' </dev/null 1> simulated_s.out 2> simulated_s.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i1 denisova_1.fastq -i2 denisova_2.fastq' </dev/null 1> denisova.out 2> denisova.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i1 denisova_1.fastq -i2 denisova_2.fastq -o1 denisova_1.uniq.s.fastq -o2 denisova_2.uniq.s.fastq -s 0' </dev/null 1> denisova_s.out 2> denisova_s.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i1 SRR12175235_1.fastq -i2 SRR12175235_2.fastq' </dev/null 1> SRR12175235.out 2> SRR12175235.err
nohup bash -c '/usr/bin/time -v ./nubeam-dedup -i1 SRR12175235_1.fastq -i2 SRR12175235_2.fastq -o1 SRR12175235_1.uniq.s.fastq -o2 SRR12175235_2.uniq.s.fastq -s 0' </dev/null 1> SRR12175235_s.out 2> SRR12175235_s.err
