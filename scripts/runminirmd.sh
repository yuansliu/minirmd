nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_0_minirmd.fastq -d 0' </dev/null 1> SRR921889_rm_0_minirmd.out 2> SRR921889_rm_0_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_1_minirmd.fastq -d 1' </dev/null 1> SRR921889_rm_1_minirmd.out 2> SRR921889_rm_1_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_2_minirmd.fastq -d 2' </dev/null 1> SRR921889_rm_2_minirmd.out 2> SRR921889_rm_2_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_3_minirmd.fastq -d 3' </dev/null 1> SRR921889_rm_3_minirmd.out 2> SRR921889_rm_3_minirmd.err

nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_0_minirmd.fastq -d 0' </dev/null 1> SRR948355_pe_rm_0_minirmd.out 2> SRR948355_pe_rm_0_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_1_minirmd.fastq -d 1' </dev/null 1> SRR948355_pe_rm_1_minirmd.out 2> SRR948355_pe_rm_1_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_2_minirmd.fastq -d 2' </dev/null 1> SRR948355_pe_rm_2_minirmd.out 2> SRR948355_pe_rm_2_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_3_minirmd.fastq -d 3' </dev/null 1> SRR948355_pe_rm_3_minirmd.out 2> SRR948355_pe_rm_3_minirmd.err

nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_0_minirmd.fastq -d 0' </dev/null 1> SRR377645_rm_0_minirmd.out 2> SRR377645_rm_0_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_1_minirmd.fastq -d 1' </dev/null 1> SRR377645_rm_1_minirmd.out 2> SRR377645_rm_1_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_2_minirmd.fastq -d 2' </dev/null 1> SRR377645_rm_2_minirmd.out 2> SRR377645_rm_2_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_3_minirmd.fastq -d 3' </dev/null 1> SRR377645_rm_3_minirmd.out 2> SRR377645_rm_3_minirmd.err

# if removing duplicates on RC strand
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_0_minirmd.rc.fastq -d 0 -r' </dev/null 1> SRR921889_rm_0_minirmd.rc.out 2> SRR921889_rm_0_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_1_minirmd.rc.fastq -d 1 -r' </dev/null 1> SRR921889_rm_1_minirmd.rc.out 2> SRR921889_rm_1_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_2_minirmd.rc.fastq -d 2 -r' </dev/null 1> SRR921889_rm_2_minirmd.rc.out 2> SRR921889_rm_2_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR921889.fastq -o SRR921889_rm_3_minirmd.rc.fastq -d 3 -r' </dev/null 1> SRR921889_rm_3_minirmd.rc.out 2> SRR921889_rm_3_minirmd.rc.err

nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_0_minirmd.rc.fastq -d 0 -r' </dev/null 1> SRR948355_pe_rm_0_minirmd.rc.out 2> SRR948355_pe_rm_0_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_1_minirmd.rc.fastq -d 1 -r' </dev/null 1> SRR948355_pe_rm_1_minirmd.rc.out 2> SRR948355_pe_rm_1_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_2_minirmd.rc.fastq -d 2 -r' </dev/null 1> SRR948355_pe_rm_2_minirmd.rc.out 2> SRR948355_pe_rm_2_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR948355_1.fastq -f SRR948355_2.fastq -o SRR948355_pe_rm_3_minirmd.rc.fastq -d 3 -r' </dev/null 1> SRR948355_pe_rm_3_minirmd.rc.out 2> SRR948355_pe_rm_3_minirmd.rc.err

nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_0_minirmd.rc.fastq -d 0 -r' </dev/null 1> SRR377645_rm_0_minirmd.rc.out 2> SRR377645_rm_0_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_1_minirmd.rc.fastq -d 1 -r' </dev/null 1> SRR377645_rm_1_minirmd.rc.out 2> SRR377645_rm_1_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_2_minirmd.rc.fastq -d 2 -r' </dev/null 1> SRR377645_rm_2_minirmd.rc.out 2> SRR377645_rm_2_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR377645.fastq -o SRR377645_rm_3_minirmd.rc.fastq -d 3 -r' </dev/null 1> SRR377645_rm_3_minirmd.rc.out 2> SRR377645_rm_3_minirmd.rc.err

nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_0_minirmd.fastq -d 0' </dev/null 1> simulated_rm_0_minirmd.out 2> simulated_rm_0_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_1_minirmd.fastq -d 1' </dev/null 1> simulated_rm_1_minirmd.out 2> simulated_rm_1_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_2_minirmd.fastq -d 2' </dev/null 1> simulated_rm_2_minirmd.out 2> simulated_rm_2_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_3_minirmd.fastq -d 3' </dev/null 1> simulated_rm_3_minirmd.out 2> simulated_rm_3_minirmd.err

nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_0_minirmd.fastq -d 0' </dev/null 1> denisova_pe_rm_0_minirmd.out 2> denisova_pe_rm_0_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_1_minirmd.fastq -d 1' </dev/null 1> denisova_pe_rm_1_minirmd.out 2> denisova_pe_rm_1_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_2_minirmd.fastq -d 2' </dev/null 1> denisova_pe_rm_2_minirmd.out 2> denisova_pe_rm_2_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_3_minirmd.fastq -d 3' </dev/null 1> denisova_pe_rm_3_minirmd.out 2> denisova_pe_rm_3_minirmd.err

nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_0_minirmd.fastq -d 0' </dev/null 1> SRR12175235_pe_rm_0_minirmd.out 2> SRR12175235_pe_rm_0_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_1_minirmd.fastq -d 1' </dev/null 1> SRR12175235_pe_rm_1_minirmd.out 2> SRR12175235_pe_rm_1_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_2_minirmd.fastq -d 2' </dev/null 1> SRR12175235_pe_rm_2_minirmd.out 2> SRR12175235_pe_rm_2_minirmd.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_3_minirmd.fastq -d 3' </dev/null 1> SRR12175235_pe_rm_3_minirmd.out 2> SRR12175235_pe_rm_3_minirmd.err

nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_0_minirmd.rc.fastq -d 0 -r' </dev/null 1> simulated_rm_0_minirmd.rc.out 2> simulated_rm_0_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_1_minirmd.rc.fastq -d 1 -r' </dev/null 1> simulated_rm_1_minirmd.rc.out 2> simulated_rm_1_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_2_minirmd.rc.fastq -d 2 -r' </dev/null 1> simulated_rm_2_minirmd.rc.out 2> simulated_rm_2_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i simulated.fastq -o simulated_rm_3_minirmd.rc.fastq -d 3 -r' </dev/null 1> simulated_rm_3_minirmd.rc.out 2> simulated_rm_3_minirmd.rc.err

nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_0_minirmd.rc.fastq -d 0 -r' </dev/null 1> denisova_pe_rm_0_minirmd.rc.out 2> denisova_pe_rm_0_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_1_minirmd.rc.fastq -d 1 -r' </dev/null 1> denisova_pe_rm_1_minirmd.rc.out 2> denisova_pe_rm_1_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_2_minirmd.rc.fastq -d 2 -r' </dev/null 1> denisova_pe_rm_2_minirmd.rc.out 2> denisova_pe_rm_2_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i denisova_1.fastq -f denisova_2.fastq -o denisova_pe_rm_3_minirmd.rc.fastq -d 3 -r' </dev/null 1> denisova_pe_rm_3_minirmd.rc.out 2> denisova_pe_rm_3_minirmd.rc.err

nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_0_minirmd.rc.fastq -d 0 -r' </dev/null 1> SRR12175235_pe_rm_0_minirmd.rc.out 2> SRR12175235_pe_rm_0_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_1_minirmd.rc.fastq -d 1 -r' </dev/null 1> SRR12175235_pe_rm_1_minirmd.rc.out 2> SRR12175235_pe_rm_1_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_2_minirmd.rc.fastq -d 2 -r' </dev/null 1> SRR12175235_pe_rm_2_minirmd.rc.out 2> SRR12175235_pe_rm_2_minirmd.rc.err
nohup bash -c '/usr/bin/time -v ./minirmd -i SRR12175235_1.fastq -f SRR12175235_2.fastq -o SRR12175235_pe_rm_3_minirmd.rc.fastq -d 3 -r' </dev/null 1> SRR12175235_pe_rm_3_minirmd.rc.out 2> SRR12175235_pe_rm_3_minirmd.rc.err
