# VarTable

Uses bam-readcount and PyVCF to assist manual curation of viral genomes.

```bash

$ vartable_report testdata/fullsample.bam.vcf --bam testdata/fullsample.bam --ref testdata/Den1__WestPac__1997.fasta --type base_caller  --mindepth 10 --minpercent 1 --out testdata/example.tsv

```

Output Info
===========

"Ref Frequency" and "Alt Frequency" are percentage values.

"Codon", "Codon Type" and translation are currently unsupported.



####For bam-redcount:

Project home:   https://github.com/genome/bam-readcount

base → the base that all reads following in this field contain at the reported position i.e. C

count → the number of reads containing the base

avg_mapping_quality → the mean mapping quality of reads containing the base

avg_basequality → the mean base quality for these reads

avg_se_mapping_quality → mean single ended mapping quality

num_plus_strand → number of reads on the plus/forward strand

num_minus_strand → number of reads on the minus/reverse strand

avg_pos_as_fraction → average position on the read as a fraction (calculated with respect to the length after clipping). This value is normalized to the center of the read (bases occurring strictly at the center of the read have a value of 1, those occurring strictly at the ends should approach a value of 0)

avg_num_mismatches_as_fraction → average number of mismatches on these reads per base

avg_sum_mismatch_qualities → average sum of the base qualities of mismatches in the reads

num_q2_containing_reads → number of reads with q2 runs at the 3’ end

avg_distance_to_q2_start_in_q2_reads → average distance of position (as fraction of unclipped read length) to the start of the q2 run

avg_clipped_length → average clipped read length of reads

avg_distance_to_effective_3p_end → average distance to the 3’ prime end of the read (as fraction of unclipped read length)



Development
============

Docker
------

```bash
$ docker-compose -f docker-compose.test.yml -p ci build 
$ docker-compose -f docker-compose.test.yml -p ci up
```

Setup
-----

```bash
$ ./install.sh
$ export PATH=$PWD/miniconda/bin/:$PATH 
$ pip install -r test-requirements.txt
```

Make sure PATH points to miniconda/bin


Test 
-----

```bash
$ nosetests test/*.py --nocapture  # don't silence stdout
```

Test functions must start or end with `test` or they will be ignored by nosetets

Debugging 
-----

```bash
$ pip install ipdbplugin
$ nosetests test/*.py --ipdb --ipdb-failure --nocapture
```
