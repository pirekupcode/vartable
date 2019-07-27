test/bin/run-linters.sh || true
nosetests test/*.py
vartable_report testdata/fullsample.bam.vcf --bam testdata/fullsample.bam --ref testdata/Den1__WestPac__1997.fasta --type base_caller  --mindepth 10 --minpercent 1 --out testdata/docker-test.tsv
python -c "assert (open('testdata/docker-test.tsv').read() == open('testdata/example.tsv').read())"
