from glob import glob
import sh
projs = glob('Projects/*')
import os
for pdir in projs:
    # pdir = glob( proj + '/*runsample*')[0] 
    # pdir = glob(proj + '/*runsample*')[0]
    cons = glob(pdir + '/*.consensus.fasta')[0]
    # print "remove: ", cons, os.remove)
    # os.remove(cons)
    bam = glob(pdir + '/' + '*.bam')[0]
    vcf = glob(pdir + '/' + '*.vcf')[0]
    ref = [f for f in glob(pdir + '/' + '*.fasta') if (f != cons)][0]

    # print dict(ref=ref, vcf=vcf, bam=bam, majority=95, minbq=30, mind=10, output=cons)
    # sh.lf_consensus(ref=ref, vcf=vcf, bam=bam, majority=95, minbq=30, mind=10, output=cons)
    print ('with options',  dict(vcf=vcf, ref=ref, bam=bam, type='base_caller', minpercent=1, mindepth=1000, out=vcf + '.vartable.tsv') )
    sh.vartable_report(vcf, ref=ref, bam=bam, type='base_caller', minpercent=1, mindepth=1000, out=vcf + '.vartable.tsv')
    # lf_consensus(ref=ref, vcf=vcf, bam=bam, majority=95, minbq=30, mind=10, output=cons)

# usage: vartable_report [-h] --type {lofreq,base_caller} --out OUT --mindepth
#                        MIND --minpercent MINP [--bam BAM] [--ref REF]
#                        [--genbank GENBANK] [--cds-rev start end]
#                        [--cds-fwd start end]
#                        vcf_path
# vartable_report: error: the following arguments are required: vcf_path, --type, --out, --mindepth, --minpercent


