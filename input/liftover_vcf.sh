#!/bin/bash

#
# LiftOver VCF file from hg19 to hg38
#

# INPUT:
VCFold="/home/stef/CLUSTER/scripts/B_mummy/testrun_mummy_vcf/input_vcf/PD_TARG2.filtered_snps_clean.vcf"
VCFnew="/home/stef/CLUSTER/scripts/B_mummy/testrun_mummy_vcf/input_vcf/PD_hg38.vcf"
CHAIN="/home/stef/CLUSTER/scripts/B_mummy/testrun_mummy_vcf/hg19ToHg38.over.chain"
REF="/home/stef/CLUSTER/reference/data/reference_seq_hg38/Homo_sapiens_assembly38.fasta.gz"


# no chain file to liftOver from b37 (which is original assembly to hg38), hence:
#grep "^#" $VCFold >header.txt
#grep -v "^#" $VCFold | sed 's/^/chr/g' >body.txt
#cat header.txt body.txt >PD_TARG_hg19.vcf

gatk --java-options "-Xmx8G" LiftoverVcf \
    -I="PD_TARG_hg19.vcf" \
    -O="$VCFnew" \
    -C="$CHAIN" \
    -R="$REF" \
    --REJECT="/home/stef/CLUSTER/scripts/B_mummy/testrun_mummy_vcf/input_vcf/PD_hg38.reject.vcf" \
    --MAX_RECORDS_IN_RAM=100000

