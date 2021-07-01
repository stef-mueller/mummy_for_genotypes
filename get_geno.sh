#!/usr/bin/env bash

# source config file which holds path for EIGEN scores
source ./config.txt
# source helper functions
source ./helper.sh

function liftOver_bed_hg38_to_b37(){

    # --- OBJECTIVE: liftover bed file with genomic position in hg38 to hg37

    # --- REQUIREMENTS:
    # liftover

    # --- INPUT:
    # liftOver: liftOver executable 
    # inBED: bed file with hg38 positions
    # CHAIN: file with translation table for genomic positions (from UCSC)
    # outBED: name of bed output file with hg37 positions 
    # unliftBED: name of bed output file with unmapped positions
   
    local liftOver=$1
    local inBED=$2
    local CHAIN=$3
    local outBED=$4
    local unliftBED=$5

    $liftOver/liftOver $inBED $CHAIN tmp.$outBED $unliftBED
    # rm "chr" at start of chromosome identifier
    sed -i 's/chr//g' tmp.$outBED
    # add last column with names (unnecessary for analysis but necessary to fulfilf plink input filter)
    awk '{print $0 , "tmp"}' tmp.$outBED > $outBED
}

function remove_samples(){

    # --- OBJECTIVE: remove samples from Impute2 data using plink2 and extract genomic region of interest
    # Impute2 gen format is without sample header so use plink2 to securely identify
    # samples to remove by name 

    # --- REQUIREMENTS:
    # plink2

    # --- INPUT:
    # CHR: chromosome 
    # GEN: compressed impute2 dosage genotypes
    # SAMPLE: Oxford sample information file, see also "https://www.cog-genomics.org/plink/2.0/formats"
    # REMOVE: text file with sample ids to remove, one sample per line
    # BED: bed file with genomic regions of interest

    local CHR=$1
    local GEN=$2
    local SAMPLE=$3
    local BED=$4
    local REMOVE=$5

    # --- STEP1: read dosage data and sample file (created on the fly in R based on
    # sample order file)
    # need --oxford-single-chr flag since impute2 files not correctly formatted
    # need to create temporary pgen file otherwise error is thrown

    if [[ -z "$REMOVE" ]];
    then
        plink2 --silent \
            --gen $GEN \
            --sample $SAMPLE \
            --oxford-single-chr $CHR \
            --extract range $BED \
            --export oxford \
            --out "temp.oxford_clean"
    else
        plink2 --silent \
            --gen $GEN \
            --sample $SAMPLE \
            --remove $REMOVE \
            --oxford-single-chr $CHR \
            --extract range $BED \
            --export oxford \
            --out "temp.oxford_clean"

    fi
}

function marker_QC(){

    # --- OBJECTIVE: flag SNPs with insufficient imputation certainty, 
    # outside given MAF range or aligned to major allele

    # --- REQUIREMENTS:
    # none

    # --- INPUT:
    # INFO: impute2 info INFO
    # infoThreshold: threshold for impute2 info score metric
    # BED: bed file with regions of interest 
    # maxMAF: maximum MAF threshold, variants with larger MAFs will be removed
    # maxMAF: minimum MAF threshold, variants with smaller MAFs will be removed


    local INFO=$1
    local infoThreshold=$2
    local BED=$3
    local maxMAF=$4
    local minMAF=$5

    if [ -z "$maxMAF" ]
    then
        local maxMAF=1
    fi

    if [ -z "$minMAF" ]
    then
        local minMAF=0
    fi

    # calculate inverse MAF thresholds for variants with wrongly defined major allele
    local OneMinusMaxMAF=$(  echo " 1- $maxMAF" | bc |  awk '{printf "%f", $0}' )
    local OneMinusMinMAF=$(  echo " 1- $minMAF" | bc |  awk '{printf "%f", $0}' )

    # --- STEP1: flag bad info scores 
    if [ -f snp_flag_info_bad.txt ];then rm snp_flag_info_bad.txt ;fi 

    local localCHR START STOP NAME
    while read localCHR START STOP NAME
    do

        awk -v thresh=$infoThreshold ' $7 < thresh  {print $2, $3 }' $INFO |\
        awk -v start=$START -v stop=$STOP '$2 >= start && $2 <= stop {print $1 }' \
        >>snp_flag_info_bad.txt

    done<$BED

    # --- STEP2: flag variants with maf outside given threshold range
    if [ -f snp_flag_maf_bad_one.txt ];then rm snp_flag_maf_bad_one.txt ;fi 
    if [ -f snp_flag_maf_bad_two.txt ];then rm snp_flag_maf_bad_two.txt ;fi 

    local localCHR START STOP NAME
    while read localCHR START STOP NAME
    do

        #first case: correct minor allele definition outside bound
        grep -vf snp_flag_info_bad.txt $INFO |\
        awk -v minMAF=$minMAF -v maxMAF=$maxMAF ' ($6 < 0.5 && $6 > maxMAF) || ($6 < 0.5 && $6 < minMAF) { print $2, $3 }' |\
        awk -v start=$START -v stop=$STOP '$2 >= start && $2 <= stop {print $1 }' \
        >>snp_flag_maf_bad_one.txt

        #second case: wrong minor allele definition outside bound
        grep -vf snp_flag_info_bad.txt $INFO |\
        awk -v OneMinusMaxMAF=$OneMinusMaxMAF -v OneMinusMinMAF=$OneMinusMinMAF ' ($6 > 0.5 && $6 < OneMinusMaxMAF) || ($6 > 0.5 && $6 > OneMinusMinMAF) { print $2, $3 }' |\
        awk -v start=$START -v stop=$STOP '$2 >= start && $2 <= stop {print $1 }' \
        >>snp_flag_maf_bad_two.txt

    done<$BED

    # combine list of variants outside MAF threshold
    cat snp_flag_maf_bad_one.txt snp_flag_maf_bad_two.txt > snp_flag_maf_bad.txt
    rm snp_flag_maf_bad_one.txt
    rm snp_flag_maf_bad_two.txt


    # --- STEP3: flag snps to flip minor + major definition;
    # at same time ensure SNPs fulfil info filter to avoid unnecessary flipping
    if [ -f snp_flag_flip_alleles.txt ];then rm snp_flag_flip_alleles.txt ;fi
    
    local localCHR START STOP NAME
    while read localCHR START STOP NAME
    do
        grep -vf snp_flag_info_bad.txt $INFO |\
        grep -vf snp_flag_maf_bad.txt |\
        awk ' $6 > 0.5 { print $2, $3 }' |\
        awk -v start=$START -v stop=$STOP '$2 >= start && $2 <= stop {print $1 }' \
        >>snp_flag_flip_alleles.txt
    done<$BED

}

function generate_MONSTER_genotypes(){

    # --- OBJECTIVE: generate minor allele dosage file as Input for MONSTER

    # --- REQUIREMENTS:
    # none

    # --- INPUT:

    # CHR: chromosome
    # cleanIMPUTEfile: impute2 dosage genotypes
    # cleanSAMPLEfile: Oxford sample information file, see also "https://www.cog-genomics.org/plink/2.0/formats"
    # SNPinfoBAD: text file with SNP IDs of markers with insufficient info score
    # SNPmafBAD: text file with SNP IDs of markers with MAF outside specified range
    # SNPflip: text file with SNP IDs of markers aligned to major and not minor allele


    local CHR=$1
    local cleanIMPUTEfile=$2
    local cleanSAMPLEfile=$3
    local SNPinfoBAD=$4
    local SNPmafBAD=$5
    local SNPflip=$6

    ### retrieve sample number
    local Nsam=$(tail -n +3 $cleanSAMPLEfile | wc -l | cut -d " " -f1)

    # --- STEP1: make mean geno file for correct minor allele defined markers
    # first remove markers flag for bad info score,
    # then remove markers which should be flipped,
    # then create new output file with first column containing new SNP identifier
    # concatenated from CHR, POS, Major allele and Minor Allele and subsequent
    # columns holding information of imputed minor allele counts, i.e. one column per sample


    cat $cleanIMPUTEfile |\
    grep -vf $SNPinfoBAD |\
    grep -vf $SNPmafBAD |\
    grep -vf $SNPflip |\
    awk -v chr=$CHR -v s=$Nsam '{ printf chr "_" $3 "_" $4 "_" $5 ; for(i=1; i<=s; i++) printf "\t" $(i*3+5)*2+$(i*3+4); printf "\n" }' >geno_one.txt

    # --- STEP2: make mean geno file for wrongly defined minor alleles
    # first select markers flag for flipping of major and minor allele,
    # then create new output file with first column containing new SNP identifier
    # concatenated from CHR, POS, Major allele and Minor Allele (already flipped) 
    # and subsequent columns holding information of imputed minor allele counts,
    #  i.e. one column per sample
    cat $cleanIMPUTEfile |\
    grep -f $SNPflip |\
    awk -v chr=$CHR -v s=$Nsam '{ printf chr "_" $3 "_" $5 "_" $4 ; for(i=1; i<=s; i++) printf "\t" $(i*3+3)*2+$(i*3+4); printf "\n" }' >geno_two.txt

    # --- STEP3: bring geno info together
    cat geno_one.txt geno_two.txt | gzip >temp_geno.txt.gz
    rm geno_one.txt
    rm geno_two.txt

    # --- STEP4: create header file
    echo -n "0	" >temp_header.txt # -n option suppresses new line character
    cut -d " " -f1 $cleanSAMPLEfile | tail -n +3 | paste -d "\t" -s >>temp_header.txt
    gzip temp_header.txt

    # --- STEP5: bring everything together
    cat temp_header.txt.gz temp_geno.txt.gz >MONSTER.geno.txt.gz

    # --- STEP6: clean up
    #rm temp*

}

function add_eigen_scores(){

    # --- OBJECTIVE: Add EIGEN scores to MONSTER snplist file

    # --- REQUIREMENTS:
    # tabix

    # --- ATTENTION: 
    # 1. Score data in hg19!!!
    # 2. Eigen Scores only availbale for SNPs -> need to remove INDELs
    # 3. Eigen and EigenPC scores can contain negative values which will break MONSTER
    
    # --- INPUT:
    # TYPE: which Eigen score format to use, possible options "Eigen", "EigenPC", "EigenPhred", "EigenPCPhred"
    # TABIX: tabix query formatted as chr:pos
    # ALT: alternative allele at specific position
    # REF: reference allele at specific position
    # EigenPath: Path and file name to tabix-ed database containing different Eigen scores

    local TYPE=$1
    local TABIX=$2
    local ALT=$3
    local REF=$4
    local EigenPATH=$5
    

    case $TYPE in
    Eigen )
        local COL=5 ;;
    EigenPC )
        local COL=6 ;;
    EigenPhred )
        local COL=7 ;;
    EigenPCPhred )
        local COL=8 ;;
    *) my_display_help "[Error] Chosen Eigen Score not defined. Please chose from following values: Eigen, EigenPC, EigenPhred, EigenPCPhred" ;;
    esac

    tabix $EigenPATH $TABIX | grep $ALT | grep $REF | cut -f $COL 

} 

function make_MONSTER_snpfile_without_SCORE(){

    # --- OBJECTIVE: create Snplist no scores are added

    # --- REQUIREMENTS:
    # none

    # --- INPUT:
    # MONSTERGeno: compressed MONSTER geno input file

    local MONSTERGeno=$1

    # --- STEP1: create file and save snplist name and set switch for inclusion
    # SNP weights
    echo -n "SNPlist01 0   " >MONSTER.snplist.txt

    # --- STEP2: extract Variant IDs from genotype file 
    zcat $MONSTERGeno |\
        cut -f1 |\
        tail -n +2 >temp_variants.txt

    # --- STEP3: remove INDELs
    ./make_MONSTER_snplist_remove_indel.R -i "temp_variants.txt"

    # --- STEP4: append snp ids to the output file
    cut -f1 "monster_snps.txt" |\
        paste -s >>MONSTER.snplist.txt

}

function make_MONSTER_snpfile_with_SCORE(){

    # --- OBJECTIVE: create Snplist with added pathogenicity scores for each variant

    # --- REQUIREMENTS:
    # none

    # --- INPUT:
    # MONSTERGeno: compressed MONSTER geno input file
    # score: which Eigen score format to use, possible options "Eigen", "EigenPC", "EigenPhred", "EigenPCPhred"

    local MONSTERGeno=$1
    local score=$2

    # --- STEP1: create file and save snplist name and set switch for inclusion
    # SNP weights
    echo -n "SNPlist01 1   " >MONSTER.snplist.txt

     # --- STEP2: extract Variant IDs from genotype file 
    zcat $MONSTERGeno |\
        cut -f1 |\
        tail -n +2 >temp_variants.txt

    # --- STEP3: remove INDELs, creates output "monster_snps.txt"
    ./make_MONSTER_snplist_remove_indel.R -i "temp_variants.txt"
    

    # --- STEP4: append snp ids to the output file
    cut -f1 "monster_snps.txt" |
        paste -s >>MONSTER.snplist.txt

    # --- STEP5: create second line containing weights
    echo -n "SNPlist01 1   " >>MONSTER.snplist.txt
    
    local VAR TABIX ALT REF
    while read VAR TABIX ALT REF
    do
        echo -n "$TABIX $ALT " >>temp.weights.txt
        echo $(add_eigen_scores $score $TABIX $ALT $REF $EigenPath) >>temp.weights.txt

    done<"monster_snps.txt"


    # --- STEP7: append weigths in second line
    cut -d " " -f3 temp.weights.txt |\
        paste -s >>MONSTER.snplist.txt


}
