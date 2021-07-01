#!/usr/bin/env bash

#ADDED: refactored functions in separate scripts for better readability
source ./helper.sh
source ./check_input.sh
source ./get_geno.sh
source ./catch_error.sh
source ./config.txt

# A wrapper script to automate genome-wide burden testing using MONSTER.
# For more information on the applied method see: http://www.stat.uchicago.edu/~mcpeek/software/MONSTER/

# For more information on how to use the wrapper, please see the README.md file or check the repository:
# https://github.com/wtsi-team144/burden_testing

version="V2 [2019_12_20]"
# Get the date and time
today=$(date "+%Y.%m.%d_%H_%M") 

# The variant selector script, that identifies region fo interest: transcript + eQTLS and regulatory regions
regionSelector=my_Burden_testing.pl

# Folder with all the scripts:
export scriptDir="mummy_for_genotypes"

myMONSTER=$(which myMONSTER)

missing_cutoff=1 # Above the set missingness threshold, variants will be excluded. Below the missing genotype will be imputed.
imputation_method='-A' # The default imputation method is BLUP, slowest, but the most accurate. For other options, see documentation.

# Testing if MONSTER is in the path, exiting if not:
if [[ -z "${myMONSTER}" ]]; then echo "[Error] myMONSTER is not in the path. Exiting."; exit; fi


# --- Defining default values ---------------------------------------------------------------------

export infoThreshold=0.7 # By default the info metric imputation quality threshold is 0.9
#export chunkNo=1 # By default we are processing the first chunk.
export maxMAF=1 # By default this is the upper minor allele frequency.
export minMAF=0 # By default this is the lower minor allele frequency.
export MONSTEReig="" # set empty variable to use as switch to check whether file given as input

# --- Capture command line options ----------------------------------------------------------------

# Help message is printed if no command line arguments has been passed:
if [ $# == 0 ]; then my_display_help; fi

# Looping through all command line options:
OPTIND=1
while getopts ":h:L:c:d:p:P:K:I:b:i:g:m:n:s:l:e:x:k:t:o:f:w:j:r:S:E:b:" optname; do
    case "$optname" in
      # Gene list related parameters:
        "L") geneListFile=${OPTARG} ;;
        "c") export chunkNo=${OPTARG} ;;

      # Genotype input options:
        "I" ) INDEX=${OPTARG} ;;
        "r" ) sampleRemove=${OPTARG} ;;
        "S" ) sampleFile=${OPTARG} ;;

      # MONSTER input files:
        "p" ) phenotype=${OPTARG} ;;
        "P" ) phenotypeFile=${OPTARG} ;;
        "K" ) kinshipFile=${OPTARG} ;;
        "E" ) MONSTEReig=${OPTARG} ;;

      # Wrapper option:
        "b" ) keepTemp=${OPTARG};;
        "i" ) infoThreshold=${OPTARG};;
        "o" ) outDir=${OPTARG};;

      # variant filter parameters:
        "g") gencode=${OPTARG} ;;
        "m") maxMAF=${OPTARG} ;;
        "n") minMAF=${OPTARG} ;;
        "s") score=${OPTARG} ;;
        "l") overlap=${OPTARG} ;;
        "e") gtex=${OPTARG} ;;

      # Other parameters:
        "w") rootDir=${OPTARG} ;;
        "h") my_display_help ;;
        "?") my_display_help "[Error] Unknown option $OPTARG" ;;
        ":") my_display_help "[Error] No argument value for option $OPTARG" ;;
        *) my_display_help "[Error] Unknown error while processing options" ;;
    esac
done

# --- checking input files - if any of the tests fails, the script exits.--------------------------
check_input

# Updating working dir, and create folder:
folder="${currentGene}_output_${today}"
workingDir="${rootDir}/${folder}"
mkdir -p "${workingDir}"

# Retrieve selected gene from gene set:
awk -v cn="${chunkNo}" -v cs="1" 'NR > (cn-1)*cs && NR <= cn*cs' ${geneListFile} > ${workingDir}/input_gene.list

# print current gene
printf "%s\n" "#" 
printf "%s\t%s\n" "# Mummy running on:"  "$(cat ${workingDir}/input_gene.list)"
printf "%s\n" "#"

# --- Reporting parameters ------------------------------------------------------
report_input > "${workingDir}/${currentGene}_input_output.log"

# Entering working directory:
cd ${workingDir}

# Generate bed file with region of interest (ROI)
${scriptDir}/${regionSelector}  --build 38 --input input_gene.list --output gene_set_output ${commandOptions} --verbose > output.log

# liftover ROI.bed from nice hg38 to hg19 (genome build of imputed data)
BEDhg38="${workingDir}/roi.bed"
CHAIN="$scriptDir/hg38ToHg19.over.chain"

liftOver_bed_hg38_to_b37 $liftoverPath $BEDhg38 $CHAIN "roi_b37.bed" "roi_unlift.bed"

# identify which impute chunk files to use
BEDb37="roi_b37.bed"

${scriptDir}/identify_impute2_chunks_needed.R -i $INDEX -b $BEDb37

# extract region of interest from impute chunk file and remove samples to be excluded from analysis
CHR=$(head -n1 $BEDb37 | cut -f1)

# --- test how many impute2 chunks necessary to handle --------------------------------------------
NumberImputeChunks=$(wc -l impute2needed.txt | cut -d " " -f1)

# --- remove samples, reformat impute2 files and perform marker QC --------------------------------
### first case: multiple imputed chunks have to be handled
if [[ $NumberImputeChunks -gt 1 ]];
then

    # loop over impute chunks
    for i in $(seq 1 $NumberImputeChunks);
    do

        # remove pre-selected samples (due to insuffiecent quality) from genotype data
        IMPUTEfile=$(sed -n "${i} p" impute2needed.txt | cut -f1)
        IMPUTEfileNAME=$(basename $IMPUTEfile)

        # copy IMPUTEfile to TMPDIR
        isinTMP $IMPUTEfile "$TMPDIR/input"

        remove_samples $CHR "$TMPDIR/input/$IMPUTEfileNAME" $sampleFile $BEDb37 $sampleRemove

        # perform marker QC (remove markers with insuffiecent info metric and identify
        # markers where A1 is not defined as the minor allele

        INFOfile=${IMPUTEfile/.gz/_info}
        marker_QC $INFOfile $infoThreshold $BEDb37 $maxMAF $minMAF

        # rename files to prevent overwriting
        find -maxdepth 1 -regex ".*/temp.oxford_clean.[a-z]*" -exec mv {} {}_$i \;
        find -maxdepth 1 -regex ".*/snp_flag_[._a-z]*" -exec mv {} {}_$i \;
    
    done

    # catenate the genotype files into one file
    find -maxdepth 1 -regex ".*/temp.oxford_clean.gen_[1-9]*" -exec cat {} \; >"temp.oxford_clean.gen"
    # catenate the snp_flag files into one file
    find -maxdepth 1 -regex ".*/snp_flag_info_bad.txt_[1-9]*" -exec cat {} \; >"snp_flag_info_bad.txt"
    find -maxdepth 1 -regex ".*/snp_flag_maf_bad.txt_[1-9]*" -exec cat {} \; >"snp_flag_maf_bad.txt"
    find -maxdepth 1 -regex ".*/snp_flag_flip_alleles.txt_[1-9]*" -exec cat {} \; >"snp_flag_flip_alleles.txt"
    # rename one of the identical sample files
    mv "temp.oxford_clean.sample_1" "temp.oxford_clean.sample"
    #clean up
    rm temp.oxford_clean.gen_*
    rm snp_flag_info_bad.txt_*
    rm snp_flag_flip_alleles.txt_*
    rm snp_flag_maf_bad.txt_*

### second case: only one impute chunk
else

    # remove pre-selected samples (due to insuffiecent quality) from genotype data
    IMPUTEfile=$(cut -f1 impute2needed.txt)
    IMPUTEfileNAME=$(basename $IMPUTEfile)

    # copy IMPUTEfile to TMPDIR
    isinTMP $IMPUTEfile "$TMPDIR/input"

    remove_samples $CHR "$TMPDIR/input/$IMPUTEfileNAME" $sampleFile $BEDb37 $sampleRemove

    # perform marker QC (remove markers with insuffiecent info metric and identify
    # markers where A1 is not defined as the minor allele

    INFOfile=${IMPUTEfile/.gz/_info}
    marker_QC $INFOfile $infoThreshold $BEDb37 $maxMAF $minMAF
fi

### CATCH-ERROR: check there are SNP in ROI
if [[ ! -e temp.oxford_clean.gen ]]
then
      failed "No_imputed_variants_in_ROI"
      exit 1
fi


# --- finally generate the MONSTER genotypes ------------------------------------------------------
cleanIMPUTEfile="temp.oxford_clean.gen"
cleanSAMPLEfile="temp.oxford_clean.sample"
SNPinfoBAD="snp_flag_info_bad.txt"
SNPmafBAD="snp_flag_maf_bad.txt"
SNPflip="snp_flag_flip_alleles.txt"

printf "%s\n" "" "#" "# Now making MONSTER genos!" "#"  ""

generate_MONSTER_genotypes $CHR $cleanIMPUTEfile $cleanSAMPLEfile $SNPinfoBAD $SNPmafBAD $SNPflip

### CATCH-ERROR: check at least 3 variants left after QC
NUMVAR=$(zcat MONSTER.geno.txt.gz | wc -l )
if [[ $NUMVAR -lt 4 ]]
then
      failed "Less_than_3_imputed_variants_left_for_analysis"
      exit 1
fi

### CATCH-ERROR: more than 5000 variants left after QC
NUMVAR=$(zcat MONSTER.geno.txt.gz | wc -l )
if [[ $NUMVAR -gt 5000 ]]
then
      failed "Less_than_5000_imputed_variants_for_analysis"
      exit 1
fi

# --- create MONSTER SNPfile ---------------------------------------------------------------------
MONSTERgeno="MONSTER.geno.txt.gz"

if [[ $score = "None" ]]
then
    # SNPset file without scores
    make_MONSTER_snpfile_without_SCORE $MONSTERgeno

else
    # SNPset file with Eigen scores
    make_MONSTER_snpfile_with_SCORE $MONSTERgeno $score

    ### CATCH-ERROR: check SNP weights have been assigned, there are 
    #                genomic regions in input EIGEN score database which
    #                are misisng
    NUMuniqWEIGTHS=$(cut -d " " -f3 temp.weights.txt | head | sort | uniq | wc -l )
    if [[ $NUMuniqWEIGTHS -lt 2 ]]
    then
          failed "SNP_weights_not_correctly_assigned"
          exit 1
    fi
fi


# --- copy other input files to TMPDIR
isinTMP "$IMONSTERkinship" "$TMPDIR/input"
isinTMP "$MONSTEReig" "$TMPDIR/input"

echo "---------"
echo "[BEFORE MUMMY EXECUTION] TMPDIR before mummy execution"
ls -lhR $TMPDIR
echo "---------"

# --- run MONSTER ---------------------------------------------------------------------------------
MONSTERkinship=$kinshipFile
MONSTERsnp="MONSTER.snplist.txt"
MONSTERpheno=$phenotypeFile

### decompress geno file temporarily
gunzip $MONSTERgeno
MONSTERgenoUNZIP=${MONSTERgeno/.gz/}

# if MONSTER.eig file given use it to speed up computation
if [[ -e $MONSTEReig ]]
then

  myMONSTER \
    -k $MONSTERkinship \
    -p $MONSTERpheno \
    -m 1 \
    -g $MONSTERgenoUNZIP  \
    -s $MONSTERsnp \
    -e $MONSTEReig \
    ${imputation_method}

else 

  myMONSTER \
    -k $MONSTERkinship \
    -p $MONSTERpheno \
    -m 1 \
    -g $MONSTERgenoUNZIP  \
    -s $MONSTERsnp \
    ${imputation_method}

fi

# --- Logging -------------------------------------------------------------------------------------
log_data >> "${workingDir}/${currentGene}_input_output.log"
# cp log file to outdir
cp "${workingDir}/${currentGene}_input_output.log" ${outDir}

# --- Clean up largest temp file ------------------------------------------------------------------
rm ${workingDir}/temp.oxford_clean.gen

# --- Clean up further unnecessary temp files
if [[ "$keepTemp" != "yes" ]]
then
  clean_up
fi

# --- Move data to result directory ----------------------------------------------------------------
# first all files being created while running
mv ${workingDir} ${outDir}
cd "${outDir}"
# create tar archive to reduce disk space and file numbers
tar --remove-files -czf "${folder}.tar.gz" "${folder}"

# create a folder to move tar archives to
if [[ ! -d "TARS" ]]
then
  mkdir "TARS"
fi

# move tar archives to TAR folder
mv "${folder}.tar.gz" "TARS"