#!/usr/bin/env bash

# helper functions for my_MONSTER 

function my_display_help {
    echo "Das war wohl nix :("
}

# --- Report Input variables ----------------------------------------------------------------------
function report_input {

    # --- Reporting parameters ------------------------------------------------------
    echo "##" echo "## Genome-wide Monster wrapper version ${version}"
    printf "%s\t%s\t%s\n" "Info" "Date" "$(date "+%Y.%m.%d_%H:%M:%S")"
    printf "%s\t%s\t%s\n" "Info" "Output folder" "${outDir}"
    printf "%s\t%s\t%s\n" "GeneList" "GeneListFile" "${geneListFile}"
    printf "%s\t%s\t%s\n" "GeneList" "NumberGeneChunks" "$(wc -l $geneListFile | cut -d ' ' -f1)"
    printf "%s\t%s\t%s\n" "GeneList" "CurrentGeneChunks" "${chunkNo}"
    printf "%s\t%s\t%s\n" "GeneList" "CurrentGene" "$(cat ${workingDir}/input_gene.list)"
    printf "%s\t%s\t%s\n" "Input" "GenotypeIndexFile" "${INDEX}"

    printf "%s\t%s\t%s\n" "VariantFiltering" "GencodeFeatures" "${gencode}"
    printf "%s\t%s\t%s\n" "VariantFiltering" "GtexFeatures" "${gtex}"
    printf "%s\t%s\t%s\n" "VariantFiltering" "OverlappingFeatures" "${overlap}"
    printf "%s\t%s\t%s\n" "VariantFiltering" "UpperMAFThreshold" "${maxMAF}"
    printf "%s\t%s\t%s\n" "VariantFiltering" "LowerMAFThreshold" "${minMAF}"
    printf "%s\t%s\t%s\n" "VariantFiltering" "InfoThreshold" "${infoThreshold}"

    printf "%s\t%s\t%s\n" "VariantWeighting" "Score" "${score}"

    printf "%s\t%s\t%s\n" "MonsterInput" "ImputationMethod" "${imputation_method:-BLUP}"
    printf "%s\t%s\t%s\n" "MonsterInput" "KinshipFile" "${kinshipFile}"
    printf "%s\t%s\t%s\n" "MonsterInput" "Phenotype" "${phenotype}"
    printf "%s\t%s\t%s\n" "MonsterInput" "PhenoCovarFile" "${phenotypeFile}"
    printf "%s\t%s\t%s\n" "MonsterInput" "SampleNumber" "$(wc -l ${phenotypeFile} | cut -d ' ' -f1)"
    printf "%s\t%s\t%s\n" "MonsterInput" "MONSTEReig" "${MONSTEReig}"


}

# --- Log more data -------------------------------------------------------------------------------
function log_data {

    # --- Logging results and generated data ----------------------------------------

    printf "%s\t%s\t%s\n" "Input" "ImputeChunks" "$(cat ${workingDir}/impute2needed.txt | cut -f1 | paste -s | sed 's/\t/,/')"
    printf "%s\t%s\t%s\n" "SNPinfo" "uniqueBasesInRoi" "$(cat ${workingDir}/roi.bed | awk '{for(i=$2;i<$3;i++) print $1"\t"i}' | sort | uniq | wc -l )"
    printf "%s\t%s\t%s\n" "SNPinfo" "SNPinROI" "$(wc -l ${workingDir}/temp.oxford_clean.gen | cut -d ' ' -f1)"
    printf "%s\t%s\t%s\n" "SNPinfo" "SNPfailingInfoThreshold" "$(wc -l ${workingDir}/snp_flag_info_bad.txt | cut -d ' ' -f1)"
    printf "%s\t%s\t%s\n" "SNPinfo" "SNPfailingMAFThreshold" "$(wc -l ${workingDir}/snp_flag_maf_bad.txt | cut -d ' ' -f1)"
    
    printf "%s\t%s\t%s\n" "MonsterOutput" "SampleNumber" "$(cut -f 2 ${workingDir}/MONSTER.out | tail -n1)"
    printf "%s\t%s\t%s\n" "MonsterOutput" "VariantNumber" "$(cut -f 3 ${workingDir}/MONSTER.out | tail -n1)"
    printf "%s\t%s\t%s\n" "MonsterOutput" "rho" "$(cut -f 4 ${workingDir}/MONSTER.out | tail -n1)"
    printf "%s\t%s\t%s\n" "MonsterOutput" "PValue" "$(cut -f 5 ${workingDir}/MONSTER.out | tail -n1)"


}

# --- Clean up ------------------------------------------------------------------------------------

function if_exists_delete () {

    file=$1

    if [[ -f $file ]];then
        rm $file
    fi
    return 

}

function clean_up {

    if_exists_delete "${workingDir}/temp.oxford_clean.gen"
    if_exists_delete "${workingDir}/impute2needed.txt"
    if_exists_delete "${workingDir}/input_gene.list"
    if_exists_delete "${workingDir}/MONSTER.geno.txt.gz"
    if_exists_delete "${workingDir}/output.log"
    if_exists_delete "${workingDir}/roi_b37.bed"
    if_exists_delete "${workingDir}/roi_unlift.bed"
    if_exists_delete "${workingDir}/snp_flag_flip_alleles.txt"
    if_exists_delete "${workingDir}/snp_flag_info_bad.txt"
    if_exists_delete "${workingDir}/snp_flag_maf_bad.txt"
    if_exists_delete "${workingDir}/temp_geno.txt.gz"
    if_exists_delete "${workingDir}/temp_header.txt.gz"
    if_exists_delete "${workingDir}/temp.oxford_clean.gen"
    if_exists_delete "${workingDir}/temp.oxford_clean.log"
    if_exists_delete "${workingDir}/temp.oxford_clean.sample"
    if_exists_delete "${workingDir}/tmp_filtered_regions.bed"
    if_exists_delete "${workingDir}/tmp.roi_b37.bed"
    if_exists_delete "${workingDir}/temp_variants.txt"
    if_exists_delete "${workingDir}/temp.weights.txt"
    if_exists_delete "${workingDir}/monster_snps.txt"

}

# --- function to check whether file already in TMPDIR, if not file will be copied there
function isinTMP(){
    
    FILEPATH=$1
    DEST=$2
    FILE=$(basename $FILEPATH)

      if [[ ! -f "$DEST/$FILE" ]]
    then
        cp $FILEPATH $DEST
    fi
}


#TODO write correct display_help
