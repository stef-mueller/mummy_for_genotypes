#!/usr/bin/env bash

#
# Objective: Submission script for mummy pipeline based on genotype data
#

# source helper functions
source ./helper.sh

# --- derive variables from submit script
export currentGene=$1
SUBSET=$2

InputDir="<path/to/input/directory>"
GENE_LIST="$InputDir/ordered_gene_list.txt"
CHUNK_COUNT=$(wc -l $GENE_LIST | cut -d " " -f1)

# --- copy files to TMPDIR
isinTMP "$InputDir/$SUBSET.MONSTER.pheno_and_covar.sorted.txt" "$TMPDIR/input"
isinTMP "$InputDir/$SUBSET.samples_to_remove.remove" "$TMPDIR/input"
isinTMP "$InputDir/onco.sample" "$TMPDIR/input"

# --- run mummy
mummy_for_genotypes/my_MONSTERgenome-wideV2.sh \
    -g "gene,UTR,CDS,exon,transcript" \
    -e "promoter,enhancer,TF_bind" \
    -l "promoter,enhancer,TF_bind" \
    -c $currentGene \
    -I "$InputDir/index_impute2_genos.txt" \
    -K "$InputDir/MONSTER.kin.subset$SUBSET.txt" \
    -P "$TMPDIR/input/$SUBSET.MONSTER.pheno_and_covar.sorted.txt" \
    -r "$TMPDIR/input/$SUBSET.samples_to_remove.remove" \
    -S "$TMPDIR/input/onco.sample" \
    -p "BreastCancer" \
    -s "EigenPCPhred" \
    -L ${GENE_LIST} \
    -E "$InputDir/MONSTER.$SUBSET.eig" \
    -i "0.7" \
    -m "0.05" \
    -n "0" \
    -o "$TMPDIR/output" \
    -b "no" \
    -w "$TMPDIR"

# --- cp output back to project space
cp $TMPDIR/output/*log  "</path/to/result/directory>"
cp $TMPDIR/output/TARS/*gz  "</path/to/result/directory>"

# --- check wether output in failed
if ls $TMPDIR/output/failed/FAILED* 1> /dev/null 2>&1;
then
    cp $TMPDIR/output/failed/FAILED*  "</path/to/result/directory"
fi

echo "---------"
echo "[AFTER MUMMY] TMPDIR after mummy"
ls -lhR $TMPDIR

# --- clean up outputdir
rm -r $TMPDIR/output/*

echo "---------"
echo "[AFTER CLEAN] TMPDIR after cleanup"
ls -lhR $TMPDIR

# --- input explanation:
##########################################


# -c: [integer] currentGeneChunk: ie. on which line of the geneListFile to operate on
# -L: [path] geneListFile: text file with gene symbols, one gene per line
# -I: [path] INDEX: path to genotype index file
# -r: [path] sampleRemove: path to file with samples FID and IID to remove due to QC thresholds
# -S: [path] sampleFile: path to sample file
# -p: [string] phenotype: studied phenotype
# -P: [path] phenotypeFile: path to MONSTER phenotype and covariate file
# -K: [path] kinshipFile:path to MONSTER kinship file
# -E: [path] MONSTEReig: path to MONSTER.eig file, a binary file that stores the eigenvalues and eigenvectors of the kinship
#            matrix, speeds up computational time tremendously
# -b: [string] keepTemp: keep temporary files? chose "yes" to keep temporary files
# -i: [numeric] infoThreshold: INFO score threshold for marker QC of imputed variant, number in range of 0-1
# -o: [path] outDir: path to output directory
# -g: [string] gencode: chosen Gencode feature, one or more of "gene,UTR,CDS,exon,transcript" 
# -m: [numeric] maxMAF: upper minor allele frequency threshold for imputed markers, number in range of 0-1
# -n: [numeric] minMAF: lower minor allele frequency threshold for imputed markers, number in range of 0-1
# -s: [string] score: which EIGEN scores to use? chose one of "None, Eigen, EigenPC, EigenPhred, EigenPCPhred"
# -l: [string] overlap: chosen regulatory features associated by overlap, one or more of "promoter,enhancer,TF_bind" 
# -e: [string] gtex: chosen regulatory features associated by GTEx feature, one or more of "promoter,enhancer,TF_bind"
# -w: [string] rootDir: chosen regulatory features associated by GTEx feature, one or more of "promoter,enhancer,TF_bind"