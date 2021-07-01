#!/usr/bin/env bash

# function to check input files

function check_input () {

    # Kinship file (Is it set? Does it exist?):
    if [[ -z "${kinshipFile}" ]]; then
        my_display_help "[Error] Kinship file has to be specified!";
    elif [[ ! -e "${kinshipFile}" ]]; then
        echo "[Error] Kiship file could not be opened: $kinshipFile";
        exit;
    fi                                                                                                                                      

    # Gene list file (Is it set? Does it exist?):
    if [[ ! -e "${geneListFile}" ]]; then
        echo "[Error] Gene list file could not be opened: $geneListFile";
        exit;
    fi

    # Working directory (Is it set? Does it exist?):
    if [[ -z "${rootDir}" ]]; then
        rootDir=$(pwd)
    elif [[ ! -d "${rootDir}" ]]; then
        echo "[Error] The directory does not exists: $rootDir";
        exit;
    fi

    # GENCODE -expecting a list of feature names separated by a comma.
    if [[ ! -z "${gencode}" ]]; then
        commandOptions="${commandOptions} --GENCODE ${gencode}"
    fi

    # GTEx - expecting a list of feature names separeted by comma.
    if [[ ! -z "$gtex" ]]; then
        commandOptions="${commandOptions} --GTEx ${gtex}"
    fi

    # Overlap - expecting a list of features separeated by comma.
    if [[ ! -z "${overlap}" ]]; then
        commandOptions="${commandOptions} --overlap ${overlap}"
    fi

    # MAF - expecting a float between 0 and 0.5
    if [[ ! -z "$MAF" ]]; then
        commandOptions="${commandOptions} --maf ${MAF}"
    fi

    # IF loftee is set, only those variants are included in the test that were predicted to be LC or HC lof variants by loftee:
    if [[ ! -z "$loftee" ]]; then
        commandOptions="${commandOptions} --loftee "
    fi

    # IF lofteeHC is set, only those variants are included in the test that were predicted to be HC lof variants by loftee:
    if [[ ! -z "$lofteeHC" ]]; then
        commandOptions="${commandOptions} --lofteeHC "
    fi

    # If lof is set, only variants with severe consequences will be selected.
    if [[ ! -z "$lof" ]]; then
        commandOptions="${commandOptions} --lof "
    fi

    # Score - If score is not given we apply no score. Otherwise we test the submitted value:
    # Accepted scores:
    if [[ ! -z "${score}" ]]; then
        score="${score^^}"
        case "${score}" in
            EIGEN )        score="Eigen";;
            EIGENPC )      score="EigenPC";;
            EIGENPHRED )   score="EigenPhred";;
            EIGENPCPHRED ) score="EigenPCPhred";;
            CADD )         score="CADD";;
            None )         score="None";;
            * )            score="None";;
        esac
    else
        echo "[Warning] Submitted score is not recognized! Accepted scores: CADD, Eigen, EigenPC, EigenPhred, EigenPCPhred or None."
        echo "[Warning] No scores are being applied."
        score="None"
    fi

    # Only adding score to command line if score is requested:
    if [[ "${score}" != "None" ]]; then
        commandOptions="${commandOptions} --score ${score}";
    fi

    # If Eigen score is applied, we shift the scores by 1, if no other value is specified:
    if [[ ("${score}" == "Eigen") && (-z "${scoreshift}" ) ]]; then scoreshift=1; fi

    # Exons might be extended with a given number of residues:
    if [[ ! -z "${xtend}" ]]; then
        commandOptions="${commandOptions} --extend ${xtend}"
    fi

    # Setting score cutoff, below which the variant will be removed from the test:
    if [[ ! -z "${cutoff}" ]]; then
        commandOptions="${commandOptions} --cutoff ${cutoff}"
    fi

    # Setting score shift, a number that will be added to every scores (MONSTER does not accept scores <= 0!!):
    if [[ ! -z "${scoreshift}" ]]; then
        commandOptions="${commandOptions} --shift ${scoreshift}"
    fi

    # Checking if phenotype is provided and if that phenotype file:
    if [[ -z "${phenotype}" ]]; then
        echo "[Error] Phenotype was not set! Exiting."
        exit 1
    fi
}