#!/usr/bin/env bash


function failed (){


    # --- OBJECTIVE: general error handling function

    # --- REQUIREMENTS:
    # none

    # --- INPUT:
    # message: error message and reason for failing


    message="${1}"

    echo "[Error] Gene has failed: ${currentGene}"
    echo "[Error] Gene has failed: ${currentGene}" >&2 # Reporting to the error log as well.
    echo "[Error] ${message}"
    echo "[Error] folder moved to failed/FAILED_${folder}"
    echo ""

    # Add error message to log file:
    printf "%s\t%s\n" "ERROR" "${message}" >> "${workingDir}/${currentGene}_input_output.log"

    # Generating failed directory 
    mkdir -p "${outDir}/failed"
    # update folder names to FAILED_${currentGene}_output_${today}
    NewFolder="FAILED_${currentGene}_output_${today}"

    # move folder to FAILED folder
    mv ${workingDir} ${outDir}/failed/${NewFolder}
    # create tar archive to reduce disk space and file numbers
    cd ${outDir}/failed/
    cp ${NewFolder}/${currentGene}_input_output.log ./FAILED_${currentGene}_input_output.log
    tar --remove-files -czf "${NewFolder}.tar.gz" "${NewFolder}"

    return 0
}

