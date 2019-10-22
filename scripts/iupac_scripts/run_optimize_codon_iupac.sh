#!/bin/bash

declare -a iupac=("A" "C" "G" "T" "W" "S" "M" "K" "R" "Y" "B" "D" "H" "V" "N")

aa_input=$1 # input file for aa profiles
outdir=$2   # output directory
TYPE=$3     # type of optimize program. should be "allstop", "nostop", "tag2stop", "tag2gln"

mkdir -p $outdir

for code1 in "${iupac[@]}"; do
    for code2 in "${iupac[@]}"; do
        for code3 in "${iupac[@]}"; do
            param="${code1}${code2}${code3}"
            nt_input="IUPAC_input/nt_${param}.csv"
	    param_outdir="${outdir}/${param}"
	    mkdir -p ${param_outdir}
	    Rscript optimize_codon_${TYPE}.R ${aa_input} ${nt_input} ${param_outdir}
        done
    done
done
