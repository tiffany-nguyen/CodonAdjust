#!/bin/bash

declare -a iupac=("A" "C" "G" "T" "W" "S" "M" "K" "R" "Y" "B" "D" "H" "V" "N")

aa_input=$1 # input file for aa profiles
nt_indir=$2 # input directory for IUPAC initial nt frequencies
cf_val=$3   # control factor (0 < cf_val < 1)
outdir=$4   # output directory
TYPE=$5     # type of optimize option. should be "allstop", "nostop", "tag2stop", "tag2gln"

mkdir -p $outdir

# get number of input AA profiles
col_num=$(head -n 1 ${aa_input} | awk -F, '{print NF}')
LEN=$((col_num - 1))

for code1 in "${iupac[@]}"; do
    for code2 in "${iupac[@]}"; do
        for code3 in "${iupac[@]}"; do
            param="${code1}${code2}${code3}"
            nt_input="${nt_indir}/nt_${param}.csv"
	    param_outdir="${outdir}/${param}"
	    mkdir -p ${param_outdir}
	    Rscript optimize_codon_${TYPE}_cf.R ${aa_input} ${nt_input} ${cf_val} ${param_outdir}
        done
    done
done

# get best results from all IUPAC codes
python get_best_iupac_codes.py ${outdir} ${LEN}
