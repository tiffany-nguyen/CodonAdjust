#!/bin/bash

declare -a iupac=("A" "C" "G" "T" "W" "S" "M" "K" "R" "Y" "B" "D" "H" "V" "N")

TYPE=$1     # type of optimize option. should be "allstop", "nostop", "tag2stop", "tag2gln"

aa_input=$2 # input file for aa profiles
nt_indir=$3 # input directory for IUPAC initial nt frequencies
outdir=$4   # output directory
ths=$5      # rouding option

mkdir -p $outdir

# get number of input AA profiles
col_num=$(head -n 1 ${aa_input} | awk -F, '{print NF}')
LEN=$((col_num - 1))

# run optimization
for code1 in "${iupac[@]}"; do
    for code2 in "${iupac[@]}"; do
        for code3 in "${iupac[@]}"; do
            param="${code1}${code2}${code3}"
            nt_input="${nt_indir}/nt_${param}.csv"
	    param_outdir="${outdir}/${param}"
	    mkdir -p ${param_outdir}
	    Rscript optimize_codon_${TYPE}.R ${aa_input} ${nt_input} ${param_outdir} ${ths}
        done
    done
done

# get best results from all IUPAC codes
python get_best_iupac_codes.py ${outdir} ${LEN}
