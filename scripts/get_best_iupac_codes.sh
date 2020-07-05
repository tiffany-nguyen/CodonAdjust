#!/bin/bash

declare -a iupac=("A" "C" "G" "T" "W" "S" "M" "K" "R" "Y" "B" "D" "H" "V" "N")

indir=$1   # directory that stores optimization results of IUPAC codes
len=$2     # the number of input AA profiles (len = 14 for aa_input in sample)

mkdir -p ${indir}/BEST

for ((i=1; i<=${len}; i++)); do
    for code1 in "${iupac[@]}"; do
	for code2 in "${iupac[@]}"; do
            for code3 in "${iupac[@]}"; do
		param="${code1}${code2}${code3}"
		opt_MSE="${indir}/${param}/MSE_opt.all.csv"
		best_MSE=1
		curr_MSE=``
		mkdir -p ${param_outdir}
		Rscript optimize_codon_${TYPE}.R ${aa_input} ${nt_input} ${param_outdir}
            done
	done
    done
done
