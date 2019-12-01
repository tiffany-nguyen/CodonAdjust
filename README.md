# CodonAdjust
CodonAdjust is a free software to optimize a nucleotide composition mimicking a certain AA profile. CodonAdjust formulates the optimization of a nucleotide composition as a non-convex optimization problem which minimizes the squared error between the desired and the calculated AA profiles. We provide CodonAdjust with four different options, which have various customization in practical scenarios such as selecting or avoiding specific codons.

# Required Packages ############################
* R 3.6.x
* IPOPT  >=3.12
* HSL

# Install ######################################
	1. Get IPOPT code from
		https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.13.tgz

	2. Get HSL Archive code from 
		http://hsl.rl.ac.uk/ipopt
		Download coinhsl-archive-2014.01.17.tar.gz
	
	3. Unpack the downloaded files
		tar zxvf Ipopt-3.12.13.tgz
		tar zxvf coinhsl-archive-2014.01.17.tar.gz
		
	4. Run below to install IPOPT
		mv coinhsl-archive-2014.01.17 ./Ipopt-3.12.13/ThirdParty/HSL/coinhsl
		cd ./Ipopt-3.12.13
		mkdir build
		cd build
		../configure 2>&1 | tee configure.log
		make 2>&1 | tee make.log
		make test 2>&1 | tee make.test.log
		make install 2>&1 | tee make.install.log
	
	5. Install IPOPT package to R
		install.packages('~/path/to/Ipopt-3.12.13/build/Ipopt/contrib/RInterface', repos=NULL, type='source')

# Usage
## Options without control factor
	1. For "allstop" option
		Rscript optimize_codon_allstop.R aa_input nt_input outdir
	2. For "nostop" option
		Rscript optimize_codon_nostop.R aa_input nt_input outdir
	3. For "tag2stop" option
		Rscript optimize_codon_tag2stop.R aa_input nt_input outdir
	4. For "tag2gln" option
		Rscript optimize_codon_tag2gln.R aa_input nt_input outdir

	where,
	- aa_input: input file of targeted amino acids profiles, where each profile is written in a column.
	- nt_input: input file of initial nucleotide frequencies for each aa profile. Each of the initial 
	nucleotide frequency is written in a column.
	- outdir  : output directory to save the optimized results.

	See "sample" folder for an example of aa_input and nt_input.
	
## Options with control factor
	1. For "allstop" option
		Rscript optimize_codon_allstop_cf.R aa_input nt_input cf_val outdir
	2. For "nostop" option
		Rscript optimize_codon_nostop_cf.R aa_input nt_input cf_val outdir
	3. For "tag2stop" option
		Rscript optimize_codon_tag2stop_cf.R aa_input nt_input cf_val outdir
	4. For "tag2gln" option
		Rscript optimize_codon_tag2gln_cf.R aa_input nt_input cf_val outdir
		
	where,
	- cf_val  : control factor, which is used to guarantee a certain rate for all targeted amino acids.
	This should be a positive number smaller than 1.0.
	- other parameters (aa_input, nt_input, outdir): described in "Options without control factor"
	
## Option for global optimization search
	1. Move to scripts/iupac_codes folder, and run below command to prepare nucleotide input 
	corresponding to all 3375 IUPAC codes.
		python mk_iupac_nt_input.py len
	where len is the number of amino acids in the input AA profile.
	The script uses 3375 IUPAC codes in iupac_code.tar.gz as input.
	Use tar -xzvf iupac_code.tar.gz to decompress this file before running python script.
	The output will be stored in IUPAC_input folder.
	2. Run optimize program in the scripts folder to find global optimization.
		bash run_optimize_codon_iupac.sh aa_input nt_indir outdir TYPE
	where,
	- nt_indir: specifies path to the IUPAC_input folder generated in step 1. 
	- TYPE:  specifies the type of optimize option to use. 
	It should be "allstop", "nostop", "tag2stop", "tag2gln".
	- other parameters (aa_input, nt_input, outdir): described in "Options without control factor"
		
## Example
* Rscript optimize_codon_allstop.R sample/aa_input.csv sample/nt_input.csv allstop_output
* Rscript optimize_codon_allstop_cf.R sample/aa_input.csv sample/nt_input.csv 0.1 allstop_output
* bash run_optimize_codon_iupac.sh sample/aa_input.csv iupac_codes/IUPAC_input allstop

## Output sample
* Output for optimize_codon_*option*.R, and optimize_codon_*option*_cf.R

Below is an output sample when using optimize_codon_allstop.R
![output_sample](/img/CodonAdjust_output_sample.png)

	- allstop.n.optimize.out:
		Output the optimizing process & its results for input AA profile number n.
	- MSE_init.all.csv:
		MSE between initial AAs (calculated from nt_input) and the input AAs.
	- nt_opt.all.csv:
		Optimized nucleotide frequencies for all AA profiles.
	- nt_opt.all_rounded.csv:
		Values smaller than a threshold of 10^(-15) in nt_opt.all.csv are rounded to 0,
		and output to this file.
	- aa_opt.all.csv:
		Optimized AAs calculated from the optimized nucleotide frequencies.
	- aa_opt.all_rounded.csv:
		Values smaller than a threshold of 10^(-15) in aa_opt.all.csv are rounded to 0,
		and output to this file.
	- MSE_opt.all.csv:
		MSE between optimized AAs and the input AAs.
	- MSE_opt.all_rounded.csv:
		Values smaller than a threshold of 10^(-15) in MSE_opt.all.csv are rounded to 0,
		and output to this file.

* Output for option with IUPAC code

	Optimization results for each IUPAC code will be output to a subfolder in the outdir with  IUPAC code as folder name. The best IUPAC code for each input AA profile is searched from all 3,375 IUPAC codes, be output to the **BEST** folder in the outdir.

	Below is a sample output for IUPAC code with allstop.
	![iupac_sample](/img/CodonAdjust_iupac_sample.png)
		
# Reference
* T.D.N, Y.S, T.K, "CodonAdjust: a software for in silico design of a mutagenesis library with specific amino acid profiles", *submitted*.
