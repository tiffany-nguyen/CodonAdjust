# CodonAdjust
CodonAdjust is a free software to optimize a nucleotide composition mimicking a certain AA profile. CodonAdjust formulates the optimization of a nucleotide composition as a non-convex optimization problem which minimizes the squared error between the desired and the calculated AA profiles. We provide CodonAdjust with four different programs, which have various customization in practical scenarios such as selecting or avoiding specific codons.

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
## Programs without control factor
	1. For allstop program
		Rscript optimize_codon_allstop.R aa_input nt_input outdir
	2. For nostop program
		Rscript optimize_codon_nostop.R aa_input nt_input outdir
	3. For tag2stop program
		Rscript optimize_codon_tag2stop.R aa_input nt_input outdir
	4. For tag2gln program
		Rscript optimize_codon_tag2gln.R aa_input nt_input outdir

	where,
	- aa_input: input file of targeted amino acids profiles, where each profile is written in a column.
	- nt_input: input file of initial nucleotide frequencies for each aa profile. Each of the initial 
	nucleotide frequency is written in a column.
	- outdir  : output directory to save the optimized results.

	See "sample" folder for an example of aa_input and nt_input.
	
## Programs with control factor
	1. For allstop program
		Rscript optimize_codon_allstop.R aa_input nt_input cf_val outdir
	2. For nostop program
		Rscript optimize_codon_nostop.R aa_input nt_input cf_val outdir
	3. For tag2stop program
		Rscript optimize_codon_tag2stop.R aa_input nt_input cf_val outdir
	4. For tag2gln program
		Rscript optimize_codon_tag2gln.R aa_input nt_input cf_val outdir
		
	where,
	- cf_val  : control factor, which is used to guarantee a certain rate for all targeted amino acids.
	This should be a positive number smaller than 1.0.
	
## Programs for IUPAC code
* Prepare nucleotide input corresponding to all 3375 IUPAC codes.
  python mk_iupac_nt_input.py len
  where len is the number of amino acids in the input AA profile.
  The script uses 3375 IUPAC codes in iupac_code.tar.gz as input.
  The output will be stored in IUPAC_input folder.
* Run optimize program over all 3375 nt input.
  bash run_optimize_codon_iupac.sh aa_input outdir TYPE
  where TYPE specifies the type of optimize programe to use. It should be "allstop", "nostop", "tag2stop", "tag2gln".
	
## Example
* Rscript optimize_codon_allstop.R sample/aa_input.csv sample/nt_input.csv allstop_output
* Rscript optimize_codon_allstop.R sample/aa_input.csv sample/nt_input.csv 0.1 allstop_output
* bash run_optimize_codon_iupac.sh sample/aa_input.csv allstop_output allstop
	
# Reference
* Paper coming soon.

# Support
* If you have any questions regarding to CodonAdjust, please reach me at tiffany.nguyen[at]aist.go.jp
