# CodonAdjust
========
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
	1. For allstop program
		Rscript optimize_codon_allstop.R aa_input nt_input outdir
	2. For nostop program
		Rscript optimize_codon_nostop.R aa_input nt_input outdir
	3. For tag2stop program
		Rscript optimize_codon_tag2stop.R aa_input nt_input outdir
	4. For tag2gln program
		Rscript optimize_codon_tag2gln.R aa_input nt_input outdir
		
	where,
	* aa_input: is an input file of targeted aa profiles, where each profile is written in a column.
	* nt_input: is an input file of initial nt frequencies for each aa profile. Each of the initial nt 
	is written in a column.
	* outdir: output directory to save the optimized results.
	See "sample" folder for a sample of aa_input and nt_input.
	
# Support
* Feel free to reach me at tiffany.nguyen[at]aist.go.jp
