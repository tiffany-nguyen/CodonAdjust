CodonAdjust
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

	6. 
# Usage
