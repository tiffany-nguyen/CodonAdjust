library(ipoptr)

# beta:
# nucleotide frequencies in the spike-in codon.
# represented as a 12-dimensional vactor.
# these are the variables to be optimized.
#
#     A        C        G        T
#   +------------------------------------
# 1 | beta[1]  beta[2]  beta[3]  beta[4]
# 2 | beta[5]  beta[6]  beta[7]  beta[8]
# 3 | beta[9]  beta[10] beta[11] beta[12]
#
# beta[4*(i-1)+j] 
# represents a frequency at position i of nucleotide j
# where the numbering of nucleotides is the alphabetical order
# i.e.
# A(j=1), C(j=2), G(j=3), T(j=4)

# beta_init:
# initial values of beta at the beginning of optimization.
# use uniform nucleotide frequencies.
#beta_init <- rep(0.25, 12)
# alternatively, use the values calculated by Nguyen-san's procedure
#beta_init <- c(0,0,0.2348,0.7652,0.06,0.1748,0.6126,0.1526,0.210866667,0.394566667,0,0.394566667)

# beta_lb, beta_ub:
# lower and upper bounds of beta
beta_lb <- rep(0.0, 12)
beta_ub <- rep(1.0, 12)

# Tr:
# translation table from codons to AAs.
# represented as a 64-dimensional vector.
# these are used as parameters during optimization.
#
# Tr[i] = j
# represents that codon i is translated into AA j
# where the numbering of codons and AAs is the alphabetical order
# i.e. 
# AAA(i=1), AAC(i=2), AAG(i=3), ..., TTT(i=64)
# *(j=0), A(j=1), C(j=2), D(j=3), ..., Y(j=20)
# amber stop codon TAG(i=51) is translated into Q(j=14)
Tr <- c(9,12,9,12,17,17,17,17,15,16,15,16,8,8,11,8,14,7,14,7,13,13,13,13,15,15,15,15,10,10,10,10,4,3,4,3,1,1,1,1,6,6,6,6,18,18,18,18,0,20,14,20,16,16,16,16,0,2,19,2,10,5,10,5)

# Pa:
# AA frequencies in a machine-learning predction list
# represented as a 20-dimensional vector
# these are used as parameters during optimization.
#
# Pa[i]
# represents a frequency of AA i
# A(i=1), C(i=2), D(i=3), ..., Y(i=20)
#
# for example, use AA frequencies at the pos 1 in wCys_top10k_loop2.csv
#Pa <- c(1748,6126,600,0,0,0,0,0,0,1526,0,0,0,0,0,0,0,0,0,0)
#Pa <- Pa / sum(Pa)

# therefore, the parameter vector has 64+20=84 dimension in total.
#params <- c(Tr, Pa)


# objective function
eval_f <- function(beta, params) {
  Tr <- params[1:64]
  Pa <- params[65:84]
  Pa_calc <- rep(0.0, 20)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
        Pa_calc[Tr[i]] <- Pa_calc[Tr[i]] + beta[b1]*beta[4+b2]*beta[8+b3]
      }
    }
  }
  delta <- Pa_calc - Pa

  return(sum(delta^2))
}


# gradient
eval_grad_f <- function(beta, params) {
  Tr <- params[1:64]
  Pa <- params[65:84]
  Pa_calc <- rep(0.0, 20)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
        Pa_calc[Tr[i]] <- Pa_calc[Tr[i]] + beta[b1]*beta[4+b2]*beta[8+b3]
      }
    }
  }
  delta <- Pa_calc - Pa

  grad <- rep(0.0, length(beta))
  d <- matrix(0.0, nrow=12, ncol=20)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	d[b1, Tr[i]] <- d[b1, Tr[i]] + beta[4+b2]*beta[8+b3]
      }
    }
  }

  for (b2 in 1:4) {
    for (b1 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	d[4+b2, Tr[i]] <- d[4+b2, Tr[i]] + beta[b1]*beta[8+b3]
      }
    }
  }

  for (b3 in 1:4) {
    for (b1 in 1:4) {
      for (b2 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	d[8+b3, Tr[i]] <- d[8+b3, Tr[i]] + beta[b1]*beta[4+b2]
      }
    }
  }

  for (b1 in 1:4) {
    grad[b1] <- sum(2*delta*d[b1,])
  }
  for (b2 in 1:4) {
    grad[4+b2] <- sum(2*delta*d[4+b2,])
  }
  for (b3 in 1:4) {
    grad[8+b3] <- sum(2*delta*d[8+b3,])
  }

  return(grad);
}


# constraint
eval_g <- function(beta, params) {
  return( c(beta[1] + beta[2]  + beta[3]  + beta[4], 
            beta[5] + beta[6]  + beta[7]  + beta[8],
            beta[9] + beta[10] + beta[11] + beta[12],
            beta[4] * (beta[5] + beta[7]) * beta[9]) )
}

# lower and upper bounds of constraints
constraint_lb <- c(1.0, 1.0, 1.0, 0.0)
constraint_ub <- c(1.0, 1.0, 1.0, 0.0)

# Jacobian of constraint
eval_jac_g_structure <- list( c(1,2,3,4), 
                              c(5,6,7,8), 
                              c(9,10,11,12), 
                              c(4,5,7,9) )
eval_jac_g <- function(beta, params) {
  return( c( rep(1.0, 4), 
             rep(1.0, 4),
             rep(1.0, 4),
             (beta[5]+beta[7])*beta[9], beta[4]*beta[9], beta[4]*beta[9], beta[4]*(beta[5]+beta[7]) ) )
}


# Hessian of Lagrangian
eval_h_structure <- list( c(1), 
                          c(1,2), 
                          c(1,2,3), 
                          c(1,2,3,4), 
                          c(1,2,3,4,5), 
                          c(1,2,3,4,5,6), 
                          c(1,2,3,4,5,6,7), 
                          c(1,2,3,4,5,6,7,8), 
                          c(1,2,3,4,5,6,7,8,9), 
                          c(1,2,3,4,5,6,7,8,9,10), 
                          c(1,2,3,4,5,6,7,8,9,10,11), 
                          c(1,2,3,4,5,6,7,8,9,10,11,12) ) 
eval_h <- function(beta, obj_factor, hessian_lambda, params) {
  Tr <- params[1:64]
  Pa <- params[65:84]
  Pa_calc <- rep(0.0, 20)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
        Pa_calc[Tr[i]] <- Pa_calc[Tr[i]] + beta[b1]*beta[4+b2]*beta[8+b3]
      }
    }
  }
  delta <- Pa_calc - Pa

  d <- matrix(0.0, nrow=12, ncol=20)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	d[b1, Tr[i]] <- d[b1, Tr[i]] + beta[4+b2]*beta[8+b3]
      }
    }
  }

  for (b2 in 1:4) {
    for (b1 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	d[4+b2, Tr[i]] <- d[4+b2, Tr[i]] + beta[b1]*beta[8+b3]
      }
    }
  }

  for (b3 in 1:4) {
    for (b1 in 1:4) {
      for (b2 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	d[8+b3, Tr[i]] <- d[8+b3, Tr[i]] + beta[b1]*beta[4+b2]
      }
    }
  }

# second-order gradient regarding two different positions 4 x 4 x 3 = 48
# second-order gradient regarding the same position (1+2+3+4) x 3 = 30
  hess <- rep(0.0, 78)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      dd <- rep(0.0, 20)
      for (aa in 1:20) {
        dd[aa] <- 2*d[b1, aa]*d[4+b2, aa]
      }
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	dd[Tr[i]] <- dd[Tr[i]] + 2*delta[Tr[i]]*beta[8+b3]
      }
      cdx <- b1
      rdx <- 4+b2
      hess[cdx + rdx*(rdx-1)/2] <- obj_factor * sum(dd)
    }
  }

  for (b1 in 1:4) {
    for (b3 in 1:4) {
      dd <- rep(0.0, 20)
      for (aa in 1:20) {
        dd[aa] <- 2*d[b1, aa]*d[8+b3, aa]
      }
      for (b2 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	dd[Tr[i]] <- dd[Tr[i]] + 2*delta[Tr[i]]*beta[4+b2]
      }
      cdx <- b1
      rdx <- 8+b3
      hess[cdx + rdx*(rdx-1)/2] <- obj_factor * sum(dd)
    }
  }

  for (b2 in 1:4) {
    for (b3 in 1:4) {
      dd <- rep(0.0, 20)
      for (aa in 1:20) {
        dd[aa] <- 2*d[4+b2, aa]*d[8+b3, aa]
      }
      for (b1 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
	dd[Tr[i]] <- dd[Tr[i]] + 2*delta[Tr[i]]*beta[b1]
      }
      cdx <- 4+b2
      rdx <- 8+b3
      hess[cdx + rdx*(rdx-1)/2] <- obj_factor * sum(dd)
    }
  }

  for (b1 in 1:4) {
    for (b1_ in b1:4) {
      dd <- rep(0.0, 20)
      for (aa in 1:20) {
        dd[aa] <- 2*d[b1, aa]*d[b1_, aa]
      }
      cdx <- b1
      rdx <- b1_
      hess[cdx + rdx*(rdx-1)/2] <- obj_factor * sum(dd)
    }
  }

  for (b2 in 1:4) {
    for (b2_ in b2:4) {
      dd <- rep(0.0, 20)
      for (aa in 1:20) {
        dd[aa] <- 2*d[4+b2, aa]*d[4+b2_, aa]
      }
      cdx <- 4+b2
      rdx <- 4+b2_
      hess[cdx + rdx*(rdx-1)/2] <- obj_factor * sum(dd)
    }
  }

  for (b3 in 1:4) {
    for (b3_ in b3:4) {
      dd <- rep(0.0, 20)
      for (aa in 1:20) {
        dd[aa] <- 2*d[8+b3, aa]*d[8+b3_, aa]
      }
      cdx <- 8+b3
      rdx <- 8+b3_
      hess[cdx + rdx*(rdx-1)/2] <- obj_factor * sum(dd)
    }
  }

# terms for constraint 4
# beta[4] * (beta[5] + beta[7]) * beta[9]
  cdx <- 4
  rdx <- 5
  hess[cdx + rdx*(rdx-1)/2] <- hess[cdx + rdx*(rdx-1)/2] + hessian_lambda[4] * beta[9]
  cdx <- 4
  rdx <- 7
  hess[cdx + rdx*(rdx-1)/2] <- hess[cdx + rdx*(rdx-1)/2] + hessian_lambda[4] * beta[9]
  cdx <- 4
  rdx <- 9
  hess[cdx + rdx*(rdx-1)/2] <- hess[cdx + rdx*(rdx-1)/2] + hessian_lambda[4] * (beta[5] + beta[7])
  cdx <- 5
  rdx <- 9
  hess[cdx + rdx*(rdx-1)/2] <- hess[cdx + rdx*(rdx-1)/2] + hessian_lambda[4] * beta[4]
  cdx <- 7
  rdx <- 9
  hess[cdx + rdx*(rdx-1)/2] <- hess[cdx + rdx*(rdx-1)/2] + hessian_lambda[4] * beta[4]

  return(hess)
}


########################################################################
## 
## 
########################################################################

read_f <- function(nt_file, aa_file) {
  nt_data <- read.csv(file=nt_file, header=FALSE)
  aa_data <- read.csv(file=aa_file, header=FALSE)
  return(list("nt_data" = nt_data, "aa_data" = aa_data))
}

read_nt <- function(nt_file) {
  nt_data <- read.csv(file=nt_file, header=FALSE)
  return(list("nt_data" = nt_data))
}

read_aa <- function(aa_file) {
  aa_data <- read.csv(file=aa_file, header=FALSE)
  return(list("aa_data" = aa_data))
}

calc_opt_f <- function(beta, Tr) {
  Pa_calc <- rep(0.0, 20)

  for (b1 in 1:4) {
    for (b2 in 1:4) {
      for (b3 in 1:4) {
        i <- 16*(b1-1) + 4*(b2-1) + b3
        if (Tr[i] == 0) next
        Pa_calc[Tr[i]] <- Pa_calc[Tr[i]] + beta[b1]*beta[4+b2]*beta[8+b3]
      }
    }
  }

  m_Pa_calc <- matrix(Pa_calc, nrow=20, ncol=1)
  rownames(m_Pa_calc) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  return(m_Pa_calc)
}

# read input data
args <- commandArgs(trailingOnly = TRUE)
aa_input = args[1]
iupac_dir = args[2]
outdir = args[3]
ths = as.numeric(args[4]) # specify threshold for rounding. This should be a positive integer

nt_dir = paste(iupac_dir, "/nt_", sep="")
iupac = c("A", "C", "G", "T", "W", "S", "M", "K", "R", "Y", "B", "D", "H", "V", "N")

aa_input_data = read_aa(aa_input)
row_num <- nrow(aa_input_data$aa_data) # this should be 20
col_num <- ncol(aa_input_data$aa_data) - 1

best_init_iupac_MSE  = rep(1, col_num)
best_init_iupac_code = rep("", col_num)


###### identify iupac codes which have best init MSE #########
for (code1 in iupac) {
    for (code2 in iupac) {
    	for (code3 in iupac) {
	    param = paste(code1, code2, code3, sep = "")
	    print(param)
	    nt_input = paste(nt_dir, param, ".csv", sep = "")
	    nt_input_data = read_nt(nt_input) 
	    for (i in 1:col_num) {
	        beta_init <- t(nt_input_data$nt_data[, i])
		Pa <- t(aa_input_data$aa_data[, (i+1)])
    		Pa <- Pa / sum(Pa)
    		params <- c(Tr, Pa)

		# updating best init MSE
		tmp_MSE = eval_f(beta_init, params)/row_num
		#print("current MSE: ")
		#print(tmp_MSE)
    		if (tmp_MSE < best_init_iupac_MSE[i] ) {
		   best_init_iupac_MSE[i] = tmp_MSE
		   best_init_iupac_code[i] = param
		}
	    }      
	}
    }
}

###### run CodonAdjust with the above iupac codes ########
# initialize the output
# initialize the output
output_nt_freq = c() # optimal nt frequency
output_nt_init = c() # init nt frequency
output_aa_freq = c() # optimal aa frequency
output_aa_init = c() # init aa frequency
output_MSE_all = c() # optimal MSE
output_MSE_init = c() # initial MSE

output_nt_freq_rounded = c() # rounded optimal nt frequency
output_aa_freq_rounded = c() # rounded optimal aa frequency
output_nt_init_rounded = c() # rounded init nt frequency
output_aa_init_rounded = c() # rounded init aa frequency
output_MSE_all_rounded = c() # rounded optimal MSE

largest_deviations = c() # largest deviations of aa profile
largest_deviations_rounded = c() # largest deviations of rounded aa profile

largest_deviations_init = c() # largest deviations of aa profile with init
largest_deviations_init_rounded = c() # largest deviations of rounded aa profile with init

# create output directory if not exists
dir.create(file.path(outdir), showWarnings = FALSE)

for (i in 1:col_num) {
    nt_input = paste(nt_dir, best_init_iupac_code[i], ".csv", sep = "")
    nt_input_data = read_nt(nt_input) 
    beta_init <- t(nt_input_data$nt_data[, i])
    Pa <- t(aa_input_data$aa_data[, (i+1)])
    Pa <- Pa / sum(Pa)
    params <- c(Tr, Pa)

    tmp_MSE <- eval_f(beta_init, params)
    if (tmp_MSE == 0) next
    init_aa <- calc_opt_f(beta_init, Tr)

    outfile=paste(outdir, "/tag2gln", ".", i, ".optimize.out", sep = "")

    opts <- list("print_level"=5,
                "file_print_level"=5,
#               "derivative_test"="first-order",
                "derivative_test"="second-order",
                "output_file"=outfile)
    res <- ipoptr( x0=beta_init,
                   eval_f=eval_f,
                   eval_grad_f=eval_grad_f,
                   lb=beta_lb,
                   ub=beta_ub,
                   eval_g=eval_g,
                   eval_jac_g=eval_jac_g,
                   constraint_lb=constraint_lb,
                   constraint_ub=constraint_ub,
                   eval_jac_g_structure=eval_jac_g_structure,
                   eval_h=eval_h,
                   eval_h_structure=eval_h_structure,
                   opts=opts,
                   params=params)


    # process optimal results
    beta_opt <- res$solution
    f_opt <- calc_opt_f(beta_opt, Tr)

    aimed_aa <- params[65:84]
    max_dev <- max(abs(f_opt - aimed_aa))
    max_dev_init <- max(abs(init_aa - aimed_aa))

    largest_deviations = c(largest_deviations, max_dev)
    largest_deviations_init = c(largest_deviations_init, max_dev_init)

    output_nt_freq = c(output_nt_freq, beta_opt)
    output_aa_freq = c(output_aa_freq, f_opt)

    output_nt_init = c(output_nt_init, beta_init)
    output_aa_init = c(output_aa_init, init_aa)

    output_MSE_init = c(output_MSE_init, tmp_MSE/row_num)
    output_MSE_all = c(output_MSE_all, eval_f(beta_opt, params)/row_num)

    # rounding process
    beta_opt_rounded <- round(beta_opt, ths)
    f_opt_rounded <- calc_opt_f(beta_opt_rounded, Tr)
    beta_init_rounded <- round(beta_init, ths)
    init_aa_rounded <- calc_opt_f(beta_init_rounded, Tr)

    max_dev_rounded <- max(abs(f_opt_rounded - aimed_aa))
    max_dev_init_rounded <- max(abs(init_aa_rounded - aimed_aa))

    largest_deviations_rounded = c(largest_deviations_rounded, max_dev_rounded)
    largest_deviations_init_rounded = c(largest_deviations_init_rounded, max_dev_init_rounded)

    output_nt_freq_rounded = c(output_nt_freq_rounded, beta_opt_rounded)
    output_aa_freq_rounded = c(output_aa_freq_rounded, f_opt_rounded)

    output_nt_init_rounded = c(output_nt_init_rounded, beta_init_rounded)
    output_aa_init_rounded = c(output_aa_init_rounded, init_aa_rounded)

    output_MSE_all_rounded = c(output_MSE_all_rounded, eval_f(beta_opt_rounded, params)/row_num)
}

##############print out optimal nt frequency#######################
nt_csv=paste(outdir, "/", "nt_opt.all.csv", sep = "")
nt_rounded=paste(outdir, "/", "nt_opt.all_rounded.csv", sep = "")

nt_mat <- matrix(output_nt_freq, nrow = col_num, byrow = TRUE)
nt_mat <- t(nt_mat)
rownames(nt_mat) <- c("A1", "C1", "G1", "T1", "A2", "C2", "G2", "T2", "A3", "C3", "G3", "T3")
write.table(nt_mat, file = nt_csv, sep = ",", col.names = FALSE, quote=FALSE)


nt_mat_rounded <- matrix(output_nt_freq_rounded, nrow = col_num, byrow = TRUE)
nt_mat_rounded <- t(nt_mat_rounded)
rownames(nt_mat_rounded) <- c("A1", "C1", "G1", "T1", "A2", "C2", "G2", "T2", "A3", "C3", "G3", "T3")
write.table(nt_mat_rounded, file = nt_rounded, sep = ",", col.names = FALSE, quote=FALSE)

##############print out optimal nt frequency#######################


##############print out init nt frequency#######################
nt_csv=paste(outdir, "/", "nt_init.all.csv", sep = "")
nt_rounded=paste(outdir, "/", "nt_init.all_rounded.csv", sep = "")

nt_mat <- matrix(output_nt_init, nrow = col_num, byrow = TRUE)
nt_mat <- t(nt_mat)
rownames(nt_mat) <- c("A1", "C1", "G1", "T1", "A2", "C2", "G2", "T2", "A3", "C3", "G3", "T3")
write.table(nt_mat, file = nt_csv, sep = ",", col.names = FALSE, quote=FALSE)


nt_mat_rounded <- matrix(output_nt_init_rounded, nrow = col_num, byrow = TRUE)
nt_mat_rounded <- t(nt_mat_rounded)
rownames(nt_mat_rounded) <- c("A1", "C1", "G1", "T1", "A2", "C2", "G2", "T2", "A3", "C3", "G3", "T3")
write.table(nt_mat_rounded, file = nt_rounded, sep = ",", col.names = FALSE, quote=FALSE)

##############print out init nt frequency#######################


##############print out optimal aa frequency#######################
aa_csv=paste(outdir, "/", "aa_opt.all.csv", sep = "")
aa_rounded=paste(outdir, "/", "aa_opt.all_rounded.csv", sep = "")

aa_mat <- matrix(output_aa_freq, nrow = col_num, byrow = TRUE)
aa_mat <- t(aa_mat)
rownames(aa_mat) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
write.table(aa_mat, file = aa_csv, sep = ",", col.names = FALSE, quote=FALSE)


aa_mat_rounded <- matrix(output_aa_freq_rounded, nrow = col_num, byrow = TRUE)
aa_mat_rounded <- t(aa_mat_rounded)
rownames(aa_mat_rounded) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
write.table(aa_mat_rounded, file = aa_rounded, sep = ",", col.names = FALSE, quote=FALSE)

##############print out optimal aa frequency#######################

##############print out init aa frequency#######################
aa_csv=paste(outdir, "/", "aa_init.all.csv", sep = "")
aa_rounded=paste(outdir, "/", "aa_init.all_rounded.csv", sep = "")

aa_mat <- matrix(output_aa_init, nrow = col_num, byrow = TRUE)
aa_mat <- t(aa_mat)
rownames(aa_mat) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
write.table(aa_mat, file = aa_csv, sep = ",", col.names = FALSE, quote=FALSE)


aa_mat_rounded <- matrix(output_aa_init_rounded, nrow = col_num, byrow = TRUE)
aa_mat_rounded <- t(aa_mat_rounded)
rownames(aa_mat_rounded) <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
write.table(aa_mat_rounded, file = aa_rounded, sep = ",", col.names = FALSE, quote=FALSE)

##############print out init aa frequency#######################


##############print out optimal MSE#######################
MSE_csv=paste(outdir, "/", "MSE_opt.all.csv", sep = "")
MSE_rounded=paste(outdir, "/", "MSE_opt.all_rounded.csv", sep = "")

MSE_mat <- matrix(output_MSE_all, nrow = col_num, byrow = TRUE)
MSE_mat <- t(MSE_mat)
rownames(MSE_mat) <- c("opt_MSE")
write.table(MSE_mat, file=MSE_csv, sep = ",", col.names=FALSE, quote=FALSE)

MSE_mat_rounded <- matrix(output_MSE_all_rounded, nrow = col_num, byrow = TRUE)
MSE_mat_rounded <- t(MSE_mat_rounded)
rownames(MSE_mat_rounded) <- c("opt_MSE_rounded")
write.table(MSE_mat_rounded, file=MSE_rounded, sep = ",", col.names=FALSE, quote=FALSE)

##############print out optimal MSE#######################


##############print out initial MSE#######################
MSE_init_csv=paste(outdir, "/", "MSE_init.all.csv", sep = "")

MSE_mat_init <- matrix(output_MSE_init, nrow = col_num, byrow = TRUE)
MSE_mat_init <- t(MSE_mat_init)

rownames(MSE_mat_init) <- c("init_MSE")
write.table(MSE_mat_init, file=MSE_init_csv, sep = ",", col.names=FALSE, quote=FALSE)
##############print out initial MSE#######################


##############print out largest deviations#######################
largest_deviations_csv=paste(outdir, "/", "largest_deviations.csv", sep = "")
largest_deviations_rounded_csv=paste(outdir, "/", "largest_deviations_rounded.csv", sep = "")

deviations_mat <- matrix(largest_deviations, nrow = col_num, byrow = TRUE)
deviations_mat <- t(deviations_mat)
rownames(deviations_mat) <- c("largest deviations")
write.table(deviations_mat, file=largest_deviations_csv, sep = ",", col.names=FALSE, quote=FALSE)

deviations_mat_rounded <- matrix(largest_deviations_rounded, nrow = col_num, byrow = TRUE)
deviations_mat_rounded <- t(deviations_mat_rounded)
rownames(deviations_mat_rounded) <- c("largest deviations rounded")
write.table(deviations_mat_rounded, file=largest_deviations_rounded_csv, sep = ",", col.names=FALSE, quote=FALSE)
##############print out largest deviations#######################

##############print out init largest deviations#######################
largest_deviations_csv=paste(outdir, "/", "largest_deviations_init.csv", sep = "")
largest_deviations_rounded_csv=paste(outdir, "/", "largest_deviations_init_rounded.csv", sep = "")

deviations_mat <- matrix(largest_deviations_init, nrow = col_num, byrow = TRUE)
deviations_mat <- t(deviations_mat)
rownames(deviations_mat) <- c("largest deviations_init")
write.table(deviations_mat, file=largest_deviations_csv, sep = ",", col.names=FALSE, quote=FALSE)

deviations_mat_rounded <- matrix(largest_deviations_init_rounded, nrow = col_num, byrow = TRUE)
deviations_mat_rounded <- t(deviations_mat_rounded)
rownames(deviations_mat_rounded) <- c("largest deviations_init rounded")
write.table(deviations_mat_rounded, file=largest_deviations_rounded_csv, sep = ",", col.names=FALSE, quote=FALSE)
##############print out init largest deviations#######################
