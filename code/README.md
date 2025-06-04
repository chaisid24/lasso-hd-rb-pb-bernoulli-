# Different code files which were used in the simulation study

The data-gen-code-all.R file generates all Monte-Carlo data sets, for both choices of error distribution and for both homoskedastic and heteroskedastic errors. Detailed description has been provided in the data-gen-code-all.R file to explain the functioning of the code. The following two output files (main data sets used for our Monte-Carlo iterations),

yx-all-n-150-p-500-true-beta-5.Rdata, and 

yx-all-n-300-p-500-true-beta-5.Rdata

obtained by running this code have been stored in the data directory in this github repository.

The rb-pb-oracle-homoskedastic-code-new.R file computes RB based, PB based, Oracle based and Debiased Lasso based CIs for target (true non-zero regression) parameters. As the name suggests, this code is only valid for the homoskedastic errors case. This file contains very detailed comments which explains how the code functions. These comments also help to comprehend the other code file for the heteroskedastic case.

The pb-oracle-heteroskedastic-code-new.R file contains code for PB based, Oracle based and Debiased Lasso based CIs in the heteroskedastic errors case for target (true non-zero regression) parameters. The main difference from the homoskedastic case is the change in construction of PB based pivotal statistics. Otherwise, it is similar to the rb-pb-oracle-homoskedastic-code-new.R code file for the homoskedastic case. The explanations/comments provided in rb-pb-oracle-homoskedastic-code-new.R are helpful in understanding this code. 

At present we do not advise running this code for the (n = 150, p = 500) case, as the higher error levels in the heteroskedastic case, and the relatively small true signal strength can create problems in some Monte Carlo iterations. The code works perfectly for the (n = 300, p = 500) case.

The plot-code-new-homoskedastic-rb-pb-oracle-case.R file is code for generating the plots shown in Figures 1, 2 (in the main article), and Figures 1 and 2 in the supplementary materials. The output (.Rdata file) generated from rb-pb-oracle-homoskedastic-code-new.R needs to be kept in the same directory. 3 and 4

The plot-code-new-heteroskedastic-pb-case.R is code for plotting Figures 3 and 4 in the supplementary materials file. The output from pb-oracle-heteroskedastic-code-new.R needs to be stored in the same directory.

