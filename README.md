# stPCA


***********************************************************************************************************

1. System requirements

a. The codes can be run within MATLAB environment and R/Rstudio on any operating system.

b. We implemented the codes with MATLAB R2022b, R on Mac/Unbuntu.

c. No non-standard hardware is required.

***********************************************************************************************************

2. Installation guide

a. The MATLAB codes can be run directly without installation.

b. The required R packages are listed in the first few lines of code.

***********************************************************************************************************

3. Sample codes

I. For coupled-Lorenz model simulation, run "main_Lorenz_stPCA_PCA";

II. For multiple-nodes model simulation, run "main_MultipleNodes_stPCA".

III. For the MIMIC application, run "process_multiple_admissions.R" to preprocesses the original table, and run "mimic_DR_stPCA.m" to get 1d representation.

IV. For the singlecell application, run "pseudotime_ordering_singlecells.m" to order the single cells, and run "singlecell_DR_stPCA.m" to get 1d representation.

V. For the T2D application, run "colas_DR_stPCA.m" to get 1d representation.
