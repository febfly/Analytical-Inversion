# Analytical Inversion
R script to compute analytical inversion of TROPOMI data to quantify methane emissions from the Permian Basin
Yuzhong Zhang, et al.: Quantifying methane emissions from the largest oil-producing basin in the United States from space, Science Advances, 6, eaaz5120, 2020. https://doi.org/10.1126/sciadv.aaz5120

## script_proc_simu.R
Process GEOS-Chem output and match with TROPOMI observations

## script_emis_based_jac.R
Convert relative jacobian matrix (ppbv CH4/100% perturbation of emission) to absolute jacobian matrix (ppbv CH4/(kg/m2/s) change of emission)

## script_inverse.R
Compute Bayesian inversion to optimize emission spatial distributions. Computation is performed to achieve monthly emissions.
