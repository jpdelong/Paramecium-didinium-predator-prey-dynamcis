In this repository, we store the raw data files and data analysis files used in the stability selection paper by DeLong et al.

Files names as "Did_par_##_all.csv" are population density data files for the outcrossed lines numbered ## = 9, 33, 34, 56, 57. The columns indicate the sampling day of the experiment, and data are reported as cells / mL for each replicate and each species. The average densities for paramecia and didinium are also included.

The file named "Mixed_genotype_dishes.csv" contains the population count/density data for the mixed-genotype replicates. The columns indicate the sampling day of the experiment, and data are reported as cells / mL for each replicate and each species. The average densities for paramecia and didinium are also included.

The file named "MorphDataCombined.csv" contains the video phenotyping data of paramecia through time by each replicate. The columns indicate replicate, sampling day, file names, and a range of traits. Lengths and distances are um. Speeds are um per second.

The Julia file "plotting_dynamics_PAR_DID_[date].jl" is the main analysis file for the whole study. This script reads in the data, plots it, and runs some statistical tests. It calls the following functions:
1. "end_points.jl" calculates the ...
2. "plot_gt_dynamics2.jl" is a convenience function for plotting the same thing repeatedly.
3. "plot_morph_dynamics.jl" is a convenience function for plotting the same thing repeatedly.
