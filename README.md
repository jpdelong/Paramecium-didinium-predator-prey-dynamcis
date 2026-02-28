In this repository, we store the raw data files (in the data subfolder) and data analysis files (in the code subfolder) used in the stability selection paper by DeLong et al.

Files named as "Did_par_##_all.csv" are population density data files for the outcrossed lines numbered ## = 9, 33, 34, 56, 57. The columns indicate the sampling day of the experiment, and data are reported as cells / mL for each replicate and each species. The average densities for paramecia and didinium are also included.

The file named "Mixed_genotype_dishes.csv" contains the population count/density data for the mixed-genotype replicates. The columns indicate the sampling day of the experiment, and data are reported as cells / mL for each replicate and each species. The average densities for paramecia and didinium are also included.

The file named "MorphDataCombined.csv" contains the video phenotyping data of paramecia through time by each replicate. The columns indicate replicate, sampling day, file names, and a range of traits. Lengths and distances are um. Speeds are um per second.

The Julia file "Workflow_dynamics_PAR_DID_02_28_26.jl" is the main analysis file for the whole study. This script reads in the data, plots it, and runs some statistical tests. It is the main (and probably only) script that a user would need to run, unless one wants to change how the plotting looks. It calls the following functions, in this order:
1. "stats_on_reps.jl" calculates the features of the dynamic time series of each replicate trial
2. "plot_morph_dynamics.jl" is a convenience function for plotting the same thing repeatedly.
3. "plot_gt_dynamics.jl" is a convenience function for plotting the same thing repeatedly.
