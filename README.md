# A Rank-Based Sequential Test of Independence

This repository contains replication material for the preprint


The files `sequential_tests.cpp` and `r_wrappers.R` contain the functions to
run the tests. They additionally contain another method for achieving uniform
marginals based on constrained maximum likelihood, which does not appear in the
article. For the BET, we used `BET_functions_JASA0223.r`, which is from the
supplementary materials of *Zhang, K. (2019). BET on independence. Journal of
the American Statistical Association, 114(528), 1620-1637.* Data for the
simulation example is generated with `simulation_examples.R`.

The simulations were run on a HPC cluster with SLURM; `simulations_cluster.R`
produces a single simulation, and `collect_simulation_results.R` collects
the results. Rejection rates for the BET are obtained from `simulations_BET.R`.
The file `ville_gap.R` estimates the corrected rejection thresholds
when stopping at a given sample size N, and produces the corresponding figure
in the paper. 

Plots and tables are generated with the following files: `plot_simulations.R`,
`illustration_theoretical_properties.R` and `illustration_sequential_BET`
generate the figures that illustrate the methods; `plot_simulation_results.R`
gives the rejection rates for the sequential tests; and `compare_against_bet.R`
the table for the comparison with the BET.

