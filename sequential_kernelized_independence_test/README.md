# Code for the sequential kernelized independence test (SKIT)

This folder contains code for the simulations with the
[sequential kernelized independence test](https://arxiv.org/abs/2212.07383) (SKIT).

The code is copied from the repository
https://github.com/a-podkopaev/Sequential-Kernelized-Independence-Testing/tree/main.
The function  `simulations_skit.py` is new and contains a function to run the
SKIT on the simulation examples from the paper.

The code is structured in such a way that the data is first generated in R
with the files `simulation_examples.R` and `simulate_data_skit.R` to
have the same pseudo-random numbers as for the methods implemented in R. Then
the SKIT is applied with python code, and the results are exported.

To run the simulations, upload the contents of this directory to a HPC cluster
and run `sbatch array_task.sh`, possibly adjusting the `#SBATCH` options to the
specific cluster settings. To run single simulations locally, use
`bash run_locally.sh`.
  
