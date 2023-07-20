# import modules and functions

import numpy as np, pandas as pd, os
from utils.testing import SeqIndTester, compute_hyperparam

# sequential kernelized independence test for our examples

def seq_kern_independence(X: np.ndarray,
                          Y: np.ndarray,
                          lmbd_type: str,
                          n0: int):
    """
    Runs the sequential kernelized independence test and returns the test martingale.
    
    Parameters
    ----------
    X: array_like
        the X observations
    Y: array_like
        the Y observations
    lmbd_type: either 'ONS' or 'aGRAPA'
        method for choosing betting fraction
    n0: even integer
        initial sample size for choosing kernel bandwidth with median heuristic,
        for observations until n0 the bandwidths are 1.
        
    Details
    -------
    Parameters are fixed to their default values in SeqIndTester, apart from the following.
    
    truncation_level:
        0.5 for lmbd_type = 'ONS', 0.9 for lmbd_type = 'aGRAPA'
    kernel_param_x:
        output from median heuristic on first n0 samples
    kernel_param_y:
        output from median heuristic on first n0 samples

    The test martingale for the first n0 observations is set to 1.compute_hyperparam(x).
    
    Returns
    -------
    test_martingale:np.ndarray
        the test martingale
                        
    
    """
    # container for test martingale
    n = X.size
    n_pairs = n//2
    test_martingale = np.ones(n_pairs)

    # specificatinos of the tester
    tester = SeqIndTester()
    tester.lmbd_type = lmbd_type
    if lmbd_type == 'ONS':
        tester.truncation_level = 0.5
    else:
        tester.truncation_level = 0.9

    
    # run the test
    n_pairs = n//2
    tester.kernel_param_x = 1
    tester.kernel_param_y = 1
    for cur_pair in range(1, n0//2, 1):
       tester.process_pair(X[(2*cur_pair):(2*(cur_pair + 1))], Y[(2*cur_pair):(2*(cur_pair + 1))],
                               X[:(2*cur_pair)], Y[:(2*cur_pair)])
       test_martingale[cur_pair] = tester.wealth
       test_martingale[cur_pair + 1] = tester.wealth

    tester.kernel_param_x = compute_hyperparam(X[:n0])
    tester.kernel_param_y = compute_hyperparam(Y[:n0])
    
    for cur_pair in range(n0//2, n_pairs - 1):
        tester.process_pair(X[(2*cur_pair):(2*(cur_pair + 1))], Y[(2*cur_pair):(2*(cur_pair + 1))],
                                X[:(2*cur_pair)], Y[:(2*cur_pair)])
        test_martingale[cur_pair] = tester.wealth
        test_martingale[cur_pair + 1] = tester.wealth
    
    # return test martingale
    return test_martingale
        
    
# import the data from csv file

task_id = os.getenv('SLURM_ARRAY_TASK_ID')
l = np.ceil(int(task_id)/1000)
in_path = 'simulation_examples_' + task_id + '.csv'
df = pd.read_csv(in_path, sep = ',', na_values = [])

# run the test and export the result

n0 = 20
lmbd_type = os.getenv('lmbd_type')

sims = df['sim'].unique()
nsims = sims.size
out = []

for s in sims:
    data = df.loc[df['sim'] == s, ['x', 'y']].to_numpy()
    X = data[:,0]
    Y = data[:,1]
    results = pd.DataFrame(seq_kern_independence(X, Y, lmbd_type, n0),
                           columns = ['martingale'])
    results['l'] = l
    results['sim'] = s
    out.append(results)

out = pd.concat(out)
out.reset_index(inplace = True, drop = True)
out_path = 'simulations_' + lmbd_type + '_' + str(n0) + '_' + task_id + '.csv'
out.to_csv(out_path, sep = ',', index = False)
