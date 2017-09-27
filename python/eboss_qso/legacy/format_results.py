import argparse
import pandas as pd
import numpy as np
from nbodykit.dataset import DataSet
from nbodykit.lab import *
import os

def main(ns):
    
    # results
    _, ext = os.path.splitext(ns.result_file)
    if ext == '.dat':
        r = DataSet.from_plaintext(['k'], ns.result_file)
        attrs = r.attrs
    elif ext == '.json':
        res = ConvolvedFFTPower.load(ns.result_file)
        r = res.poles
        attrs = res.attrs
    else:
        raise ValueError("extension should be .dat or .json")
    
    # compute the shot noise
    names = ['ra', 'dec', 'z', 'nbar', 'weight_fkp', 'weight', 'weight_sys']
    data = pd.read_csv(ns.data_file, delim_whitespace=True, engine='c', header=None, names=names)
    
    # if all close pairs are correct (S2 > S1)
    S1 = (data['weight_fkp']**2*data['weight']*data['weight_sys']).sum()
    
    # if all close pairs are wrong (S2 > S1)
    S2 = (data['weight_fkp']**2*data['weight']**2).sum()
        
    # compute the shot noises
    if ext == '.dat': 
        assert np.allclose(attrs['S_data'], S2)
        Pshot2 = (attrs['S_ran'] + S1) / attrs['A_ran']
        Pshot1 = (attrs['S_ran'] + S2) / attrs['A_ran']
    else:
        S1 /= attrs['randoms.A']
        S2 /= attrs['randoms.A']
        assert np.allclose(attrs['data.S'], S2)
        Pshot2 = attrs['randoms.S'] + S1
        Pshot1 = attrs['randoms.S'] + S2
    
    # output data
    kcen = r.coords['k']
    idx = kcen <= ns.kmax
    kmean = r['k'][idx]
    P0 = r['power_0'][idx].real
    P2 = r['power_2'][idx].real
    N = r['modes'][idx]
    out = np.vstack([kcen[idx], kmean, P0-Pshot1, P0-Pshot2, P2, N]).T
    
    # write
    with open(ns.output, 'wb') as ff:
        
        ff.write(b"# k_cen k_mean P0-Pshot1 P0-Pshot2 P2 N_modes\n")
        ff.write(b"# Pshot1 = %.6f (assuming 100%% CP accuracy)\n" %Pshot1)
        ff.write(b"# Pshot2 = %.6f (assuming 0%% CP accuracy)\n" %Pshot2)
        ff.write(b"# alpha = %.6f\n" %attrs['alpha'])
        if ext == '.dat':
            ff.write(b"# normalization = %.6f\n" %attrs['A_ran'])
            ff.write(b"# N_data (unweighted) = %d\n" %attrs['N_data'])
            ff.write(b"# N_data (weighted) = %d\n" %attrs['W_data'])
            ff.write(b"# N_randoms (unweighted) = %d\n" %attrs['N_ran'])
            ff.write(b"# N_randoms (weighted) = %d\n" %attrs['W_ran'])
        else:
            ff.write(b"# normalization = %.6f\n" %attrs['randoms.A'])
            ff.write(b"# N_data (unweighted) = %d\n" %attrs['data.N'])
            ff.write(b"# N_data (weighted) = %d\n" %attrs['data.W'])
            ff.write(b"# N_randoms (unweighted) = %d\n" %attrs['randoms.N'])
            ff.write(b"# N_randoms (weighted) = %d\n" %attrs['randoms.W'])
        
        np.savetxt(ff, out, fmt=['%.4f', '%.7e', '%.7e', '%.7e', '%.7e', '%d'])
    
    

if __name__ == '__main__':
    
    desc = "format the results into a format suggested by Hector"
    parser = argparse.ArgumentParser(description=desc)
    
    h = 'the result file to load'
    parser.add_argument('result_file', type=str, help=h)
    
    h = 'the data file to load'
    parser.add_argument('data_file', type=str, help=h)
    
    h = 'the output file'
    parser.add_argument('-o', '--output', type=str, help=h)
    
    h = 'the kmax value to use'
    parser.add_argument('--kmax', type=float, default=0.4, help=h)
    
    ns = parser.parse_args()
    main(ns)
    
    