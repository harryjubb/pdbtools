#!/usr/bin/python
'''
'''

import sys

from collections import OrderedDict

import numpy as np
from scipy.stats import hmean

def gini(x):
    '''
    `http://www.ellipsix.net/blog/2012/11/the-gini-coefficient-for-distribution-inequality.html`
    '''
    # REQUIRES ALL VALUES IN X TO BE ZERO OR POSITIVE NUMBERS,
    # OTHERWISE RESULTS ARE UNDEFINED
    x = np.array(x)
    n = len(x)
    s = x.sum()
    r = np.argsort(np.argsort(-x)) # calculates zero-based ranks
    return 1 - (2.0 * (r*x).sum() + s)/(n*s)

if __name__ == '__main__':
    
    olg_filename = sys.argv[1]
    summary_tsv_filename = olg_filename + '_sum.tsv'
    per_atom_tsv_filename = olg_filename + '.tsv'
    
    data = OrderedDict()
    
    # PARSE OLG PDB FILE
    with open(olg_filename, 'rb') as fo:
        
        for line in fo:
            
            line = line.strip()
            
            if line.startswith('HETATM'):
                
                chain = line[21:22]
                resnum = line[22:26].strip()
                inscode = line[26:27].strip()
                atom_name = line[11:15].strip()
                
                pymol_string = '/{}/{}{}/{}'.format(chain,
                                                    resnum,
                                                    inscode,
                                                    atom_name)
                
                rinacc = float(line[73:78])
                
                data[pymol_string] = rinacc
    
    # SUMMARY STATISTICS
    all_rinaccs = data.values()
    
    # MIN
    min_rinacc = min(all_rinaccs)
    
    # MAX
    max_rinacc = max(all_rinaccs)
    
    # MEAN
    mean_rinacc = np.mean(all_rinaccs)
    
    # SD
    sd_rinacc = np.std(all_rinaccs)
    
    # MEDIAN
    med_rinacc = np.median(all_rinaccs)
    
    # GINI
    gini_rinacc = gini(all_rinaccs)
    
    # HARMONIC MEAN
    harmonic_mean_rinacc = hmean(all_rinaccs)
    
    # SUMMARY OUTPUT
    output = [min_rinacc, max_rinacc, mean_rinacc, sd_rinacc, med_rinacc,
              gini_rinacc, harmonic_mean_rinacc]
    
    with open(summary_tsv_filename, 'wb') as fo:
        fo.write('{}\n'.format('\t'.join([str(x) for x in output])))
    
    # PER-ATOM OUTPUT
    with open(per_atom_tsv_filename, 'wb') as fo:
        
        for k, v in data.iteritems():
            fo.write('{}\t{}\n'.format(k, v))
    
    
    