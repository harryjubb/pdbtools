#!/usr/bin/python
'''
'''

import sys

from collections import OrderedDict

if __name__ == '__main__':
    
    atom_sasa_filename = sys.argv[1]
    per_atom_tsv_filename = atom_sasa_filename + '.tsv'
    
    data = OrderedDict()
    
    # PARSE ATOM SASA PDB FILE
    with open(atom_sasa_filename, 'rb') as fo:
        
        for line in fo:
            
            line = line.strip()
            
            if line.startswith('HETATM') or line.startswith('ATOM'):
                
                chain = line[21:22]
                resnum = line[22:26].strip()
                inscode = line[26:27].strip()
                atom_name = line[11:15].strip()
                
                pymol_string = '/{}/{}{}/{}'.format(chain,
                                                    resnum,
                                                    inscode,
                                                    atom_name)
                
                solv_acc = float(line[60:66])
                
                data[pymol_string] = solv_acc
    
    # PER-ATOM OUTPUT
    with open(per_atom_tsv_filename, 'wb') as fo:
        
        for k, v in data.iteritems():
            fo.write('{}\t{}\n'.format(k, v))
    
    
    