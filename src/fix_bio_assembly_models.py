# IMPORTS
import logging
import operator
import os
import string
import sys

from collections import OrderedDict

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Select
from Bio.PDB.Polypeptide import PPBuilder

class SelectFirstModel(Select):
    '''
    '''
    
    def __init__(self, selected_structure):
    
        self.selected_structure = selected_structure
    
    def accept_model(self, model):
        
        if model is self.selected_structure[0]:
            return 1
        
        return 0

# MAIN
if __name__ == '__main__':
    
    chain_id_pool = list(string.ascii_uppercase + string.digits + string.ascii_lowercase)
    
    # LOAD STRUCTURE
    pdb_path = sys.argv[1]
    pdb_noext, pdb_ext = os.path.splitext(pdb_path)
    pdb_ext = pdb_ext.replace('.', '')
    
    pdb_parser = PDBParser()
    structure = pdb_parser.get_structure(os.path.split(os.path.splitext(pdb_path)[0])[1], pdb_path)
    
    # QUIT IF ONLY ONE MODEL
    #if len(structure) == 1:
    #    print 'Only one model found in {}. Exiting.'.format(pdb_path)
    #    sys.exit()
    
    # GET FIRST MODEL
    first_model = structure[0]
    first_model_chain_ids = [chain.id for chain in first_model.child_list]
    
    # REMOVE CHAIN IDS FROM THE FIRST MODEL FROM THE POTENTIAL POOL
    for chain_id in first_model_chain_ids:
        chain_id_pool.remove(chain_id)
    
    # FOR EVERY OTHER MODEL
    for model in structure:
        
        if model is first_model:
            continue
    
        for chain in model.child_list:
            
            # RENAME THE CHAIN
            chain.id = chain_id_pool.pop(0)
            
            # PUT THE CHAIN INTO THE FIRST MODEL
            first_model.add(chain)
    
    # WRITE OUT THE FIRST MODEL
    io = PDBIO()
    io.set_structure(structure)
    
    io.save('.'.join((pdb_noext, 'fixed_models', pdb_ext)), select=SelectFirstModel(structure))
    
    
    