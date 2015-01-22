# IMPORTS
import logging
import operator
import os
import sys

from collections import OrderedDict

from Bio.PDB import PDBParser
from Bio.PDB import NeighborSearch
from Bio.PDB.Polypeptide import PPBuilder

from numpy.linalg import norm

# CONSTANTS
COVALENT_RADII = {
	'H': 0.31, 'HE': 0.28, 'LI': 1.28, 'BE': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'NE': 0.58, 'NA': 1.66,
	'MG': 1.41, 'AL': 1.21, 'SI': 1.11, 'P': 1.07, 'S': 1.05, 'CL': 1.02, 'AR': 1.06, 'K': 2.03, 'CA': 1.76, 'SC': 1.7, 'TI': 1.6,
	'V': 1.53, 'CR': 1.39, 'MN': 1.39, 'FE': 1.32, 'CO': 1.26, 'NI': 1.24, 'CU': 1.32, 'ZN': 1.22, 'GA': 1.22, 'GE': 1.2, 'AS': 1.19,
	'SE': 1.2, 'BR': 1.2, 'KR': 1.16, 'RB': 2.2, 'SR': 1.95, 'Y': 1.9, 'ZR': 1.75, 'NB': 1.64, 'MO': 1.54, 'TC': 1.47, 'RU': 1.46,
	'RH': 1.42, 'PD': 1.39, 'AG': 1.45, 'CD': 1.44, 'IN': 1.42, 'SN': 1.39, 'SB': 1.39, 'TE': 1.38, 'I': 1.39, 'XE': 1.4, 'CS': 2.44,
	'BA': 2.15, 'LA': 2.07, 'CE': 2.04, 'PR': 2.03, 'ND': 2.01, 'PM': 1.99, 'SM': 1.98, 'EU': 1.98, 'GD': 1.96, 'TB': 1.94, 'DY': 1.92,
	'HO': 1.92, 'ER': 1.89, 'TM': 1.9, 'YB': 1.87, 'LU': 1.87, 'HF': 1.75, 'TA': 1.7, 'W': 1.62, 'RE': 1.51, 'OS': 1.44, 'IR': 1.41,
	'PT': 1.36, 'AU': 1.36, 'HG': 1.32, 'TL': 1.45, 'PB': 1.46, 'BI': 1.48, 'PO': 1.4, 'AT': 1.5, 'RN': 1.5, 'FR': 2.6, 'RA': 2.21,
	'AC': 2.15, 'TH': 2.06, 'PA': 2.0, 'U': 1.96, 'NP': 1.9, 'PU': 1.87, 'AM': 1.8, 'CM': 1.69, 'BK': 1.6, 'CF': 1.6, 'ES': 1.6,
	'FM': 1.6, 'MD': 1.6, 'NO': 1.6, 'LR': 1.6, 'RF': 1.6, 'DB': 1.6, 'SG': 1.6, 'BH': 1.6, 'HS': 1.6, 'MT': 1.6, 'DS': 1.6, 'RG': 1.6,
	'CN': 1.6, 'UUT': 1.6, 'FL': 1.6, 'UUP': 1.6, 'LV': 1.6, 'UUH': 1.6, 'UUH': 1.6
}

MAX_COV_RADIUS = max(COVALENT_RADII.values())
MAX_COV_BOND = MAX_COV_RADIUS * 2

PDB_LINE_TEMPLATE = '{record: <6}{serial: >5} {atom_name: ^4}{altloc: ^1}{resname: ^3} {chain_id: ^1}{resnum: >4}{icode: ^1}   {x: >8.3f}{y: >8.3f}{z: >8.3f}{occ: >6.2f}{tfac: >6.2f}          {element: >2}{charge: >2}'

# FILENAME MUNGING
pdb_path = sys.argv[1]
pdb_noext, pdb_ext = os.path.splitext(pdb_path)
pdb_ext = pdb_ext.replace('.', '')

# LOAD THE PDB
pdb_parser = PDBParser()
structure = pdb_parser.get_structure(os.path.split(os.path.splitext(pdb_path)[0])[1], pdb_path)

# EXTRACT THE MODEL
model = structure[0]

# ASSIGN COVALENT RADII
for atom in model.get_atoms():
    
    try:
        atom.covrad = COVALENT_RADII[atom.element.strip().upper()]
    except:
        print 'Covalent radius could not be determined for atom {}'.format(atom)

# DETERMINE POLYPEPTIDES AND CHAIN BREAKS
ppb = PPBuilder()
polypeptides = ppb.build_peptides(model, aa_only=False)

# GET ALL POLYPEPTIDE RESIDUES IN THE MODEL
polypeptide_residues = set([])

for polypeptide in polypeptides:
    for residue in polypeptide:
        polypeptide_residues.add(residue)

# GET NON-POLYPEPTIDE HETEROATOM RESIDUES
heteroresidues = []

for residue in model.get_residues():
    
    if 'H_' in residue.get_full_id()[3][0] and residue not in polypeptide_residues:
        heteroresidues.append(residue)

# DISCARD COVALENTLY BOUND RESIDUES
ns = NeighborSearch(list(model.get_atoms()))

non_cov_heteroresidues = []

for residue in heteroresidues:
    
    #print '\n', residue
    residue_is_cov = False
    
    for atom in residue.child_list:
        
        nearby_atoms = ns.search(atom.coord, MAX_COV_BOND)
                
        if nearby_atoms:
            for nearby_atom in nearby_atoms:
                
                # NO SELFIES!
                # BECAUSE OF COURSE THE ATOM WILL COV-BONDED IN IT'S OWN RESIDUE...
                if nearby_atom in residue.child_list:
                    continue
                
                sum_cov_radii = atom.covrad + nearby_atom.covrad
                distance = norm(atom.coord - nearby_atom.coord)
                
                if distance <= sum_cov_radii:
                    
                    #print atom, nearby_atom, sum_cov_radii, distance
                    residue_is_cov = True
                    break
        
        if residue_is_cov:
            break
    
    if not residue_is_cov:
        non_cov_heteroresidues.append(residue)

# LIMIT TO > 7 HEAVY ATOM LIGANDS
filtered_heteroresidues = []

for residue in non_cov_heteroresidues:
    
    if len(residue.child_list) > 7:
        filtered_heteroresidues.append(residue)

# WRITE OUT THE LIGANDS AS PDB FILES
with open(pdb_noext + '.liglist', 'wb') as llfo:
	for residue in filtered_heteroresidues:
	    
	    chain_id = residue.get_parent().id
	    resname = residue.resname
	    resnum = residue.get_full_id()[3][1]
	    inscode = residue.get_full_id()[3][2]
	    
	    if inscode == ' ':
		inscode = 'i'
	    
	    ligand_filename = '{}_{}_{}_{}_{}.{}'.format(pdb_noext,
							 chain_id,
							 resname,
							 resnum,
							 inscode,
							 pdb_ext)
	    
	    # WRITE OUT THE EXTRACTED LIGAND
	    llfo.write(ligand_filename + '\n')
	    with open(ligand_filename, 'wb') as fo:
		
		for atom in residue.child_list:

                    # FIX ATOM NAME BUG
                    if len(atom.name) == 3:
                        atom.name = ' ' + atom.name
				
		    output_line = PDB_LINE_TEMPLATE.format(record='HETATM',
							   serial=atom.serial_number,
							   atom_name=atom.name,
							   altloc=' ',
							   resname=residue.resname,
							   chain_id=residue.get_parent().id,
							   resnum=residue.get_id()[1],
							   icode=residue.get_id()[2],
							   x=atom.coord[0],
							   y=atom.coord[1],
							   z=atom.coord[2],
							   occ=1.00,
							   tfac=atom.bfactor,
							   element=atom.element,
							   charge='')
	
		    fo.write('{}\n'.format(output_line))
		
		# WRITE OUT THE PDB WITHOUT THE EXTRACTED LIGAND (AKA THE RECEPTOR)
		receptor_filename = '{}_{}_{}_{}_{}_receptor.{}'.format(pdb_noext,
																chain_id,
																resname,
																resnum,
																inscode,
																pdb_ext)
		
		# RECEPTOR WITH NO WATERS
		receptor_filename_dry = receptor_filename.replace('_receptor', '_receptor_dry')
		
		# WHOLE COMPLEX WITH NO WATERS
		complex_filename_dry = receptor_filename.replace('_receptor', '_all_dry')
		
		with open(receptor_filename, 'wb') as fo, open(receptor_filename_dry, 'wb') as dfo, open(complex_filename_dry, 'wb') as cdfo:
		
			for receptor_residue in model.get_residues():
				
				for atom in receptor_residue.child_list:
				
					# PDB OUTPUT
					# ATOM SERIALS ARE RENUMBERED FROM 1
					# ALTLOCS ARE ALWAYS BLANK
					# CHARGES ARE ALWAYS BLANK(?)
					# OCCUPANCIES ARE ALWAYS 1.00
					output_line = PDB_LINE_TEMPLATE.format(record='ATOM',
														   serial=atom.serial_number,
														   atom_name=atom.name,
														   altloc=' ',
														   resname=receptor_residue.resname,
														   chain_id=receptor_residue.get_parent().id,
														   resnum=receptor_residue.get_id()[1],
														   icode=receptor_residue.get_id()[2],
														   x=atom.coord[0],
														   y=atom.coord[1],
														   z=atom.coord[2],
														   occ=1.00,
														   tfac=atom.bfactor,
														   element=atom.element,
														   charge='')
					
					# RECEPTOR WITHOUT LIGAND BUT WITH WATERS
					if receptor_residue != residue:
						fo.write('{}\n'.format(output_line))
						
						# RECEPTOR WITHOUT LIGAND OR WATERS
						if receptor_residue.resname != 'HOH':
							dfo.write('{}\n'.format(output_line))
					
					# WHOLE COMPLEX WITHOUT WATERS
					if receptor_residue.resname != 'HOH':
						cdfo.write('{}\n'.format(output_line))
