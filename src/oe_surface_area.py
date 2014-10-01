import os
import sys

from openeye.oechem import *
from openeye.oespicoli import *

if __name__ == '__main__':
    
    pdb_file = sys.argv[1]
    pdb_path_noext = os.path.splitext(pdb_file)[0]
    pdb_noext, pdb_ext = os.path.splitext(os.path.split(pdb_file)[1])
    
    surface_area_outfile = pdb_path_noext + '.sasa'
    volume_outfile = pdb_path_noext + '.vol'
    per_atom_sasa_pdb = pdb_path_noext + '.atom_sasa'
    
    # READ IN MOLECULE
    ifs = oemolistream()
    #ifs.SetFormat(OEFormat_PDB)
    
    mol = OEMol()
    
    if ifs.open(pdb_file):
        
        OEReadMolecule(ifs, mol)
        
    else:
        OEThrow.Fatal("Unable to open file {}".format(pdb_file))
    
    # MAKE ACCESSIBLE SURFACE
    OEAssignBondiVdWRadii(mol)
    
    surf = OESurface()
    OEMakeAccessibleSurface(surf, mol)
    
    mol_surface_area = OESurfaceArea(surf)
    mol_surface_volume = OESurfaceVolume(surf)
    
    # WRITE OUT SURFACE PROPERTIES
    with open(surface_area_outfile, 'wb') as fo:
        fo.write(str(mol_surface_area))
        
    with open(volume_outfile, 'wb') as fo:
        fo.write(str(mol_surface_volume))
    
    # SPLIT THE SURFACE BY ATOMS
    split_surf = OESurface()
    OESplitSurfaceByAtoms(split_surf, surf, mol)
    
    # ADAPTED FROM `http://docs.eyesopen.com/toolkits/spicoli/python/OESpicoliFunctions.html#OESpicoli::OECalculateTriangleAreas`
    areas = OEFloatArray(split_surf.GetNumTriangles())
    OECalculateTriangleAreas(split_surf, areas)

    atomareas = [0.0] * mol.GetMaxAtomIdx()
    
    for i in xrange(split_surf.GetNumTriangles()):
        tri = split_surf.GetTriangle(i)
        
        a1 = split_surf.GetAtomsElement(tri[0])
        a2 = split_surf.GetAtomsElement(tri[1])
        a3 = split_surf.GetAtomsElement(tri[2])
        
        atomareas[a1] += areas[i]/3.0
        atomareas[a2] += areas[i]/3.0
        atomareas[a3] += areas[i]/3.0
    
    # ASSIGN PER-ATOM SURFACE AREA TO B-FACTOR
    for atom in mol.GetAtoms():
        
        residue = OEAtomGetResidue(atom)
        residue.SetBFactor(atomareas[atom.GetIdx()])
    
    # WRITE OUT A PDB WITH THE B-FACTORS SET TO PER-ATOM SASA
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_PDB)
    
    if ofs.open(per_atom_sasa_pdb):
        OEWriteMolecule(ofs, mol)
    else:
        OEThrow.Fatal('Unable to write per-atom SASA pdb: {}'.format(per_atom_sasa_pdb))
    