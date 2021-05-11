from Bio.PDB import *
from functions import *

pdb_id = '00020e001'
structure = PDBParser(QUIET=True).get_structure(pdb_id, "data/PED{}.pdb".format(pdb_id))

print(radius_gyration(structure[0]['A']))


