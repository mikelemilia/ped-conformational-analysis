import csv
from Bio.PDB import *
from functions import *

folder = "data" #Questo deve diventare input da riga di comando
pdb_id_list = extract_filenames("data", "pdb")
# ['00020e001', '00020e002', '00020e003', '00020e004', '00020e005']

for pdb_id in pdb_id_list:

    input = "{}/{}.pdb".format(folder, pdb_id)
    print(input)
    structure = PDBParser(QUIET=True).get_structure(pdb_id, input)

    for model in structure:

        # capire come far decidere la cartella di output
        output = "data/features_files/{}/{}_{}.csv".format(pdb_id, pdb_id, str(model.id))
        features = {'N': 0, 'RG': 0, 'ASA': [], 'SS': [], 'DIST': []}

        features['N'] = len(list(model['A'].get_residues()))
        features['RG'] = radius_gyration(model['A'])

        dssp = DSSP(model, input, dssp="binx/dssp/mkdssp")  # WARNING Check the path of mkdssp

        for ss in dssp:
            # dssp index, amino acid, secondary structure, relative ASA, phi, psi,
            # NH_O_1_relidx, NH_O_1_energy, O_NH_1_relidx, O_NH_1_energy,
            # NH_O_2_relidx, NH_O_2_energy, O_NH_2_relidx, O_NH_2_energy
            features['ASA'].append(ss[3])

        dssp_dict = dict(dssp)

        features['SS'] = SSRama(model, dssp_dict)

        features['DIST'] = get_distance_matrix(list(model.get_residues()))

        save_features_csv(output, features)
