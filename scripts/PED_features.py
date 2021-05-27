import numpy as np
import pandas
from Bio.PDB import Superimposer, PDBParser
from Model_features import extract_vectors_feature, Model_Features
from utils import *

class PED_features():

    def __init__(self, path):
        ped_names = extract_filenames(path, 'csv')
        self.ped_ids = []
        for ped_name in ped_names:
            self.ped_ids.append(ped_name.removesuffix('_features'))
        self.ped_name = (self.ped_ids[0].split('e'))[0]

        paths = []
        for i in range(len(self.ped_ids)):
            paths.append('data/features_files/' + ped_names[i] + '.csv') # TODO: sistemare input

        self.models_features = []
        for path in paths:
            models = Model_Features('data', self.ped_name)
            self.models_features.append(models.extract(path))

        self.num_residues = int(self.models_features[0][0, 1])

        self.ped_features = []

    def compare(self):

        for k in range(len(self.models_features)):

            ped_dict = {
                'PED_ID': self.ped_ids[k],
                'RD': self.models_features[k][:, 2],
                'EN': self.compute_entropy(k),
                'MED_ASA': self.compute_median_asa(k),
                'MED_RMSD': self.compute_median_rmsd(k),
                'MED_DIST': self.compute_median_dist(k),
                'STD_DIST': self.compute_std_dist(k)
            }

            x = []
            for key in ped_dict.keys():
                if key in ['PED_ID']:
                    x.append(k)
                else:
                    for e in ped_dict[key]:
                        x.append(e)
            self.ped_features.append(x)

        return self.ped_features

    def compute_entropy(self, k):
        ss = extract_vectors_feature(self.models_features[k], 'SS')
        entropies = []
        for i in range(ss.shape[1]):
            unique, counts = np.unique(ss[:, i], return_counts=True)
            probs = counts / np.sum(counts)
            entropy = -np.sum(probs * np.log(probs))
            entropies.append(entropy)
        return entropies

    def compute_median_asa(self, k):
        asa = extract_vectors_feature(self.models_features[k], 'ASA')
        return np.median(asa, axis=0)

    def compute_median_rmsd(self, k):
        super_imposer = Superimposer()
        structure_rmsd_fragments = []  # RMSD, no_models X no_fragments X fragment_size
        window_size = 9

        ped_id = self.ped_ids[k]
        structure = PDBParser(QUIET=True).get_structure(ped_id, "data/{}.pdb".format(ped_id))
        ref_model = [atom for atom in structure[0].get_atoms() if atom.get_name() == "CA"]  # TODO: capire!!!

        for i, model in enumerate(structure):
            if i > 0:
                model_rmsd = []  # RMSD, no_fragment X fragment_size
                alt_model = [atom for atom in model.get_atoms() if atom.get_name() == "CA"]  # coords of the model

                # Iterate fragments
                for start in range(len(ref_model) - window_size):
                    end = start + window_size
                    ref_fragment = ref_model[start:end]
                    alt_fragment = alt_model[start:end]

                    # Calculate rotation/translation matrices
                    super_imposer.set_atoms(ref_fragment, alt_fragment)
                    # print(super_imposer.rms, super_imposer.rotran)

                    # Rotate-translate coordinates
                    alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
                    alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
                    alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]

                    # Calculate RMSD: https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
                    ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
                    dist = ref_fragment_coord - alt_fragment_coord
                    # rmsd_fragment = np.sqrt(np.sum(dist * dist) / window_size)  # Total RMSD of the fragment. Identical to super_imposer.rms
                    rmsd_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment

                    model_rmsd.append(rmsd_res)

                structure_rmsd_fragments.append(model_rmsd)

        # Calculate average RMSD per position
        structure_rmsd_fragments = np.array(structure_rmsd_fragments)
        structure_rmsd_fragments = np.average(structure_rmsd_fragments, axis=0)
        structure_rmsd_fragments = np.pad(structure_rmsd_fragments,
                                          ((0, 0), (0, structure_rmsd_fragments.shape[0])))

        # Roll the fragments one by one (add heading zeros)
        for i, row in enumerate(structure_rmsd_fragments):
            structure_rmsd_fragments[i] = np.roll(row, i)

        # Calculate average along columns of overlapping fragments (average RMSD per residue)
        return np.average(structure_rmsd_fragments, axis=0)

    def compute_median_dist(self, k):
        dist = extract_vectors_feature(self.models_features[k], 'DIST')
        return np.median(dist, axis=0)

    def compute_std_dist(self, k):
        dist = extract_vectors_feature(self.models_features[k], 'DIST')
        return np.std(dist, axis=0, dtype='float64')

    def save(self, output):

        with open(output, 'w') as f:
            for ped in self.ped_features:
                f.write("%d" % ped[0])  # index of the model
                for i in range(1, len(ped)):
                    f.write(",%f" % ped[i])
                f.write("\n")
        print("{}_features.csv saved".format(self.ped_name))

    def extract(self, path):

        df = pandas.read_csv(path, index_col=None, header=None)

        self.ped_features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self.ped_features[row, :] = np.array(df.iloc[row])

        return self.ped_features

    def get_ped_name(self):
        return self.ped_name