import numpy as np
import pandas
from Bio.PDB import Superimposer, PDBParser
from ModelFeatures import extract_vectors_model_feature, ModelFeatures
from utils import *


def extract_vectors_ped_feature(residues, conformations, key=None, features=None, peds=None, indexes=False, index_slices=False):
    begin = end = -1
    residues = int(residues)
    conformations = int(conformations)

    slices = []

    if key == 'PED_ID' or index_slices:
        begin = 0
        end = 1
        slices.append(slice(begin, end))

    if key == 'RD' or index_slices:
        begin = 1
        end = conformations + 1
        slices.append(slice(begin, end))

    if key == 'EN' or index_slices:
        begin = conformations + 1
        end = conformations + residues + 1
        slices.append(slice(begin, end))

    if key == 'MED_ASA' or index_slices:
        begin = conformations + residues + 1
        end = conformations + 2 * residues + 1
        slices.append(slice(begin, end))

    if key == 'MED_RMSD' or index_slices:
        begin = conformations + 2 * residues + 1
        end = conformations + 3 * residues + 1
        slices.append(slice(begin, end))

    if key == 'MED_DIST' or index_slices:
        begin = conformations + 3 * residues + 1
        end = int(conformations + 3 * residues + 1 + residues * (residues - 1) / 2)
        slices.append(slice(begin, end))

    if key == 'STD_DIST' or index_slices:
        begin = int(conformations + 3 * residues + 1 + residues * (residues - 1) / 2)
        end = None
        slices.append(slice(begin, end))

    if begin == -1:
        return None

    if index_slices:
        return slices

    if indexes is True or features is None:
        return begin, end

    if peds is None:
        return features[:, begin:end]
    else:
        if isinstance(peds, int):
            return np.array(features[peds][begin:end])
        else:
            return features[peds, begin:end]


class PedFeatures:

    def __init__(self, folder, ped_name):
        # Extract all the files PEDxxxxxexxx_features.csv
        model_folder = folder + '/model_features/'
        ped_names = extract_filenames(model_folder, ped_name, 'csv')

        # Extract PEDxxxxxexxx from filenames
        self._ped_ids = []
        for name in ped_names:
            self._ped_ids.append(name.split('_')[0])

        # PEDxxxxx name
        self._ped_name = ped_name

        # Build all the paths to the features files and extract them
        self._models_features = []
        confs = []
        for i in range(len(self._ped_ids)):
            path = model_folder + ped_names[i] + '.csv'
            models = ModelFeatures('data', self._ped_name)
            self._models_features.append(models.extract(path))
            confs.append(self._models_features[i].shape[0])

        # Prepare the variables for the subsequent analysis
        self._num_residues = int(self._models_features[0][0, 1])
        self._num_conformations = max(confs)
        self._ped_features = []

        # Folders
        self._data_folder = folder
        self._folder = folder + '/ped_features/'
        os.makedirs(self._folder, exist_ok=True)
        self._file = ped_name + '_features.csv'

    def choice_maker(self):

        if os.path.exists(self._folder + self._file):
            print('\nLoading features comparison...')
            self._ped_features = self.extract(self._folder + self._file)
        else:
            print('\nComparing features...')
            self.compare()
            self.save(self._folder + self._file)

        return self._ped_features

    def compare(self):

        for k in range(len(self._models_features)):

            ped_dict = {
                'PED_ID': self._ped_ids[k],
                'RD': self._models_features[k][:, 2],
                'EN': self.compute_entropy(k),
                'MED_ASA': self.compute_median_asa(k),
                'MED_RMSD': self.compute_median_rmsd(k),
                'MED_DIST': self.compute_median_dist(k),
                'STD_DIST': self.compute_std_dist(k)
            }

            if len(ped_dict['RD']) < self._num_conformations:
                ped_dict['RD'] = np.concatenate((ped_dict['RD'], np.zeros(self._num_conformations-len(ped_dict['RD']))))

            x = []
            for key in ped_dict.keys():
                if key in ['PED_ID']:
                    x.append(k)
                else:
                    for e in ped_dict[key]:
                        x.append(e)
            self._ped_features.append(x)

        return self._ped_features

    def compute_entropy(self, k):

        ss = extract_vectors_model_feature(residues=self._num_residues, key='SS', features=self._models_features[k])
        entropies = []
        for i in range(ss.shape[1]):
            unique, counts = np.unique(ss[:, i], return_counts=True)
            counts = counts[1:]
            probs = counts / np.sum(counts)
            entropy = -np.sum(probs * np.log(probs))
            entropies.append(entropy)
        return entropies

    def compute_median_asa(self, k):

        asa = extract_vectors_model_feature(residues=self._num_residues, key='ASA', features=self._models_features[k])
        return np.median(asa, axis=0)

    def compute_median_rmsd(self, k):

        super_imposer = Superimposer()
        structure_rmsd_fragments = []  # RMSD, no_models X no_fragments X fragment_size
        window_size = 9

        ped_id = self._ped_ids[k]
        structure = PDBParser(QUIET=True).get_structure(ped_id, "{}/{}.pdb".format(self._data_folder, ped_id))
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
        structure_rmsd_fragments = np.pad(structure_rmsd_fragments, ((0, 0), (0, structure_rmsd_fragments.shape[0])))

        # Roll the fragments one by one (add heading zeros)
        for i, row in enumerate(structure_rmsd_fragments):
            structure_rmsd_fragments[i] = np.roll(row, i)

        # Calculate average along columns of overlapping fragments (average RMSD per residue)
        return np.average(structure_rmsd_fragments, axis=0)

    def compute_median_dist(self, k):

        dist = extract_vectors_model_feature(residues=self._num_residues, key='DIST', features=self._models_features[k])
        return np.median(dist, axis=0)

    def compute_std_dist(self, k):

        dist = extract_vectors_model_feature(residues=self._num_residues, key='DIST', features=self._models_features[k])
        return np.std(dist, axis=0, dtype='float64')

    def save(self, output):

        with open(output, 'w') as f:
            for ped in self._ped_features:
                f.write("%d" % ped[0])  # index of the model
                for i in range(1, len(ped)):
                    f.write(",%f" % ped[i])
                f.write("\n")
        print("{}.csv saved".format(self._ped_name))

    def extract(self, path):

        df = pandas.read_csv(path, index_col=None, header=None)

        self._ped_features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self._ped_features[row, :] = np.array(df.iloc[row])

        return self._ped_features

    @property
    def num_residues(self):
        return self._num_residues

    @property
    def num_conformations(self):
        return self._num_conformations
