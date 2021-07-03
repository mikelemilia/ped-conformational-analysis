import math
import matplotlib.pyplot as plt
import numpy as np
import pandas
import scipy.stats
import seaborn
import os

from Bio.PDB import Superimposer, PDBParser
from ModelFeatures import extract_vectors_model_feature, ModelFeatures
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import *
from utils import extract_filenames


def extract_vectors_ped_feature(residues, conformations, key=None, peds=None, features=None, indexes=False, index_slices=False):
    """
    This function allows you to extract information of the ped features from the data structure. In particular allows:
    - all rows or a specific subset of them, containing a certain feature (i.e., RD, EN, MED_ASA, etc ...)
    - the interval extremes for a certain features (i.e., RD, EN, MED_ASA, etc ...)
    - all the feature intervals as slices
    :param residues: number of residues in the model
    :param conformations: maximum number of conformations available
    :param key: the key of the feature or None if considering all of them, default: False
    :param peds: the ped id or None if considering all of them, default: False
    :param features: matrix of features or None if extracting only the indexes, default: False
    :param indexes: return (begin, end) indexes of a feature if it's True, default: False
    :param index_slices: return all the intervals of the features if it's True, default: False
    :return: begin/end, slices or features
    """

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

    begin = int(begin)
    if end is not None:
        end = int(end)

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

        self._ped_name = ped_name   # PEDxxxxx name
        self._data_folder = folder  # Folder where PED pdb files are contained
        self._folder = 'data/ped-features/'             # create ped-features always inside data folder
        self._output_folder = 'output/ped-features'     # create ped-features always inside output folder
        self._file = ped_name + '_features.csv'         # Filename for the features file
        self._file_path = self._folder + self._file     # Path for the features file

        os.makedirs(self._folder, exist_ok=True)
        os.makedirs(self._output_folder, exist_ok=True)

        # List of pdb id for the current ped (of type PEDxxxxxexxx)
        self._ped_ids = extract_filenames(self._data_folder, ped_name, 'pdb')

        # Initialization for the subsequent analysis
        self._num_residues = 0                          # Number of residues analyzed
        self._num_conformations = 0                     # Maximum number of conformations present in the ensembles
        self._models_features = []                      # Set of matrices, one for each model features file
        self._ped_features = []                         # Matrix with the features extracted from each ensemble in the rows

    def load_models_files(self):
        """
        Function to compute or extract the model features (relative to task 1) for the pdb files associated to the
        currently considered PED.
        """

        print('\n\t- Looking for model features...')
        conformations = []

        # For each ped id, apply the first task computation for the feature extraction (load if present or compute
        # if not present).
        for i, ped_id in enumerate(self._ped_ids):
            models = ModelFeatures(self._data_folder, ped_id)
            self._models_features.append(models.choice_maker(ped_id+' '))
            conformations.append(self._models_features[i].shape[0])

        # Save the number of residues and maximum number of conformations for the ensembles considered
        self._num_residues = int(self._models_features[0][0, 1])
        self._num_conformations = max(conformations)

    def choice_maker(self):
        """
        This function allows to check if the file containing the features for the selected PED has already been
        generated and stored. If so, it loads these features, otherwise it calculates the correspondent matrix
        and saves the file.
        """

        # Compute the features files corresponding to the first task
        self.load_models_files()

        if os.path.exists(self._file_path):

            # If file existing, load it
            print('\n\t- Loading PED features...')
            self.extract(self._file_path)

        else:

            # If file not existing, compute all the features and save it
            print('\n\t- Computing PED features...')
            self.compute()
            self.save(self._file_path)

        # Cast of the obtained matrix and check the presence of nan (eventually converted to 0)
        self._ped_features = np.array(self._ped_features)
        for i in range(self._ped_features.shape[0]):
            for j in range(self._ped_features.shape[1]):
                if math.isnan(self._ped_features[i, j]):
                    self._ped_features[i, j] = 0

    # ------------------------------------------------------------------------------------------------------------------
    # Ped features computation

    def compute(self):
        """
        This function allows you to insert ped features into a data structures. It's used to compare
        two different ensembles.
        """

        for k in range(len(self._models_features)):

            # Extract the features for each ensemble, cast it to a row and append it to the entire features matrix
            ped_dict = {
                'PED_ID': self._ped_ids[k],                 # PED ID
                'RD': self._models_features[k][:, 2],       # Radius of gyration for each conformation in the ensemble
                'EN': self.compute_entropy(k),              # Secondary structure entropy for each position
                'MED_ASA': self.compute_median_asa(k),      # Median solvent accessibility for each position
                'MED_RMSD': self.compute_median_rmsd(k),    # Median RMSD for each position
                'MED_DIST': self.compute_median_dist(k),    # Median distance of each pair of equivalent positions
                'STD_DIST': self.compute_std_dist(k)        # Standard deviation of the distance of each pair of equivalent positions
            }

            # Zero-padding of the radius of gyration if the ensembles contains a lower number of conformations with
            # respect to the maximum one
            if len(ped_dict['RD']) < self._num_conformations:
                ped_dict['RD'] = np.concatenate((ped_dict['RD'], np.zeros(self._num_conformations-len(ped_dict['RD']))))

            # Conversion of the dictionary to a vector and appending
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
        """
        This function computes the entropy of a specific ped using the secondary structures of its models
        :param k: model id
        :return: vector of entropies
        """

        # Extract the secondary structure vector for each k-th models features matrix
        ss = extract_vectors_model_feature(residues=self._num_residues, key='SS', features=self._models_features[k])

        # Compute entropy based on its probabilistic definition, extracting unique values and counts for each residue
        entropies = []
        for i in range(ss.shape[1]):
            unique, counts = np.unique(ss[:, i], return_counts=True)
            counts = counts[1:]
            probs = counts / np.sum(counts)
            entropy = -np.sum(probs * np.log(probs))
            entropies.append(entropy)

        return entropies

    def compute_median_asa(self, k):
        """
        This function computes the median asa of a specific ped using the asa of its models
        :param k: model id
        :return: vector of median asa
        """

        # Extract the asa vector for each k-th models features matrix
        asa = extract_vectors_model_feature(residues=self._num_residues, key='ASA', features=self._models_features[k])

        return np.median(asa, axis=0)

    def compute_median_rmsd(self, k):
        """
        This function computes the median of rmsd of a specific ped thanks to the application of the Superimposer
        :param k: model id
        :return: vector of median rmsd
        """

        super_imposer = Superimposer()
        structure_rmsd_fragments = []  # no_models X no_fragments X fragment_size
        window_size = 9

        # Get the current structure
        ped_id = self._ped_ids[k]
        structure = PDBParser(QUIET=True).get_structure(ped_id, "{}/{}.pdb".format(self._data_folder, ped_id))
        ref_model = [atom for atom in structure[0].get_atoms() if atom.get_name() == "CA"]

        for i, model in enumerate(structure):
            if i > 0:
                model_rmsd = []  # no_fragment X fragment_size
                alt_model = [atom for atom in model.get_atoms() if atom.get_name() == "CA"]  # coords of the model

                # Iterate fragments and calculate the correspondent RMSD thanks to the super_imposer operation
                for start in range(len(ref_model) - window_size):
                    end = start + window_size
                    ref_fragment = ref_model[start:end]
                    alt_fragment = alt_model[start:end]

                    # Calculate rotation/translation matrices
                    super_imposer.set_atoms(ref_fragment, alt_fragment)

                    # Rotate-translate coordinates
                    alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
                    alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
                    alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]

                    # Calculate RMSD
                    ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
                    dist = ref_fragment_coord - alt_fragment_coord
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
        """
        This function computes the median distances of a specific ped using the distance matrices of
        its models
        :param k: model id
        :return: vector of median distances
        """

        # Extract the dist vectors (linearized matrix) for each k-th models features matrix
        dist = extract_vectors_model_feature(residues=self._num_residues, key='DIST', features=self._models_features[k])

        return np.median(dist, axis=0)

    def compute_std_dist(self, k):
        """
        This function computes the standard deviation of distances in a specific ped using the distance matrices of
        its models
        :param k: ped id
        :return: vector of the standard deviations of distances
        """

        # Extract the dist vectors (linearized matrix) for each k-th models features matrix
        dist = extract_vectors_model_feature(residues=self._num_residues, key='DIST', features=self._models_features[k])

        return np.std(dist, axis=0, dtype='float64')

    def extract(self, path):
        """
        This function allows to extract features from the saved features file
        :param path: the path of features file
        """

        # Read of the file as dataframe and initialization of the matrix
        df = pandas.read_csv(path, index_col=None, header=None)
        self._ped_features = np.full((df.shape[0], df.shape[1]), None)

        # Scan of the dataframe rows and save in the data matrix
        for row in range(0, df.shape[0]):
            self._ped_features[row, :] = np.array(df.iloc[row])

    def save(self, output):
        """
        This function saves the models features matrix as a csv file
        :param output: name of output file
        """

        with open(output, 'w') as f:
            # Scan the row (models) and write each value (separated with commas) as float
            for ped in self._ped_features:
                f.write("%d" % ped[0])  # Index of the model
                for i in range(1, len(ped)):
                    f.write(",%f" % ped[i])
                f.write("\n")   # New line
        print("\t- {} saved".format(self._file))

    # ------------------------------------------------------------------------------------------------------------------
    # Ped features analysis

    def global_metric(self, x, y):
        """
        This function implements a specific metric for the ped features to compute the distance
        between two ped features vectors
            - rg (absolute difference of the means)
            - entropy (Chebyshev distance)
            - median asa (euclidean distance)
            - median rmsd (euclidean distance),
            - median distance (complementary of correlation)

        :param x: features of one ped
        :param y: features of one ped
        :return: distance between x and y
        """

        # Extract the index slices for each feature contained in the vector - to subsequently divide x and y in their
        # corresponding features vectors
        indexes = extract_vectors_ped_feature(residues=self._num_residues, conformations=self._num_conformations,
                                              index_slices=True)

        # Radius of gyration - not considering the zeros added for the padding
        x_rg = x[indexes[1]]
        x_rg_nozero = x_rg[x_rg != 0]
        y_rg = y[indexes[1]]
        y_rg_nozero = y_rg[y_rg != 0]
        # Mean of RG values for each input and absolute difference of them
        rd = np.abs(np.mean(x_rg_nozero) - np.mean(y_rg_nozero))

        # Entropy (Chebyshev distance)
        en = chebyshev(x[indexes[2]], y[indexes[2]])

        # Median ASA (euclidean distance)
        med_asa = euclidean(x[indexes[3]], y[indexes[3]])

        # Median RMSD (euclidean distance)
        med_rmsd = euclidean(x[indexes[4]], y[indexes[4]])

        # Median Distance matrix (cosine distance)
        med_dist = cosine(np.array(x[indexes[5]], dtype='float32'), np.array(y[indexes[5]], dtype='float32'))

        m = rd + en + med_asa + med_rmsd + med_dist

        return m

    def global_dendrogram(self):
        """
        This function visualizes and saves the weighted distance between pair of ensembles with a dendrogram: it is
        built from the correspondent linkage matrix (with a 'complete' approach and the global metrics)
        """

        print('\t- Plotting global dendrogram...')

        # Build linkage matrix and the corresponding labels from the IDs
        linkage_matrix = linkage(np.array(self._ped_features), 'complete', metric=self.global_metric)
        labels = np.array([pid.split('e')[1] for pid in self._ped_ids])

        # Plot the dendrogram
        dendrogram(linkage_matrix, labels=labels)
        plt.title('Global Dendrogram for {}'.format(self._ped_name))
        plt.savefig('{}/{}_dendrogram.png'.format(self._output_folder, self._ped_name))
        plt.show()

    def global_heatmap(self):
        """
        This function visualizes and saves the weighted distance between pair of ensembles with a heatmap: it is
        built exploiting the global metrics between each pair of ensembles
        """

        print('\t- Plotting global heatmap...')

        # Compute the distance matrix between each pair of ensembles
        dist = np.zeros((len(self._ped_features), len(self._ped_features)))
        for i in range(dist.shape[0]):
            for j in range(dist.shape[1]):
                dist[i, j] = self.global_metric(self._ped_features[i], self._ped_features[j])
                dist[j, i] = self.global_metric(self._ped_features[i], self._ped_features[j])

        # Generate the heatmap corresponding to the distance matrix
        labels = [pid.split('e')[1] for pid in self._ped_ids]
        seaborn.heatmap(dist, xticklabels=labels, yticklabels=labels)
        plt.title('Global Heatmap for {}'.format(self._ped_name))
        plt.savefig('{}/{}_heatmap.png'.format(self._output_folder, self._ped_name))
        plt.show()

    def local_metric(self):
        """
        This function compares features of ensembles to access the variability of each residue characteristics
         using a customized local metric
        """

        print('\t- Plotting local metric...')

        # Retrieve features (entropy, median asa, median rmsd, std dist - to be reshaped to matrix for each ensemble)
        # for each ensemble from the features matrix

        entropy = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='EN',
                                              features=self._ped_features)

        med_asa = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='MED_ASA',
                                              features=self._ped_features)

        med_rmsd = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='MED_RMSD',
                                               features=self._ped_features)

        std_dist = []
        for k in range(np.array(self._ped_features).shape[0]):
            temp_std_dist = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='STD_DIST',
                                                        features=self._ped_features, peds=k)
            # Convert dist from flatten to matrix
            std_dist_k = np.zeros((self._num_residues, self._num_residues))
            idx = np.triu_indices(self._num_residues, k=1)
            std_dist_k[idx] = temp_std_dist
            # Let the magic happen... be symmetric
            std_dist_k = std_dist_k + std_dist_k.T
            std_dist.append(std_dist_k)

        # Conversion list to array
        entropy = np.array(entropy, dtype='float64')
        med_asa = np.array(med_asa, dtype='float64')
        med_rmsd = np.array(med_rmsd, dtype='float64')
        std_dist = np.array(std_dist, dtype='float64')

        # For each feature, feature variability vectors are determined for each feature (one value for each residue)
        var_entropy = []
        var_med_asa = []
        var_med_rmsd = []
        var_std_dist = []

        # Scanning each element of the sequence and for it, consider a window for the variability computation
        window_size = 9
        for residue in range(self._num_residues):
            start = max(0, residue - window_size)
            end = min(self._num_residues - 1, residue + window_size)

            # Mean variability within the window around the residue of interest for entropy, med asa and med rmsd
            var_entropy.append(np.mean(np.std(entropy[:, start:end], axis=0)))
            var_med_asa.append(np.mean(np.std(med_asa[:, start:end], axis=0)))
            var_med_rmsd.append(np.mean(np.std(med_rmsd[:, start:end], axis=0)))

            # Determine the trimmed mean variability within the window for std dist
            current_dist = []
            for i in range(start, end):
                current_dist.append(scipy.stats.trim_mean(np.std(std_dist[:, :, i], axis=0), proportiontocut=0.2))
            var_std_dist.append(np.mean(current_dist))

        # Stack the features variability in the columns of a matrix and compute the mean by column,
        # in order to have a single variability value for each residue
        p = np.stack((var_entropy, var_med_asa, var_med_rmsd, var_std_dist), axis=1)
        val = np.mean(p, axis=1)

        # Plot results and save them
        plt.subplots(1, 1)
        plt.plot(np.arange(self._num_residues), val, color='red', ls='--')
        plt.title('Local Metric for {}'.format(self._ped_name))
        plt.savefig('{}/{}_local.png'.format(self._output_folder, self._ped_name))
        plt.show()
