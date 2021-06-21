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


def extract_vectors_ped_feature(residues, conformations, key=None, features=None, peds=None, indexes=False, index_slices=False):
    """
    This function allows you to extract information of the model features from the data structure.
    :param residues: number of residues in the model
    :param key: the key of the feature or None
    :param models: the id of model or None
    :param features: None
    :param indexes: only ruturn begin, end of same feature if it's True, default: False
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
        # Extract all the files PEDxxxxxexxx_features.csv
        model_folder = folder + '/model_features/'
        ped_names = extract_filenames(model_folder, ped_name, ['csv'])

        # Extract PEDxxxxxexxx from filenames
        self._ped_ids = []
        for name in ped_names:
            self._ped_ids.append(name.split('_')[0])

        if len(self._ped_ids) != 0:

            # PEDxxxxx name
            self._ped_name = ped_name

            # Build all the paths to the features files and extract them
            self._models_features = []

            conformations = []
            for i in range(len(self._ped_ids)):
                path = model_folder + ped_names[i] + '.csv'
                models = ModelFeatures('data', self._ped_name)
                self._models_features.append(models.extract(path))
                conformations.append(self._models_features[i].shape[0])

            # Prepare the variables for the subsequent analysis
            self._num_residues = int(self._models_features[0][0, 1])
            self._num_conformations = max(conformations)
            self._ped_features = []

            # Folders
            self._data_folder = folder
            self._folder = folder + '/ped_features/'
            os.makedirs(self._folder, exist_ok=True)
            self._file = ped_name + '_features.csv'

    def choice_maker(self):

        if len(self._ped_ids) == 0:
            return -1

        if os.path.exists(self._folder + self._file):
            print('\n\t- Loading features for comparison...')
            self.extract(self._folder + self._file)
        else:
            print('\n\t- Computing features for comparison...')
            self.compare()
            self.save(self._folder + self._file)

        self._ped_features = np.array(self._ped_features)
        for i in range(self._ped_features.shape[0]):
            for j in range(self._ped_features.shape[1]):
                if math.isnan(self._ped_features[i, j]):
                    self._ped_features[i, j] = 0

        return 0

    def compare(self):
        """
        This function allows you to insert ped features into a data structures.
        It's used to compare two different ped.
        :return: ped features
        """

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
        """
        This function computes the entropy of a specific ped using the secondary structures of
        its models
        :param k: ped id 
        :return: vector of entropies
        """

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
        """
        This function computes the median asa of a specific ped using the asa of
        its models
        :param k: ped id
        :return: vector of median asa
        """

        asa = extract_vectors_model_feature(residues=self._num_residues, key='ASA', features=self._models_features[k])
        return np.median(asa, axis=0)

    def compute_median_rmsd(self, k):
        """
        This function computes the median of rmsd of a specific ped
        :param k: ped id
        :return: vector of median rmsd
        """

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
        """
        This function computes the median distances of a specific ped using the distance matrices of
        its models
        :param k: ped id
        :return: vector of median distances
        """

        dist = extract_vectors_model_feature(residues=self._num_residues, key='DIST', features=self._models_features[k])
        return np.median(dist, axis=0)

    def compute_std_dist(self, k):
        """
        This function computes the standard deviation of distances in a specific ped using the distance matrices of
        its models
        :param k: ped id
        :return: vector of the standard deviations of distaces
        """

        dist = extract_vectors_model_feature(residues=self._num_residues, key='DIST', features=self._models_features[k])
        return np.std(dist, axis=0, dtype='float64')

    def save(self, output):
        """
        This function saves the peds features as a file
        :param output: name of output file
        :return: True or False, If the output file of all the models features is saved correctly
        """

        with open(output, 'w') as f:
            for ped in self._ped_features:
                f.write("%d" % ped[0])  # index of the model
                for i in range(1, len(ped)):
                    f.write(",%f" % ped[i])
                f.write("\n")
        print("\t- {} saved".format(self._file))

    def extract(self, path):
        """
        This function allows you to extract ped features
        :param path: the name of ped file
        :return: ped features
        """

        df = pandas.read_csv(path, index_col=None, header=None)

        self._ped_features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self._ped_features[row, :] = np.array(df.iloc[row])

        return self._ped_features

    def global_metric(self, x, y):
        """
        This function use a specific metric and all the ped features to compute the distance
        between two input points
        Entropy (std), median asa (euclidean distance), median rmsd (euclidean distance),
        median distance (correlation)
        :param x: features of one ped
        :param y: features of one ped
        :return: distance between x and y
        """

        indexes = extract_vectors_ped_feature(residues=self._num_residues, conformations=self._num_conformations,
                                              index_slices=True)

        x_rg = x[indexes[1]]
        x_rg_nozero = x_rg[x_rg != 0]
        y_rg = y[indexes[1]]
        y_rg_nozero = y_rg[y_rg != 0]
        rd = np.abs(np.mean(x_rg_nozero) - np.mean(y_rg_nozero))

        en = np.abs(np.sum(x[indexes[2]] - y[indexes[2]]))
        med_asa = euclidean(x[indexes[3]], y[indexes[3]])
        med_rmsd = euclidean(x[indexes[4]], y[indexes[4]])
        med_dist = 1 - correlation(np.array(x[indexes[5]], dtype='float32'), np.array(y[indexes[5]], dtype='float32'))

        m = rd + en + med_asa + med_rmsd + med_dist

        return m

    def distance_matrix_med_rmsd_peds(self):
        """
        This function visualizes the distance matrix heatmap using only RMSD between pair of peds
        :return: Plot heatmap with respect RMSD
        """

        print('\t- Plotting heatmap (with only RMSD)...')

        medians = extract_vectors_ped_feature(self._num_residues, self._num_conformations, 'MED_RMSD', features=self._ped_features)
        dists = np.zeros((medians.shape[0], medians.shape[0]))

        for i in range(medians.shape[0]):
            for j in range(i + 1, medians.shape[0]):
                num_pairs = (self._num_residues * (self._num_residues - 1)) / 2
                dists[i, j] = np.sqrt(1 / num_pairs * np.sum((medians[i] - medians[j]) ** 2, axis=0))
                dists[j, i] = np.sqrt(1 / num_pairs * np.sum((medians[i] - medians[j]) ** 2, axis=0))

        seaborn.heatmap(dists)
        plt.show()

        return dists

    def global_dendrogram(self):
        """
        This function visualizes the weighted distance between pair of peds with respect global metric
        :return: Dendrogram of peds, save it into img file
        """

        print('\t- Plotting global dendrogram...')

        linkage_matrix = linkage(np.array(self._ped_features), 'complete', metric=self.global_metric)  # TODO capire se è realmente necessario il cast, in teoria no
        dendrogram(linkage_matrix)
        plt.title('Global Dendrogram for {}'.format(self._ped_name))
        plt.savefig('output/plot/{}_dendrogram.png'.format(self._ped_name))
        plt.show()

    def global_heatmap(self):
        """
        This function allows you to visualize the pairwise difference of peds
        using a customized global metric.
        :return: Plot heatmap, save it into a img file
        """

        print('\t- Plotting global heatmap...')

        dist = np.zeros((len(self._ped_features), len(self._ped_features)))

        for i in range(dist.shape[0]):
            for j in range(dist.shape[1]):
                dist[i, j] = self.global_metric(self._ped_features[i], self._ped_features[j])
                dist[j, i] = self.global_metric(self._ped_features[i], self._ped_features[j])

        seaborn.heatmap(dist)
        plt.title('Global Heatmap for {}'.format(self._ped_name))
        plt.savefig('output/plot/{}_heatmap.png'.format(self._ped_name))
        plt.show()

    def local_metric(self):
        """
        This function compares features of one ped using a customized local metric
        :return: Plots of the ped features
        """

        print('\t- Plotting local metric...')

        # Retrieve features each PED

        entropy = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='EN',
                                              features=self._ped_features)
        med_asa = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='MED_ASA',
                                              features=self._ped_features)
        med_rmsd = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='MED_RMSD',
                                               features=self._ped_features)
        med_dist = []
        std_dist = []

        # Convert dist from flatten to matrix
        for k in range(np.array(self._ped_features).shape[0]): # TODO capire se è realmente necessario il cast, in teoria no
            temp_med_dist = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='MED_DIST', features=self._ped_features, peds=k)
            temp_std_dist = extract_vectors_ped_feature(self._num_residues, self._num_conformations, key='STD_DIST', features=self._ped_features, peds=k)

            med_dist_k = np.zeros((self._num_residues, self._num_residues))
            idx = np.triu_indices(self._num_residues, k=1)
            med_dist_k[idx] = temp_med_dist

            # Let the magic happen... be symmetric
            med_dist_k = med_dist_k + med_dist_k.T

            std_dist_k = np.zeros((self._num_residues, self._num_residues))
            idx = np.triu_indices(self._num_residues, k=1)
            std_dist_k[idx] = temp_std_dist

            # Let the magic happen... be symmetric
            std_dist_k = std_dist_k + std_dist_k.T

            med_dist.append(med_dist_k)
            std_dist.append(std_dist_k)

        # Conversion list to array
        med_dist = np.array(med_dist)
        std_dist = np.array(std_dist)

        # Total values for each residue
        total_entropy = []
        total_med_asa = []
        total_med_rmsd = []
        total_med_dist = []
        total_std_dist = []

        # Scanning each element of the sequence
        for i in range(self._num_residues):
            total_entropy.append(np.std(entropy[:, i]))
            total_med_asa.append(np.std(med_asa[:, i]))  # TODO check how to normalize ASA
            total_med_rmsd.append(np.std(med_rmsd[:, i]))

            # Compute std deviations for residue i
            total_std_dist.append(scipy.stats.trim_mean(np.std(std_dist[:, :, i], axis=0), proportiontocut=0.2))

        p = np.stack((total_entropy, total_med_asa, total_med_rmsd, total_std_dist), axis=1)
        val = np.mean(p, axis=1)

        # Plot results same plot
        fig, axes = plt.subplots(1, 1, figsize=(24, 12))
        # axes[0].set_title("Plots")
        # axes[0].axhline()
        # plt.plot(np.arange(self.residues), total_entropy, color='blue', ls='--')
        # plt.plot(np.arange(self.residues), total_med_asa, color='red', ls='--')
        # plt.plot(np.arange(self.residues), total_med_rmsd, color='green', ls='--')
        # # plt.plot(np.arange(self.residues), total_med_dist, color='orange', ls='--')
        # plt.plot(np.arange(self.residues), total_std_dist, color='pink', ls='--')

        plt.plot(np.arange(self._num_residues), val, color='red', ls='--')
        plt.title('Local Metric for {}'.format(self._ped_name))
        plt.savefig('output/plot/{}_local.png'.format(self._ped_name))
        plt.show()
