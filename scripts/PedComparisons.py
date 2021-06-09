import math

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import seaborn
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import *

from PedFeatures import *


class PedComparison:

    def __init__(self, features, conformations, residues):
        self.features = np.array(features)
        self.conformations = conformations
        self.residues = residues

        for i in range(self.features.shape[0]):
            for j in range(self.features.shape[1]):
                if math.isnan(self.features[i, j]):
                    self.features[i, j] = 0

    def global_metric(self, x, y):

        indexes = extract_vectors_ped_feature(residues=self.residues, conformations=self.conformations,
                                              index_slices=True)

        rd = np.abs(np.mean(x[indexes[1]]) - np.mean(y[indexes[1]]))
        en = np.abs(np.sum(x[indexes[2]] - y[indexes[2]]))
        med_asa = euclidean(x[indexes[3]], y[indexes[3]])
        med_rmsd = euclidean(x[indexes[4]], y[indexes[4]])
        med_dist = 1 - correlation(x[indexes[5]], y[indexes[5]])

        m = rd + en + med_asa + med_rmsd + med_dist

        return m

    def distance_matrix_med_rmsd_peds(self):

        medians = extract_vectors_ped_feature(self.residues, self.conformations, 'MED_RMSD', features=self.features)
        dists = np.zeros((medians.shape[0], medians.shape[0]))

        for i in range(medians.shape[0]):
            for j in range(i + 1, medians.shape[0]):
                num_pairs = (self.residues * (self.residues - 1)) / 2
                dists[i, j] = np.sqrt(1 / num_pairs * np.sum((medians[i] - medians[j]) ** 2, axis=0))
                dists[j, i] = np.sqrt(1 / num_pairs * np.sum((medians[i] - medians[j]) ** 2, axis=0))

        seaborn.heatmap(dists)
        plt.show()

        return dists

    def global_dendrogram(self):

        indexes = extract_vectors_ped_feature(self.residues, self.conformations, 'EN', indexes=True)
        linkage_matrix = linkage(self.features[:, indexes[0]:], 'complete', metric=self.global_metric)
        dendrogram(linkage_matrix)
        plt.show()

    def global_heatmap(self):

        indexes = extract_vectors_ped_feature(self.residues, self.conformations, 'EN', indexes=True)

        dist = np.zeros((len(self.features[:, indexes[0]:]), len(self.features[:, indexes[0]:])))
        for i in range(dist.shape[0]):
            for j in range(dist.shape[1]):
                dist[i, j] = self.global_metric(self.features[i, indexes[0]:], self.features[j, indexes[0]:])
                dist[j, i] = self.global_metric(self.features[j, indexes[0]:], self.features[i, indexes[0]:])
        seaborn.heatmap(dist)
        plt.show()

    def local_metric(self):

        # Retrieve features each PED

        entropy = extract_vectors_ped_feature(self.residues, self.conformations, key='EN', features=self.features)
        med_asa = extract_vectors_ped_feature(self.residues, self.conformations, key='MED_ASA', features=self.features)
        med_rmsd = extract_vectors_ped_feature(self.residues, self.conformations, key='MED_RMSD', features=self.features)
        med_dist = []
        std_dist = []

        # Convert dist from flatten to matrix
        for k in range(self.features.shape[0]):
            temp_med_dist = extract_vectors_ped_feature(self.residues, self.conformations, key='MED_DIST', features=self.features, peds=k)
            temp_std_dist = extract_vectors_ped_feature(self.residues, self.conformations, key='STD_DIST', features=self.features, peds=k)

            med_dist_k = np.zeros((self.residues, self.residues))
            idx = np.triu_indices(self.residues, k=1)
            med_dist_k[idx] = temp_med_dist

            # Let the magic happen... be symmetric
            med_dist_k = med_dist_k + med_dist_k.T

            std_dist_k = np.zeros((self.residues, self.residues))
            idx = np.triu_indices(self.residues, k=1)
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
        for i in range(self.residues):
            total_entropy.append(np.std(entropy[:, i]))
            total_med_asa.append(np.std(med_asa[:, i]))  # TODO check how to normalize ASA
            total_med_rmsd.append(np.std(med_rmsd[:, i]))

            # Compute std deviations for residue i
            total_std_dist.append(scipy.stats.trim_mean(np.std(std_dist[:, :, i], axis=0), proportiontocut=0.2))

        p = np.stack((total_entropy, total_med_asa, total_med_rmsd, total_std_dist), axis=1)
        val = np.mean(p, axis=1)
        print(p.shape)
        print(val.shape)

        # Plot results same plot
        fig, axes = plt.subplots(1, 1, figsize=(24, 12))
        # axes[0].set_title("Plots")
        # axes[0].axhline()
        # plt.plot(np.arange(self.residues), total_entropy, color='blue', ls='--')
        # plt.plot(np.arange(self.residues), total_med_asa, color='red', ls='--')
        # plt.plot(np.arange(self.residues), total_med_rmsd, color='green', ls='--')
        # # plt.plot(np.arange(self.residues), total_med_dist, color='orange', ls='--')
        # plt.plot(np.arange(self.residues), total_std_dist, color='pink', ls='--')

        plt.plot(np.arange(self.residues), val, color='red', ls='--')

        plt.show()

        # Plot results
        # fig, axes = plt.subplots(6, 1, figsize=(24, 60))
        # axes[0].set_title("ENTROPY")
        # axes[0].axhline()
        # axes[0].plot(np.arange(self.residues), total_entropy, ls='--')
        #
        # axes[1].set_title("ASA")
        # axes[1].axhline()
        # axes[1].plot(np.arange(self.residues), total_med_asa, ls='--')
        #
        # axes[2].set_title("RMDS")
        # axes[2].axhline()
        # axes[2].plot(np.arange(self.residues), total_med_rmsd, ls='--')
        #
        # axes[4].set_title("DIST")
        # axes[4].axhline()
        # axes[4].plot(np.arange(self.residues), total_med_dist, ls='--')
        #
        # axes[5].set_title("STD_DIST")
        # axes[5].axhline()
        # axes[5].plot(np.arange(self.residues), total_std_dist, ls='--')
        #
        # # plt.savefig('output/midterm-2/result_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')
        # plt.show()
        # fig.clear()
