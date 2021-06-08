import math

import matplotlib.pyplot as plt
import numpy as np
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

    def extract_indeces_ped(self, conf):
        res = self.residues
       # print(res)

        indexes = [1, int(conf+1), int(conf + res + 1), int(conf + 2 * res + 1), int(conf + 3 * res + 1), int(conf + 3 * res + 1 + res * (res - 1) / 2)]
        #RD, EN, MED_ASA, MED_RMSD, STD_DIST

        return indexes

    def global_metric(self, x, y):

        indexes = self.extract_indeces_ped(x[0])

        en = np.sum(x[indexes[1]:indexes[2]] - y[indexes[1]:indexes[2]])  # TODO: ERA SBAGLIATO
        med_asa = euclidean(x[indexes[2]:indexes[3]], y[indexes[2]:indexes[3]])
        med_rmsd = euclidean(x[indexes[3]:indexes[4]], y[indexes[3]:indexes[4]])
        std_dist = 1 - correlation(x[indexes[4]:], y[indexes[4]:])

        m = en + med_asa + med_rmsd + std_dist

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
        linkage_matrix = linkage(self.features[:, indexes[0]:], 'single', metric=self.global_metric)
        dendrogram(linkage_matrix)
        plt.show()

    def global_heatmap(self):

        dist = np.zeros((len(self.features), len(self.features)))
        for i in range(len(self.features)):
            for j in range(len(self.features)):
                dist[i, j] = self.global_metric(self.features[i], self.features[j])
                dist[j, i] = self.global_metric(self.features[i], self.features[j])

        seaborn.heatmap(dist)
        plt.show()

    def local_metric(self):
        indexes = extract_vectors_ped_feature(self.residues, self.conformations, index_slices=True)

        val_entropy = []
        val_med_asa = []
        val_med_rmsd = []
        val_med_dist = []
        val_std_dist = []

        # Scanning each element of the sequence
        for i in range(self.residues):

            entropy_i = []
            med_asa_i = []
            med_rmsd_i = []
            med_dist_i = []
            std_dist_i = []

            # Scanning each PED
            for k in range(self.features.shape[0]):
                entropy_i.append(self.features[k][indexes[2]][i])
                med_asa_i.append(self.features[k][indexes[3]][i])
                med_rmsd_i.append(self.features[k][indexes[4]][i])
                med_dist_i.append(self.features[k][indexes[5]][i])
                std_dist_i.append(self.features[k][indexes[6]][i])

            val_entropy.append(np.std(entropy_i))
            val_med_asa.append(np.std(med_asa_i))
            val_med_rmsd.append(np.std(med_rmsd_i))
            val_med_dist.append(np.mean(med_dist_i))
            val_std_dist.append(np.mean(std_dist_i))

        # Plot results same plot
        fig, axes = plt.subplots(2, 1, figsize=(24, 12))
        axes[0].set_title("Plots")
        axes[0].axhline()
        axes[0].plot(np.arange(self.residues), val_entropy, color='blue', ls='--')
        axes[0].plot(np.arange(self.residues), val_med_asa, color='red', ls='--')
        axes[0].plot(np.arange(self.residues), val_med_rmsd, color='green', ls='--')

        plt.show()

        # Plot results
        fig, axes = plt.subplots(6, 1, figsize=(24, 60))
        axes[0].set_title("ENTROPY")
        axes[0].axhline()
        axes[0].plot(np.arange(self.residues), val_entropy, ls='--')

        axes[1].set_title("ASA")
        axes[1].axhline()
        axes[1].plot(np.arange(self.residues), val_med_asa, ls='--')

        axes[2].set_title("RMDS")
        axes[2].axhline()
        axes[2].plot(np.arange(self.residues), val_med_rmsd, ls='--')

        axes[4].set_title("DIST")
        axes[4].axhline()
        axes[4].plot(np.arange(self.residues), val_med_dist, ls='--')

        axes[5].set_title("STD_DIST")
        axes[5].axhline()
        axes[5].plot(np.arange(self.residues), val_std_dist, ls='--')

        # plt.savefig('output/midterm-2/result_{}_{}.png'.format(pdb_id, chain_id), bbox_inches='tight')
        plt.show()
        fig.clear()
