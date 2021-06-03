import math
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib.pyplot as plt
import seaborn
from PedFeatures import *


class PedComparison:

    def __init__(self, features, conformations, residues):
        self.features = features
        self.conformations = conformations
        self.residues = residues

        for i in range(self.features.shape[0]):
            for j in range(self.features.shape[1]):
                if math.isnan(self.features[i, j]):
                    self.features[i, j] = 0

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
        linkage_matrix = linkage(self.features[:, indexes[0]:], 'single')  # qui cambiare la funzione di metric
        dendrogram(linkage_matrix)
        plt.show()
