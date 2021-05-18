import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial.distance import *
from sklearn.metrics import silhouette_score
from sklearn_extra import cluster


class Clustering:

    def __init__(self, features):
        self.centroids = []

        features = np.asmatrix(features)
        self.points = features[:, 1:]

    def metrics(self, x, y):
        indexes = self.extract_indeces(x[0])

        rg = x[indexes[0]] - y[indexes[0]]
        asa = euclidean(x[indexes[1]:indexes[2]], y[indexes[1]:indexes[2]])
        ss = hamming(x[indexes[2]:indexes[3]], y[indexes[2]:indexes[3]])
        dist = 1 - correlation(x[indexes[3]:], y[indexes[3]:])

        m = rg + asa + ss + dist

        # print('Metric score: {}'.format(dist))

        return m

    def extract_indeces(self, n):
        indexes = [2, 3, int(3 + n), int(2 * n - 2 + 3)]

        return indexes

    def compute_clustering(self, k_set=range(3, 9)):

        silhouettes = []

        for k in k_set:
            kMed = cluster.KMedoids(n_clusters=k, metric=self.metrics, init='k-medoids++', max_iter=1000)
            labels = kMed.fit_predict(self.points)

            s = silhouette_score(self.points, labels, metric=self.metrics)
            silhouettes.append(s)

            # u_labels = np.unique(labels)
            # for i in u_labels:
            #     plt.scatter(points[labels == i, 0], points[labels == i, 1], label=i)
            # plt.legend()
            # plt.show()

        s_max = max(silhouettes)
        k_opt = k_set[silhouettes.index(s_max)]

        print(k_opt)
        print(s_max)

        kMed = cluster.KMedoids(n_clusters=k_opt, metric=self.metrics, init='k-medoids++', max_iter=1000)
        labels = kMed.fit_predict(self.points)
        self.centroids = kMed.medoid_indices_




    def generate_graph(self):
        # Graph
        g = nx.Graph()
        for medoid in self.centroids:
            g.add_node(medoid)

        for i in range(len(self.centroids)):
            for j in range(i + 1, len(self.centroids)):
                g.add_edge(self.centroids[i], self.centroids[j], weight=np.sum(self.points[self.centroids[i]] - self.points[self.centroids[j]]))

        nx.draw(g, with_labels=True)
        plt.show()