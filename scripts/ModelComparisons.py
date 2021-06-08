# import networkx as nx
# import numpy as np
# from Bio.PDB import PDBList, Superimposer, Selection
# from Bio.PDB.PDBParser import PDBParser
# from matplotlib import pyplot as plt
# from scipy.spatial.distance import *
# from sklearn.metrics import silhouette_score
# from sklearn_extra import cluster
# # import pymol
# # from pymol import cmd
# from matplotlib import cm, colors
#
#
# class ModelComparisons:
#
#     def __init__(self, identifier, features):
#
#         self.id = identifier
#         self.centroids = []
#         self.labels = []
#         self.features = np.asmatrix(features, dtype='float64')
#         self.points = features[:, 1:]
#
#     def metrics(self, x, y):
#
#         indexes = self.extract_indeces(x[0])
#
#         rg = x[indexes[0]] - y[indexes[0]] # TODO: SISTEMARE
#         asa = euclidean(x[indexes[1]:indexes[2]], y[indexes[1]:indexes[2]])
#         ss = hamming(x[indexes[2]:indexes[3]], y[indexes[2]:indexes[3]])
#         dist = 1 - correlation(x[indexes[3]:], y[indexes[3]:])
#
#         m = rg + asa + ss + dist
#
#         # print('Metric score: {}'.format(dist))
#
#         return m
#
#     def extract_indeces(self, n):
#
#         indexes = [2, 3, int(3 + n), int(2 * n - 2 + 3)]
#
#         return indexes
#
#     def compute_clustering(self, k_set=range(3, 9)):
#
#         silhouettes = []
#
#         for k in k_set:
#             kMed = cluster.KMedoids(n_clusters=k, metric=self.metrics, init='k-medoids++', max_iter=1000)
#             labels = kMed.fit_predict(self.points)
#
#             s = silhouette_score(self.points, labels, metric=self.metrics)
#             silhouettes.append(s)
#
#             # u_labels = np.unique(labels)
#             # for i in u_labels:
#             #     plt.scatter(points[labels == i, 0], points[labels == i, 1], label=i)
#             # plt.legend()
#             # plt.show()
#
#         s_max = max(silhouettes)
#         k_opt = k_set[silhouettes.index(s_max)]
#
#         print('Number of representative conformations: {}'.format(k_opt))
#         print('Correspondent silhouette value: {}'.format(s_max))
#
#         kMed = cluster.KMedoids(n_clusters=k_opt, metric=self.metrics, init='k-medoids++', max_iter=1000)
#         self.labels = kMed.fit_predict(self.points)
#         self.centroids = kMed.medoid_indices_
#         print('Indexes of the representative conformations: {}'.format(self.centroids))
#
#     def generate_graph(self):
#
#         # Graph
#         g = nx.Graph()
#         for medoid in self.centroids:
#             g.add_node(medoid)
#
#         for i in range(len(self.centroids)):
#             for j in range(i + 1, len(self.centroids)):
#                 g.add_edge(self.centroids[i], self.centroids[j],
#                            weight=np.sum(self.points[self.centroids[i]] - self.points[self.centroids[j]]))
#
#         options = {'node_color': 'orange', 'node_size': 700, 'width': 2}
#         nx.draw_spectral(g, with_labels=True, **options)
#         path = "output/{}_graph.png".format(self.id)
#         plt.savefig(path)
#         plt.show()
#         return g

    # def generate_pymol_img(self, g):
    #     pymol.finish_launching()  # Open Pymol
    #     # Load the structure conformations
    #     h = g.nodes().value()
    #     print(h)
    #     for j in g.nodes().value():
    #         structure = PDBParser(QUIET=True).get_structure(self.id)
    #         # Superimpose all models to the first model, fragment-by-fragment (sliding window)
    #         super_imposer = Superimposer()
    #
    #         structure_feature_fragments = []
    #         window_size = 9
    #         ref_model = [atom for atom in structure[h].get_atoms() if atom.get_name() == "CA"]  # CA of the first model
    #         ref_features = self.features[h]
    #         model_features = []
    #         alt_model = [atom for atom in structure[j].get_atoms() if atom.get_name() == "CA"]  # coords of the model
    #         alt_features = self.features[j]
    #
    #         # Iterate fragments
    #         for start in range(len(ref_model) - window_size):
    #             end = start + window_size
    #             ref_fragment = ref_model[start:end]
    #             alt_fragment = alt_model[start:end]
    #
    #             # Calculate rotation/translation matrices
    #             super_imposer.set_atoms(ref_fragment, alt_fragment)
    #
    #             # Rotate-translate coordinates
    #             # alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
    #             # alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
    #             # alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]
    #
    #             #features of structure
    #             #ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
    #             # dist = np.diff(ref_features,alt_features)
    #             # feat_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment
    #             # model_features.append(feat_res)
    #         dist = np.diff(ref_features,alt_features)
    #         feat_res = np.sqrt(np.sum(dist * dist, axis=1))
    #         structure_feature_fragments.append(feat_res)
    #
    #
    #         structure_feature_fragments = np.array(structure_feature_fragments)  # no_models X no_fragments X fragment_size
    #
    #         structure_feature_fragments = np.average(structure_feature_fragments, axis=0)  # no_fragments X fragment_size
    #         # Pad with right zeros to reach the sequence length (no_fragments + fragment_size)
    #         structure_feature_fragments = np.pad(structure_feature_fragments, ((0, 0), (0, structure_feature_fragments.shape[0])))
    #
    #         # Roll the fragments one by one (add heading zeros)
    #         for i, row in enumerate(structure_feature_fragments):
    #             structure_feature_fragments[i] = np.roll(row, i)
    #
    #         # Calculate average along columns of overlapping fragments (average RMSD per residue)
    #         structure_feature_average = np.average(structure_feature_fragments, axis=0)
    #
    #
    #         #PYMOL SCRIPT
    #
    #         cmd.load("data/pdb{}.pdb".format(j), j)  # Load from file
    #         cmd.remove("resn hoh")  # Remove water molecules
    #         cmd.hide("lines", "all")  # Hide lines
    #         cmd.show("cartoon", j)  # Show cartoon
    #         norm = colors.Normalize(vmin=min(structure_feature_average), vmax=max(structure_feature_average))
    #         for i, residue in enumerate(Selection.unfold_entities(structure[0], "R")):
    #             rgb = cm.bwr(norm(structure_feature_average[i]))
    #             # print(i, residue.id, structure_rmsd_average[i], rgb)
    #             cmd.set_color("col_{}".format(i), list(rgb)[:3])
    #             cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))
    #         cmd.png("data/pymol_image", width=2000, height=2000, ray=1)
    #         h=j