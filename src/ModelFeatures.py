import math
import networkx as nx
import numpy as np
import os
import pandas
import pymol

from Bio.PDB import PDBParser, DSSP, PPBuilder, Superimposer, PDBIO, Selection
from matplotlib import pyplot as plt, patches, colors, cm
from pymol import cgo, cmd, util
from pymol2 import PyMOL
from scipy.spatial.distance import *
from sklearn_extra import cluster
from sklearn.metrics import silhouette_score


def extract_vectors_model_feature(residues, key=None, models=None, features=None, indexes=False, index_slices=False):
    """
    This function allows you to extract information of the model features from the data structure.
    :param residues: number of residues in the model
    :param key: the key of the feature or None
    :param models: the id of model or None
    :param features: matrix of features
    :param indexes: only ruturn begin, end of same feature if it's True, default: False
    :param index_slices: return all the intervals of the features if it's True, default: False
    :return: begin/end, slices or features
    """

    begin = end = -1
    residues = int(residues)
    features = np.array(features)

    slices = []

    if key == 'RG' or index_slices:
        begin = 2
        end = 3
        slices.append(slice(begin, end))

    if key == 'ASA' or index_slices:
        begin = 3
        end = residues + 3
        slices.append(slice(begin, end))

    if key == 'SS' or index_slices:
        begin = residues + 3
        end = 2 * residues + 3
        slices.append(slice(begin, end))

    if key == 'DIST' or index_slices:
        begin = 2 * residues + 3
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

    if models is None:
        return features[:, begin:end]
    else:
        if isinstance(models, int):
            return np.array(features[models][begin:end])
        else:
            return features[models, begin:end]


class ModelFeatures:

    def __init__(self, folder, identifier):

        # Path to the PDB
        self._path = os.path.join(folder, identifier + '.pdb')  # TODO: .ent too?
        # PDB id
        self._id = identifier

        # E = beta sheet, P = polyproline I && II,
        # H = alpha-helix, L = left-handed helix
        self._ranges = [
            (-180, -180, 80, 60, 'E', 'blue'),
            (-180, 50, 80, 130, 'E', 'blue'),
            (-100, -180, 100, 60, 'P', 'green'),
            (-100, 50, 100, 130, 'P', 'green'),
            (-180, -120, 180, 170, 'H', 'red'),
            (0, -180, 180, 360, 'L', 'yellow')
        ]

        self._features = []
        self._residues = None
        self._conformations = []
        # Folder to the feature file
        # self._folder = folder + '/model-features/'  # create model-features inside the input folder
        self._folder = 'data/model-features/'         # create model-features always inside data folder
        os.makedirs(self._folder, exist_ok=True)

        self._output_folder = 'output/model-features'  # create model-features always inside output folder
        os.makedirs(self._output_folder, exist_ok=True)

        # Feature file name
        self._file = identifier + '_features.csv'

        # Clustering
        self._centroids = []
        self._labels = []

        self._max_radius = None

    def choice_maker(self):
        """
        This function allows you to check if the file containing the features of
        the selected ped has already been created. If yes, it loads the features
        to be used later, otherwise it calculates the features and saves the file.

        :return: the features vector
        """

        if os.path.exists(self._folder + self._file):
            print('\n\t- Loading features...')
            self.extract(self._folder + self._file)
        else:
            print('\n\t- Computing features...')
            self.compute()
            self.save(self._folder + self._file)

        self._conformations = len(self._features)

        radius = np.array(extract_vectors_model_feature(self._residues, key='RG', features=self._features))
        self._max_radius = float(max(radius))

    def compute(self):
        """
        This function allows to calculate the features belonging to every single
        model of the analyzed structure.
        """

        structure = PDBParser(QUIET=True).get_structure(self._id, self._path)

        self._residues = int(len(list(structure[0]['A'].get_residues())))

        for model in structure:

            features = {
                'N': self._residues,
                'RG': self.compute_gyration_radius(model['A']),
                'ASA': [],
                'SS': [],
                'DIST': []
            }

            dssp = DSSP(model, self._path, dssp="binx/dssp/mkdssp")  # WARNING Check the path of mkdssp

            for ss in dssp:
                features['ASA'].append(ss[3])

            if len(features['ASA']) < self._residues:
                features['ASA'] = np.concatenate(
                    (features['ASA'], np.zeros(self._residues - len(features['ASA']))))

            features['SS'] = self.compute_secondary_structure(model)

            dist = np.asmatrix(self.compute_distance_matrix(list(model.get_residues())))
            features['DIST'] = np.triu(dist, 1).reshape((1, -1))
            features['DIST'] = features['DIST'][features['DIST'] != 0]

            features_model = [model.id]
            for key in features.keys():

                if key in ['N', 'RG']:
                    features_model.append(features[key])

                if key in ['ASA', 'DIST']:
                    for e in features[key]:
                        features_model.append(e)

                if key in ['SS']:
                    for e in features[key]:
                        if e == 'E':
                            features_model.append(1)
                        if e == 'P':
                            features_model.append(2)
                        if e == 'H':
                            features_model.append(3)
                        if e == 'L':
                            features_model.append(4)
                        if e == '-':
                            features_model.append(0)

            self._features.append(features_model)

    def compute_gyration_radius(self, chain):
        """
        This function calculates the radius of gyration of a specific
        chain passed in input

        :param chain: chain of the model
        :return: radius of gyration of the chain
        """

        # Heavy atoms coordinates
        coord = list()
        for atom in chain.get_atoms():
            if atom.get_name()[0] in ['C', 'O', 'N', 'S']:
                coord.append(atom.get_coord())
        coord = np.array(coord)  # N X 3

        barycenter = np.sum(coord, axis=0) / coord.shape[0]  # center of mass is more correct

        # Calculate distance of each atom from the barycenter
        dist = coord - barycenter
        dist = dist * dist
        dist = np.sqrt(np.sum(dist, axis=1))

        return round(math.sqrt(np.sum(dist * dist) / len(coord)), 3)

    def compute_secondary_structure(self, model):
        """
        This function defines all the secondary structures of the model passed in input
        :param model: one model
        :return: the matrix of secondary structures
        """

        # Calculate PSI and PHI
        ppb = PPBuilder()  # PolyPeptideBuilder
        rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }

        for chain in model:
            for pp in ppb.build_peptides(chain):
                phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
                for i, residue in enumerate(pp):
                    # Convert radians to degrees and remove first and last value that are None
                    if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                        rama.setdefault(chain.id, [[], [], []])
                        rama[chain.id][0].append(residue)
                        rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
                        rama[chain.id][2].append(math.degrees(phi_psi[i][1]))
                    else:
                        rama.setdefault(chain.id, [[], [], []])
                        rama[chain.id][0].append(residue)
                        rama[chain.id][1].append(math.nan)
                        rama[chain.id][2].append(math.nan)

            if len(rama[chain.id][0]) < self._residues:
                for i in range(self._residues - len(rama[chain.id][0])):
                    rama.setdefault(chain.id, [[], [], []])
                    rama[chain.id][0].append(None)
                    rama[chain.id][1].append(math.nan)
                    rama[chain.id][2].append(math.nan)

        # Get SS from phi/psi and compare with DSSP

        ss = []
        for chain_id in rama:
            for residue, phi, psi in zip(*rama[chain_id]):
                ss_class = None
                if math.isnan(phi) and math.isnan(psi):
                    ss_class = '-'
                else:
                    for x, y, width, height, ss_c, color in self._ranges:
                        if x <= phi < x + width and y <= psi < y + height:
                            ss_class = ss_c
                            break

                ss.append(ss_class)

        return ss

    def compute_distance_matrix(self, residues):
        """
        This function computes the distance matrix. It computes the pairwise distances
        between residues
        :param residues: residues of chains
        :return: the distance matrix
        """

        # Calculate the distance matrix
        distances = []
        for residue1 in residues:
            if residue1.id[0] == " " and residue1.has_id("CA"):  # Exclude hetero/water residues
                row = []
                for residue2 in residues:
                    if residue2.id[0] == " " and residue2.has_id("CA"):  # Exclude hetero/water residues
                        row.append(residue1["CA"] - residue2["CA"])
                distances.append(row)

        np.array(distances).reshape(len(residues), -1)

        return distances

    def extract(self, output):
        """
        This function allows you to extract features from the output file
        :param output: the name of output file
        :return: feature
        """

        df = pandas.read_csv(output, index_col=None, header=None)

        self._features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self._features[row, :] = np.array(df.iloc[row])

        self._residues = int(df.iloc[0][1])

        return self._features

    def save(self, output):
        """
        This function saves the models features as a file
        :param output: name of output file
        :return: True or False, If the output file of all the models features is saved correctly
        """

        with open(output, 'w') as f:
            for model in self._features:
                f.write("%d" % model[0])  # index of the model
                for i in range(1, len(model)):
                    f.write(",%f" % model[i])
                f.write("\n")
        print("\t- Successfully saved in {}".format(output))

    # Clustering

    def metrics(self, x, y):
        """
        This function use a specific metric and all the features to compute the distance
        between two input points.
        Radius of gyration (difference), asa (euclidean distance), ss (hamming), distance (correlation)
        :param x: features of one model
        :param y: features of one model
        :return: distance between x and y
        """

        indexes = extract_vectors_model_feature(residues=self._residues, index_slices=True)

        rg = np.abs(x[indexes[0]] - y[indexes[0]]) # /self._max_radius
        asa = euclidean(x[indexes[1]], y[indexes[1]])

        dist_matrix = [[0, 1, 1, 1, 1],
                       [1, 0, 1, 1, 1],
                       [1, 1, 0, 1, 1],
                       [1, 1, 1, 0, 0.5],
                       [1, 1, 1, 0.5, 0]]
        ss = 0
        for i in range(len(x[indexes[2]])):
            ss += dist_matrix[int(x[indexes[2]][i])][int(y[indexes[2]][i])]*2

        ss = ss / len(x[indexes[2]])
        # ss = hamming(x[indexes[2]], y[indexes[2]])
        dist = cosine(x[indexes[3]], y[indexes[3]])

        metric = rg + asa + ss + dist

        return metric

    def compute_clustering(self, k_set=range(3, 9)):
        """
        This function clustering all the models in a ped using KMedoids and a customed metric function.
        It also finds the best number of clusters (representative conformations) by the silhouette score.
        :param k_set: interval of the number of clusters
        :return: the centroids id of each cluster
        """

        print("\t- Clustering...")

        silhouettes = []

        for k in k_set:
            kMed = cluster.KMedoids(n_clusters=k, metric=self.metrics, init='k-medoids++', max_iter=1000)
            pred_labels = kMed.fit_predict(self._features)

            s = silhouette_score(self._features, pred_labels, metric=self.metrics)
            silhouettes.append(s)

        s_max = max(silhouettes)
        k_opt = k_set[silhouettes.index(s_max)]

        print('Number of representative conformations: {}'.format(k_opt))
        print('Correspondent silhouette value: {}'.format(s_max))

        kMed = cluster.KMedoids(n_clusters=k_opt, metric=self.metrics, init='k-medoids++', max_iter=1000)
        self._labels = kMed.fit_predict(self._features)
        self._centroids = kMed.medoid_indices_
        print('Indexes of the representative conformations: {}'.format(self._centroids))

    def generate_graph(self):
        """
        This function generates the weighted graph of the representative conformatios extracted in
        the clustering process
        :return: the graph
        """

        g = nx.Graph()
        for medoid in self._centroids:
            g.add_node(medoid)

        for i in range(len(self._centroids)):
            for j in range(i + 1, len(self._centroids)):
                g.add_edge(self._centroids[i], self._centroids[j],
                           weight=np.sum(self._features[self._centroids[i]] - self._features[self._centroids[j]]))

        options = {'node_color': 'orange', 'node_size': 700, 'width': 2}
        nx.draw_spectral(g, with_labels=True, **options)
        path = "{}/{}_graph.png".format(self._output_folder, self._id)
        plt.savefig(path)
        plt.show()
        return g

    def pymol_metric(self, asa, ss, dist):
        """
        This function computes the distance between matching residues. It is used to generate the pymol image.
        :param asa: relative accessible surface area of conformation
        :param ss: secondary structure of conformation
        :param dist: distance matrix of conformation
        :return: distance between conformations
        """

        asa_dist = np.std(asa)

        sum_ss = 0
        for i in range(len(ss)):
            for j in range(i, len(ss)):
                sum_ss += 0 if ss[i] == ss[j] else 1

        sum_dist = 0.0
        for k in range(dist.shape[0]):
            for h in range(k, dist.shape[0]):
                sum_dist += 1 - correlation(dist[k], dist[h])

        return asa_dist + sum_ss + sum_dist

    def generate_pymol_image(self, g):
        """
        This function generates the pymol image of one ensemble. It scales the structure color with respect
        the distance among representative conformations using own pymol metric.
        :param g: graph of representative conformations of one ensemble
        :return: pymol image
        """

        # pymol.finish_launching()  # Open Pymol
        # p = PyMOL()
        # p.start()

        # Load the structure conformations
        structure = PDBParser(QUIET=True).get_structure(self._id, self._path)

        # Build the matrix with the features of the representative conformations
        representative_features = []
        first = True
        for model in g.nodes():
            if first:
                io = PDBIO()
                io.set_structure(structure[model])
                io.save("data/pymol/{}.pdb".format(self._id))
            representative_features.append(self._features[model])
        representative_features = np.matrix(representative_features)

        # Build feature matrices for representative conformations
        asa = extract_vectors_model_feature(self._residues, key='ASA', features=representative_features)
        ss = extract_vectors_model_feature(self._residues, key='SS', features=representative_features)
        temp_dist = extract_vectors_model_feature(self._residues, key='DIST', features=representative_features)

        # Reshape distance matrix (linear to square)
        dist = np.zeros((len(representative_features), self._residues, self._residues))
        for k in range(len(representative_features)):
            idx = np.triu_indices(self._residues, k=1)
            dist[k][idx] = temp_dist[k]
            dist[k] = dist[k] + dist[k].T

        # For each residue, apply the metric and save it
        residue_variability = []
        for residue in range(self._residues):
            residue_variability.append(self.pymol_metric(asa[:, residue], ss[:, residue], dist[:, residue]))

        cmd.delete("all")
        cmd.load("data/pymol/{}.pdb".format(self._id), self._id)  # Load from file
        cmd.remove("resn hoh")  # Remove water molecules
        cmd.hide("lines", "all")  # Hide lines
        cmd.show("cartoon", "all")  # Show cartoon
        norm = colors.Normalize(vmin=min(residue_variability), vmax=max(residue_variability))
        for i, residue in enumerate(Selection.unfold_entities(structure[0], "R")):
            rgb = cm.bwr(norm(residue_variability[i]))
            cmd.set_color("col_{}".format(i), list(rgb)[:3])
            cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))
        cmd.png("{}/{}_pymol.png".format(self._output_folder, self._id), width=4000, height=2000, ray=1)

        # p.stop()
