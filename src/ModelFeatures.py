import math
import networkx as nx
import numpy as np
import os
import pandas

from Bio.PDB import PDBParser, DSSP, PPBuilder, PDBIO, Selection
from matplotlib import pyplot as plt, colors, cm
from pymol import cmd
from scipy.spatial.distance import *
from sklearn_extra import cluster
from sklearn.metrics import silhouette_score


def extract_vectors_model_feature(residues, key=None, models=None, features=None, indexes=False, index_slices=False):
    """
    This function allows you to extract information of the model features from the data structure. In particular allows:
    - all rows or a specific subset of them, containing a certain feature (i.e., RG, ASA, SS, etc ...)
    - the interval extremes for a certain features (i.e., RG, ASA, SS, etc ...)
    - all the feature intervals as slices
    :param residues: number of residues in the model
    :param key: the key of the feature or None if considering all of them, default: False
    :param models: the models id or None if considering all of them, default: False
    :param features: matrix of features or None if extracting only the indexes, default: False
    :param indexes: return (begin, end) indexes of a feature if it's True, default: False
    :param index_slices: return all the intervals of the features if it's True, default: False
    :return: begin/end, slices or features
    """

    residues = int(residues)
    features = np.array(features)

    begin = end = -1
    slices = []

    # For each possible key (or if index_slices is required), extract the column of the features in the features matrix
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

    # The requested one has not been found
    if begin == -1:
        return None

    # If all the indexes have been requested, return the whole slice
    if index_slices:
        return slices

    # If the key indexes are requested or if the features matrix is not provided, return the pair (begin, end)
    if indexes is True or features is None:
        return begin, end

    # Return the features requested contained in features and only for the models requested (if not specified, all)
    if models is None:
        return features[:, begin:end]
    else:
        if isinstance(models, int):
            return np.array(features[models][begin:end])
        else:
            return features[models, begin:end]


class ModelFeatures:
    """
    Class for the computation of the models features within a single ensemble (of ID PEDxxxxxexxx)
    """

    def __init__(self, folder, identifier):

        self._path = os.path.join(folder, identifier + '.pdb')      # Path to the PDB
        self._id = identifier                           # PDB ID (PEDxxxxxexxx)
        self._file = identifier + '_features.csv'       # Feature file name
        self._folder = 'data/model-features/'           # Create model-features always inside data folder
        self._output_folder = 'output/model-features'   # Create model-features always inside output folder
        self._file_path = self._folder + self._file     # File path where to store the features

        # Create folders if not existing
        os.makedirs(self._folder, exist_ok=True)
        os.makedirs(self._output_folder, exist_ok=True)

        # Initialization for the subsequent analysis
        self._features = []         # matrix of features
        self._residues = 0          # number of residues
        self._conformations = 0     # number of conformations
        self._centroids = []        # clustering centroids
        self._labels = []           # label associated to each conformation with clustering

        # Regions of Ramachandran plot
        self._ranges = [
            (-180, -180, 80, 60, 'E', 'blue'),      # E = beta sheet
            (-180, 50, 80, 130, 'E', 'blue'),
            (-100, -180, 100, 60, 'P', 'green'),    # P = polyproline I && II
            (-100, 50, 100, 130, 'P', 'green'),
            (-180, -120, 180, 170, 'H', 'red'),     # H = alpha-helix
            (0, -180, 180, 360, 'L', 'yellow')      # L = left-handed helix
        ]
        # Distance matrix for the secondary structure values
        self._dist_matrix = [[0, 1, 1, 1, 1],
                             [1, 0, 1, 1, 1],
                             [1, 1, 0, 1, 1],
                             [1, 1, 1, 0, 0.5],
                             [1, 1, 1, 0.5, 0]]

    def choice_maker(self, model_name=''):
        """
        This function allows to check if the file containing the features for the selected PED has already been
        generated and stored. If so, it loads these features, otherwise it calculates the correspondent matrix
        and saves the file.
        :param model_name: model name, eventually to be printed
        :return: the features vector
        """

        if os.path.exists(self._file_path):

            # If file existing, load it
            print('\t- Loading {}model features...'.format(model_name))
            self.extract(self._file_path)

        else:

            # If file not existing, compute all the features and save it
            print('\t- Computing {}model features...'.format(model_name))
            self.compute()
            self.save(self._file_path)

        # Extraction of the number of residues (second column of each row - Note that it is always the same value)
        self._residues = self._features[0][1]
        # Extraction the number of conformations present inside this PED (number of rows)
        self._conformations = len(self._features)

        return self._features

    def compute(self):
        """
        This function allows to calculate the features belonging to every model of the analyzed structure, building the
        correspondent dictionary.
        """

        # Extract the structure and save the number of residues
        structure = PDBParser(QUIET=True).get_structure(self._id, self._path)
        self._residues = int(len(list(structure[0]['A'].get_residues())))

        # Loop to analyze every conformation inside PED
        for model in structure:

            # Dictionary of the correspondent features
            features = {
                'N': self._residues,    # Number of residues (needed for computation)
                'RG': self.compute_gyration_radius(model['A']),     # Radius of gyration
                'ASA': [],
                'SS': [],
                'DIST': []
            }

            # Relative accessible surface area calculated thanks to DSSP
            dssp = DSSP(model, self._path, dssp="binx/dssp/mkdssp")  # WARNING Check the path of mkdssp

            for ss in dssp:
                features['ASA'].append(ss[3])

            if len(features['ASA']) < self._residues:
                features['ASA'] = np.concatenate(
                    (features['ASA'], np.zeros(self._residues - len(features['ASA']))))

            # Secondary structure
            features['SS'] = self.compute_secondary_structure(model)

            # Distance matrix: only the upper triangular (without diagonal) is maintained, since it is the
            # meaningful and not redundant information brought
            dist = np.asmatrix(self.compute_distance_matrix(list(model.get_residues())))
            features['DIST'] = np.triu(dist, 1).reshape((1, -1))
            features['DIST'] = features['DIST'][features['DIST'] != 0]

            # Conversion of the build dictionary in a row 'features_model', then appended to the whole features matrix
            features_model = [model.id]
            for key in features.keys():

                if key in ['N', 'RG']:  # Single value just appended
                    features_model.append(features[key])

                if key in ['ASA', 'DIST']:  # Scanned array to be appended
                    for e in features[key]:
                        features_model.append(e)

                if key in ['SS']:   # Conversion of the secondary structure values to integers and appending
                    for e in features[key]:
                        if e == 'E':
                            features_model.append(1)
                        if e == 'P':
                            features_model.append(2)
                        if e == 'H':
                            features_model.append(3)
                        if e == 'L':
                            features_model.append(4)
                        if e == '-':    # Value added if the Ramachandran region is not available
                            features_model.append(0)

            self._features.append(features_model)

    def compute_gyration_radius(self, chain):
        """
        This function calculates the radius of gyration of a specific chain passed in input based on its definition
        :param chain: chain of the model
        :return: radius of gyration of the chain
        """

        # Heavy atoms coordinates and barycenter calculation
        coord = list()
        for atom in chain.get_atoms():
            if atom.get_name()[0] in ['C', 'O', 'N', 'S']:
                coord.append(atom.get_coord())
        coord = np.array(coord)  # N X 3

        barycenter = np.sum(coord, axis=0) / coord.shape[0]

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
        ppb = PPBuilder()
        rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }

        for chain in model:
            for pp in ppb.build_peptides(chain):
                phi_psi = pp.get_phi_psi_list()

                for i, residue in enumerate(pp):
                    if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                        # Conversion to degrees when the values are not None (for first and last)
                        rama.setdefault(chain.id, [[], [], []])
                        rama[chain.id][0].append(residue)
                        rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
                        rama[chain.id][2].append(math.degrees(phi_psi[i][1]))
                    else:
                        # Adding of Nan if the angles are None (for first and last)
                        rama.setdefault(chain.id, [[], [], []])
                        rama[chain.id][0].append(residue)
                        rama[chain.id][1].append(math.nan)
                        rama[chain.id][2].append(math.nan)

            # Eventual nan-padding if something goes wrong during the angle computation
            if len(rama[chain.id][0]) < self._residues:
                for i in range(self._residues - len(rama[chain.id][0])):
                    rama.setdefault(chain.id, [[], [], []])
                    rama[chain.id][0].append(None)
                    rama[chain.id][1].append(math.nan)
                    rama[chain.id][2].append(math.nan)

        # Comparison of the angles with the Ramachandran regions
        ss = []
        for chain_id in rama:
            for residue, phi, psi in zip(*rama[chain_id]):
                ss_class = None
                if math.isnan(phi) and math.isnan(psi):
                    # If nan (angles not available) insert a symbol indicating this situation
                    ss_class = '-'
                else:
                    # Determine the correspondent region and store it
                    for x, y, width, height, ss_c, color in self._ranges:
                        if x <= phi < x + width and y <= psi < y + height:
                            ss_class = ss_c
                            break

                ss.append(ss_class)

        return ss

    def compute_distance_matrix(self, residues):
        """
        This function computes the distance matrix, computing the pairwise distances between residues
        :param residues: residues of chains
        :return: the distance matrix
        """

        distances = []
        for residue1 in residues:
            if residue1.id[0] == " " and residue1.has_id("CA"):  # Exclude hetero/water residues and atoms without CA
                row = []
                for residue2 in residues:
                    if residue2.id[0] == " " and residue2.has_id("CA"):  # Exclude hetero/water residues
                        row.append(residue1["CA"] - residue2["CA"])
                distances.append(row)

        # Reshape based on the number of residues
        np.array(distances).reshape(len(residues), -1)

        return distances

    def extract(self, path):
        """
        This function allows to extract features from the saved features file
        :param path: the path of features file
        """

        # Read of the file as dataframe and initialization of the matrix
        df = pandas.read_csv(path, index_col=None, header=None)
        self._features = np.full((df.shape[0], df.shape[1]), None)

        # Scan of the dataframe rows and save in the data matrix
        for row in range(0, df.shape[0]):
            self._features[row, :] = np.array(df.iloc[row])

    def save(self, output):
        """
        This function saves the models features matrix as a csv file
        :param output: name of output file
        """

        with open(output, 'w') as f:
            # Scan the row (models) and write each value (separated with commas) as float
            for model in self._features:
                f.write("%d" % model[0])  # Index of the model
                for i in range(1, len(model)):
                    f.write(",%f" % model[i])
                f.write("\n")
        print("\t- Successfully saved in {}".format(output))

    # ----------------------------------------------------------------------------------------
    # Model features analysis

    def metrics(self, x, y):
        """
        This function implements a specific metric for all the extracted features to compute the distance
        between two model features vectors.
            rg (absolute difference), asa (euclidean distance),
            ss (normalized hamming with scoring matrix), dist (cosine distance)
        :param x: features vector of one model
        :param y: features vector of one model
        :return: distance between x and y
        """

        # Extract the index slices for each feature contained in the vector - to subsequently divide x and y in their
        # corresponding features vectors
        indexes = extract_vectors_model_feature(residues=self._residues, index_slices=True)

        # Radius of gyration (absolute difference)
        rg = np.abs(x[indexes[0]] - y[indexes[0]])
        # Relative accessible surface area (euclidean distance)
        asa = euclidean(x[indexes[1]], y[indexes[1]])

        # Secondary Structure (hamming with scoring matrix)
        ss = 0
        for i in range(self._residues):
            ss += self._dist_matrix[int(x[indexes[2]][i])][int(y[indexes[2]][i])]
        ss = ss / self._residues

        # Distance (cosine)
        dist = cosine(x[indexes[3]], y[indexes[3]])

        metric = rg + asa + ss + dist

        return metric[0]

    def compute_clustering(self, k_set=range(3, 9)):
        """
        This function clusters PED models using KMedoids approach and a customized metric function.
        It also determins the best number of clusters (representative conformations) exploiting silhouette score.
        :param k_set: interval of the possible number of clusters
        """

        print("\t- Clustering...")

        silhouettes = []

        # Scan the possible number of clusters, determine the corresponding clustering and silhouette
        for k in k_set:
            kMed = cluster.KMedoids(n_clusters=k, metric=self.metrics, init='k-medoids++', max_iter=1000)
            pred_labels = kMed.fit_predict(self._features)

            s = silhouette_score(self._features, pred_labels, metric=self.metrics)
            silhouettes.append(s)

        # Determine the optimal number of clusters as the one that maximize the silhouette value
        s_max = max(silhouettes)
        k_opt = k_set[silhouettes.index(s_max)]

        print('\n\t\tNumber of representative conformations: {}'.format(k_opt))
        print('\t\tCorrespondent silhouette value: {}'.format(s_max))

        # Clustering recomputation for the optimal number of clusters
        kMed = cluster.KMedoids(n_clusters=k_opt, metric=self.metrics, init='k-medoids++', max_iter=1000)
        self._labels = kMed.fit_predict(self._features)
        self._centroids = kMed.medoid_indices_
        print('\t\tIndexes of the representative conformations: {}'.format(self._centroids))

    def generate_graph(self):
        """
        This function generates the weighted graph of the representative conformations extracted in
        the clustering process
        :return: the graph
        """

        g = nx.Graph()

        # Add the representative conformations (centroids) as nodes
        for medoid in self._centroids:
            g.add_node(medoid)

        # Add the edges between the nodes weighted with the distance between the correspondent conformations
        labels = {}
        for i in range(len(self._centroids)):
            for j in range(i + 1, len(self._centroids)):
                dist = self.metrics(self._features[self._centroids[i]], self._features[self._centroids[j]])
                g.add_edge(self._centroids[i], self._centroids[j], weight=dist)
                labels.update({(self._centroids[i], self._centroids[j]): round(dist, ndigits=4)})

        # Print a customized graph and save it
        options = {'node_color': 'orange', 'node_size': 700, 'width': 2}
        pos = nx.spring_layout(g, weight='weight')
        nx.draw(g, with_labels=True, pos=pos, **options)
        nx.draw_networkx_edge_labels(g, pos, edge_labels=labels, font_color='red')

        path = "{}/{}_graph.png".format(self._output_folder, self._id)
        plt.savefig(path)
        plt.show()

        return g

    def pymol_metric(self, asa, ss, dist):
        """
        This function computes the variability of a residue based on the features of a window of residues around it
        (already computed). It is used to generate the pymol image.
            asa (mean of the standard deviations considering separately each residue),
            ss (normalized hamming with scoring matrix considering separately each residue),
            dist (cosine distance considering separately each residue)
        :param asa: relative accessible surface area of the conformations for the residues window under analysis
        :param ss: secondary structure of the conformations for the residues window under analysis
        :param dist: distance matrix of the conformations for the residues window under analysis
        :return: variability of the central residue
        """

        # Relative accessible surface area (mean of the standard deviations)
        asa = np.array(asa, dtype='float64')
        asa_dist = np.mean(np.std(asa, axis=0))

        # Secondary structures (normalized hamming distance)
        sum_ss = 0
        for residue in range(ss.shape[1]):
            for i in range(len(ss)):
                for j in range(i, len(ss)):
                    sum_ss += self._dist_matrix[int(ss[i, residue])][int(ss[j, residue])]
        sum_ss = sum_ss / ss.shape[1]

        # Distance matrix (cosine distance)
        sum_dist = 0.0
        for residue in range(ss.shape[1]):
            for k in range(dist.shape[0]):
                for h in range(k, dist.shape[0]):
                    sum_dist += cosine(dist[k, residue], dist[h, residue])

        return asa_dist + sum_ss + sum_dist

    def generate_pymol_image(self, g):
        """
        This function generates the pymol image representing one ensemble. It scales the structure color with respect
        the variability of the residues among representative conformations using own pymol metric.
        :param g: graph of representative conformations of one ensemble
        """

        # Load the structure conformations
        structure = PDBParser(QUIET=True).get_structure(self._id, self._path)

        # Build the matrix with the features of the representative conformations
        representative_features = []
        first = True

        for model in g.nodes():
            representative_features.append(self._features[model])
            # If not done yet, open the structure for Pymol
            if first:
                io = PDBIO()
                io.set_structure(structure[model])
                io.save("data/pymol/{}.pdb".format(self._id))

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

        # For each residue, determine the window around it and apply the metric
        residue_variability = []
        window_size = 9
        for residue in range(self._residues):
            start = max(0, residue - window_size)
            end = min(self._residues - 1, residue + window_size)
            residue_variability.append(self.pymol_metric(asa[:, start:end], ss[:, start:end], dist[:, start:end]))

        # Pymol image generation
        cmd.delete("all")
        cmd.load("data/pymol/{}.pdb".format(self._id), self._id)  # Load from file
        cmd.remove("resn hoh")  # Remove water molecules
        cmd.hide("lines", "all")  # Hide lines
        cmd.show("cartoon", "all")  # Show cartoon
        norm = colors.Normalize(vmin=min(residue_variability), vmax=max(residue_variability))

        # Changing of the residue color according to its window variability
        for i, residue in enumerate(Selection.unfold_entities(structure[0], "R")):
            rgb = cm.bwr(norm(residue_variability[i]))
            cmd.set_color("col_{}".format(i), list(rgb)[:3])
            cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))

        # Pymol image saving
        cmd.png("{}/{}_pymol.png".format(self._output_folder, self._id), width=4000, height=2000, ray=1)
