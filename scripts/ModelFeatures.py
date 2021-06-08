import math
import os

import networkx as nx
import numpy as np
import pandas
from Bio.PDB import PDBParser, DSSP, PPBuilder
from matplotlib import pyplot as plt, patches
from scipy.spatial.distance import *
from sklearn.metrics import silhouette_score
from sklearn_extra import cluster


def extract_vectors_model_feature(residues, key=None, models=None, features=None, indexes=False, index_slices=False):
    begin = end = -1
    residues = int(residues)

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
        self._path = os.path.join(folder, identifier + '.pdb')
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
        self._folder = folder + '/model_features/'
        os.makedirs(self._folder, exist_ok=True)
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
            print('\nLoading features...')
            self.extract(self._folder + self._file)
        else:
            print('\nComputing features...')
            self.compute()
            self.save(self._folder + self._file)

        self._conformations = len(self._features)

        radius = np.array(extract_vectors_model_feature(self._residues, key='RG', features=self._features))
        self._max_radius = float(max(radius))

        return self._features

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

            # print(self._id, len(features['ASA']), len(features['SS']), len(features['DIST']))

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
        # print(dist)

        return round(math.sqrt(np.sum(dist * dist) / len(coord)), 3)

    def compute_secondary_structure(self, model):

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

        df = pandas.read_csv(output, index_col=None, header=None)

        self._features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self._features[row, :] = np.array(df.iloc[row])

        self._residues = int(df.iloc[0][1])

        return self._features

    def save(self, output):

        with open(output, 'w') as f:
            for model in self._features:
                f.write("%d" % model[0])  # index of the model
                for i in range(1, len(model)):
                    f.write(",%f" % model[i])
                f.write("\n")
        print("{}_features.csv saved".format(self._id))

    # Clustering

    def metrics(self, x, y):

        indexes = extract_vectors_model_feature(residues=self._residues, index_slices=True)

        rg = np.abs(x[indexes[0]] - y[indexes[0]]) #/ self._max_radius
        asa = euclidean(x[indexes[1]], y[indexes[1]])
        ss = hamming(x[indexes[2]], y[indexes[2]])
        dist = 1 - correlation(x[indexes[3]], y[indexes[3]])

        metric = rg + asa + ss + dist

        # print('Metric: {}'.format(metric))

        return metric

    def compute_clustering(self, k_set=range(3, 9)):

        silhouettes = []

        for k in k_set:
            kMed = cluster.KMedoids(n_clusters=k, metric=self.metrics, init='k-medoids++', max_iter=1000)
            pred_labels = kMed.fit_predict(self._features)

            s = silhouette_score(self._features, pred_labels, metric=self.metrics)
            silhouettes.append(s)

            # u_labels = np.unique(labels)
            # for i in u_labels:
            #     plt.scatter(points[labels == i, 0], points[labels == i, 1], label=i)
            # plt.legend()
            # plt.show()

        s_max = max(silhouettes)
        k_opt = k_set[silhouettes.index(s_max)]

        print('Number of representative conformations: {}'.format(k_opt))
        print('Correspondent silhouette value: {}'.format(s_max))

        kMed = cluster.KMedoids(n_clusters=k_opt, metric=self.metrics, init='k-medoids++', max_iter=1000)
        self._labels = kMed.fit_predict(self._features)
        self._centroids = kMed.medoid_indices_
        print('Indexes of the representative conformations: {}'.format(self._centroids))

    def generate_graph(self):

        # Graph
        g = nx.Graph()
        for medoid in self._centroids:
            g.add_node(medoid)

        for i in range(len(self._centroids)):
            for j in range(i + 1, len(self._centroids)):
                g.add_edge(self._centroids[i], self._centroids[j],
                           weight=np.sum(self._features[self._centroids[i]] - self._features[self._centroids[j]]))

        options = {'node_color': 'orange', 'node_size': 700, 'width': 2}
        nx.draw_spectral(g, with_labels=True, **options)
        path = "output/{}_graph.png".format(self._id)
        plt.savefig(path)
        plt.show()
        return g

    # def generate_pymol_img(self, g):
    #     pymol.finish_launching()  # Open Pymol
    #     # Load the structure conformations
    #     h = g.nodes().value()
    #     print(h)
    #     for j in g.nodes().value():
    #         structure = PDBParser(QUIET=True).get_structure(self._id)
    #         # Superimpose all models to the first model, fragment-by-fragment (sliding window)
    #         super_imposer = Superimposer()
    #
    #         structure_feature_fragments = []
    #         window_size = 9
    #         ref_model = [atom for atom in structure[h].get_atoms() if atom.get_name() == "CA"]  # CA of the first model
    #         ref_features = self._features[h]
    #         model_features = []
    #         alt_model = [atom for atom in structure[j].get_atoms() if atom.get_name() == "CA"]  # coords of the model
    #         alt_features = self._features[j]
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

    def plot_secondary_structure(self, rama):

        # Plot Ramachandran SS regions
        f, axes = plt.subplots(1, len(rama), figsize=(12, 12))
        axes = np.array(axes).reshape(
            -1)  # Hack to give consistency for single/multiple suplots (-1 force reshape to infer dimensions)
        for ax, chain_id in zip(axes, rama):

            # Plot SS regions
            for x, y, width, height, _, color in self._ranges:
                ax.add_patch(patches.Rectangle(
                    (x, y),  # (x,y)
                    width,  # width
                    height,  # height
                    color=color, zorder=0))  # transparency

            # Plot points
            ax.scatter(rama[chain_id][1], rama[chain_id][2], s=6, color='black', zorder=1)

            ax.set_xlabel('phi')
            ax.set_ylabel('psi')

        plt.tight_layout()  # Remove figure padding
        plt.show()

    @property
    def residues(self):
        return self._residues

    @property
    def conformations(self):
        return self._conformations
