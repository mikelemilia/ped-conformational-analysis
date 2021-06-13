import math
import os

import networkx as nx
import numpy as np
import pandas
from Bio.PDB import PDBParser, DSSP, PPBuilder
from matplotlib import pyplot as plt, patches
from scipy.spatial.distance import *
import pymol
from pymol import cgo, cmd, util
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

    def generate_pymol_img(self, g):
        pymol.finish_launching()  # Open Pymol

        # Input 1jsu (3 chains, non-standard amino acids), 1az5 (disordered loops, chain breaks)
        pdb_id = '1jsu'

        cmd.fetch(pdb_id, pdb_id, path="data/")  # Download the PDB
        # cmd.load("data/pdb{}.ent".format(pdb_id), pdb_id)  # Load from file

        cmd.remove("resn hoh")  # Remove water molecules
        cmd.hide("lines", "all")  # Hide lines
        cmd.show("cartoon", pdb_id)  # Show cartoon
        cmd.show("sticks", "hetatm")  # Show hetero atoms as sticks
        # # cmd.spectrum(selection="all")  # Rainbow color
        util.cbc(selection="all")  # Color by chain
        util.cnc(selection="all")  # Color by atom type, but not the C atoms

        # Select and color two residues
        sele_name = "nodes"
        res1 = 'B/200/'
        res2 = 'C/52/'
        cmd.select(sele_name, '{} or {}'.format(res1, res2))
        cmd.show("spheres", sele_name)
        cmd.set('sphere_transparency', 0.5, sele_name)  # Set transparency

        # Get coordinates of two atoms
        atom1 = 'B/200/SD'
        atom2 = 'C/52/CE'
        coord1 = cmd.get_coords(atom1, 1)  # shape N * 3, where N is the number of atoms in the selection
        coord2 = cmd.get_coords(atom2, 1)
        # model = cmd.get_model(pdb_id, 1)  # Fastest way to get all atom coordinates

        # Calculate center of mass between two residues and create a new "pseudo" atom
        center_of_mass = (coord1[0] + coord2[0]) / 2
        print(coord1, coord2, center_of_mass)

        obj_name = "ps_atom"
        cmd.pseudoatom(obj_name, pos=list(center_of_mass))
        cmd.show("spheres", obj_name)
        # cmd.extract(...  # Move selected atoms to a new object
        # cmd.create(...  # Create a new object from selection

        cr, cg, cb = (1.0, 0.0, 0.0)  # RGB red

        # Create lines object
        obj = [cgo.BEGIN, cgo.LINES, cgo.COLOR, cr, cg, cb]
        obj.append(cgo.VERTEX)
        obj.extend(list(coord1[0]))
        obj.append(cgo.VERTEX)
        obj.extend(list(coord2[0]))
        obj.append(cgo.END)

        # Set the object
        obj_name = 'edges'
        cmd.load_cgo(obj, obj_name)
        cmd.set("cgo_line_width", float(3), obj_name)

        cmd.orient(pdb_id)  # Set the origin to full protein
       # cmd.extend("optAlign", ModelFeatures)

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
