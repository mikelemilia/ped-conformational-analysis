import csv
import math

import numpy as np
import pandas
from Bio.PDB import PDBParser, DSSP, PPBuilder
from matplotlib import pyplot as plt, patches


class Features:

    def __init__(self, path, identifier):
        self.path = path
        self.id = identifier

        # E = beta sheet, P = polyproline I && II,
        # H = alpha-helix, L = left-handed helix
        self.ranges = [
            (-180, -180, 80, 60, 'E', 'blue'),
            (-180, 50, 80, 130, 'E', 'blue'),
            (-100, -180, 100, 60, 'P', 'green'),
            (-100, 50, 100, 130, 'P', 'green'),
            (-180, -120, 180, 170, 'H', 'red'),
            (0, -180, 180, 360, 'L', 'yellow')
        ]

        self.features = []
        self.residues = None

        self.compute = False

    def compute_features(self):

        structure = PDBParser(QUIET=True).get_structure(self.id, self.path)

        self.residues = int(len(list(structure[0]['A'].get_residues())))

        for model in structure:

            features = {
                'N': self.residues,
                'RG': self.compute_gyration_radius(model['A']),
                'ASA': [],
                'SS': [],
                'DIST': []
            }

            dssp = DSSP(model, self.path, dssp="binx/dssp/mkdssp")  # WARNING Check the path of mkdssp

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

            self.features.append(features_model)

        return self.features

    def compute_gyration_radius(self, chain):

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

        # Get SS from phi/psi and compare with DSSP

        ss = []
        for chain_id in rama:
            for residue, phi, psi in zip(*rama[chain_id]):
                ss_class = None
                for x, y, width, height, ss_c, color in self.ranges:
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

    def extract_features(self, output):

        df = pandas.read_csv(output, index_col=None, header=None)

        self.features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self.features[row, :] = np.array(df.iloc[row])

        self.residues = int(df.iloc[0][1])

        return self.features

        # row = csv
        # features = {'N': 0, 'RG': 0, 'ASA': [], 'SS': [], 'DIST': []}
        #
        # model_id = row[0]
        # features['N'] = int(row[1])
        # features['RG'] = float(row[2])
        #
        # for elem in range(3, features['N'] + 3):
        #     features['ASA'].append(row[3 + elem])
        #
        # for elem in range(features['N'] + 3, 2 * features['N'] - 2 + 3):  # removed first and last element (None)
        #     features['SS'].append(row[elem])
        #
        # for elem in range(2 * features['N'] - 2 + 3, len(csv)):
        #     features['DIST'].append(row[elem])
        #
        # return model_id, features

    def save(self, output):
        with open(output, 'w') as f:
            for model in self.features:
                f.write("%s" % model[0])  # index of the model
                for i in range(1, len(model)):
                    f.write(",%f" % model[i])
                f.write("\n")
        print("{}_features.csv saved".format(self.id))

    def plot_secondary_structure(self, rama):
        # Plot Ramachandran SS regions
        f, axes = plt.subplots(1, len(rama), figsize=(12, 12))
        axes = np.array(axes).reshape(
            -1)  # Hack to give consistency for single/multiple suplots (-1 force reshape to infer dimensions)
        for ax, chain_id in zip(axes, rama):

            # Plot SS regions
            for x, y, width, height, _, color in self.ranges:
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

