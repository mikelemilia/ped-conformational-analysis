import csv
import math
import os
import sys

import numpy as np
import pandas
from Bio.PDB import PDBParser, DSSP, PPBuilder
from matplotlib import pyplot as plt, patches


def extract_vectors_feature(features, key, models=slice(None)):
    residues = int(features[0, 1])
    if key == 'N':
        return residues

    if key == 'RG':
        return features[models, 2]

    if key == 'ASA':
        return features[models, 3:(residues + 3)]

    if key == 'SS':
        return features[models, (residues + 3):(2 * residues + 3)]  # removed first and last

    if key == 'DIST':
        return features[models, (2 * residues + 3):]

    return None


class Model_Features:

    def __init__(self, folder, identifier):
        # Path to the PDB
        self.path = os.path.join(folder, identifier + '.pdb')
        # PDB id
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
        # Folder to the feature file
        self.folder = folder + '/model_features/'
        os.makedirs(self.folder, exist_ok=True)
        # Feature file name
        self.file = identifier + '_features.csv'

    def choice_maker(self):
        if os.path.exists(self.folder + self.file):
            print('\nLoading features...')
            feat = self.extract(self.folder + self.file)
        else:
            print('\nComputing features...')
            feat = self.compute()
            self.save(self.folder + self.file)
        return feat

    def compute(self):

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
                        if e == '-':
                            features_model.append(0)

            print(self.id, len(features['ASA']),  len(features['SS']),  len(features['DIST']))

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
                    else:
                        rama.setdefault(chain.id, [[], [], []])
                        rama[chain.id][0].append(residue)
                        rama[chain.id][1].append(math.nan)
                        rama[chain.id][2].append(math.nan)

            if len(rama[chain.id][0]) < (self.residues):
                for i in range(self.residues - len(rama[chain.id][0])):
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

    def extract(self, output):

        df = pandas.read_csv(output, index_col=None, header=None)

        self.features = np.full((df.shape[0], df.shape[1]), None)
        for row in range(0, df.shape[0]):
            self.features[row, :] = np.array(df.iloc[row])

        self.residues = int(df.iloc[0][1])

        return self.features

    def save(self, output):

        with open(output, 'w') as f:
            for model in self.features:
                f.write("%d" % model[0])  # index of the model
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
