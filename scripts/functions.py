import numpy as np
import math
import glob
import os
from Bio.PDB import PPBuilder
from matplotlib import pyplot as plt
from matplotlib import patches


def extract_filenames(folder, extension):
    names = []
    paths = glob.glob(folder+"/*."+extension)
    for path in paths:
        base = os.path.basename(path)
        file = os.path.splitext(base)[0]
        names.append(file)
    return names


def radius_gyration(chain):

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


def SSRama(structure, dssp_dict):
    # List of Ramachandran areas corresponding to different secondary structure classes
    # E = beta sheet, P = polyproline I && II,
    # H = alpha-helix, R = left-handed helix
    # (lower-left phi and psi, width, height, class, color)
    rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                      (-180, 50, 80, 130, 'E', 'blue'),
                      (-100, -180, 100, 60, 'P', 'green'),
                      (-100, 50, 100, 130, 'P', 'green'),
                      (-180, -120, 180, 170, 'H', 'red'),
                      (0, -180, 180, 360, 'L', 'yellow')]

    # Calculate PSI and PHI
    ppb = PPBuilder()  # PolyPeptideBuilder
    rama = {}  # { chain : [[residue_1, ...], [phi_residue_1, ...], [psi_residue_2, ...] ] }
    # for model in structure:
    for chain in structure:
        for pp in ppb.build_peptides(chain):

            phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]

            for i, residue in enumerate(pp):
                # print(model, chain, i, residue, phi_psi[i])

                # Convert radians to degrees and remove first and last value that are None
                if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                    rama.setdefault(chain.id, [[], [], []])
                    rama[chain.id][0].append(residue)
                    rama[chain.id][1].append(math.degrees(phi_psi[i][0]))
                    rama[chain.id][2].append(math.degrees(phi_psi[i][1]))

    # Get SS from phi/psi and compare with DSSP

    ss = [None]
    for chain_id in rama:
        for residue, phi, psi in zip(*rama[chain_id]):
            ss_class = None
            for x, y, width, height, ss_c, color in rama_ss_ranges:
                if x <= phi < x + width and y <= psi < y + height:
                    ss_class = ss_c
                    break
            ss.append(ss_class)
            # print(residue, ss_class, dssp_dict.get((chain_id, residue.id))[2], phi, psi)
    ss.append(None)

    # # Plot Ramachandran SS regions
    # f, axes = plt.subplots(1, len(rama), figsize=(12, 12))
    # axes = np.array(axes).reshape(
    #     -1)  # Hack to give consistency for single/multiple suplots (-1 force reshape to infer dimensions)
    # for ax, chain_id in zip(axes, rama):
    #
    #     # Plot SS regions
    #     for x, y, width, height, _, color in rama_ss_ranges:
    #         ax.add_patch(patches.Rectangle(
    #             (x, y),  # (x,y)
    #             width,  # width
    #             height,  # height
    #             color=color, zorder=0))  # transparency
    #
    #     # Plot points
    #     ax.scatter(rama[chain_id][1], rama[chain_id][2], s=6, color='black', zorder=1)
    #
    #     ax.set_xlabel('phi')
    #     ax.set_ylabel('psi')
    #
    # plt.tight_layout()  # Remove figure padding
    # plt.show()

    return ss


def get_distance_matrix(residues):

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
