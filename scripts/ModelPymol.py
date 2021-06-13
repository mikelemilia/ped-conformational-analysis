import os
from ModelFeatures import ModelFeatures
import numpy as np
from Bio.PDB import PDBList, Superimposer, Selection
from Bio.PDB.PDBParser import PDBParser
import pymol
from pymol import cmd
from matplotlib import cm, colors


out_file = "pymol_{}".format(ModelFeatures._id)
window_size = 9

structure = PDBParser(QUIET=True).get_structure(ModelFeatures._id, ModelFeatures._path)
# Fetch a PDB file to the current dir
pdbl = PDBList()
pdbl.retrieve_pdb_file(ModelFeatures._path, pdir='{}', file_format='pdb') # Will save to pdbXXXX.ent


# Superimpose all models to the first model, fragment-by-fragment (sliding window)
super_imposer = Superimposer()
structure_rmsd_fragments = []  # RMSD, no_models X no_fragments X fragment_size
ref_model = [atom for atom in structure[0].get_atoms() if atom.get_name() == "CA"]  # CA of the first model

# Iterate all models
for i, model in enumerate(structure):
    if i > 0:
        model_rmsd = []  # RMSD, no_fragment X fragment_size
        alt_model = [atom for atom in model.get_atoms() if atom.get_name() == "CA"]  # coords of the model

        # Iterate fragments
        for start in range(len(ref_model) - window_size):
            end = start + window_size
            ref_fragment = ref_model[start:end]
            alt_fragment = alt_model[start:end]

            # Calculate rotation/translation matrices
            super_imposer.set_atoms(ref_fragment, alt_fragment)
            # print(super_imposer.rms, super_imposer.rotran)

            # Rotate-translate coordinates
            alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
            alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
            alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]

            # Calculate RMSD
            # https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
            ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
            dist = ref_fragment_coord - alt_fragment_coord
            # rmsd_fragment = np.sqrt(np.sum(dist * dist) / window_size)  # Total RMSD of the fragment. Identical to super_imposer.rms
            rmsd_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment

            model_rmsd.append(rmsd_res)

        structure_rmsd_fragments.append(model_rmsd)


# Calculate average RMSD per position
structure_rmsd_fragments = np.array(structure_rmsd_fragments)  # no_models X no_fragments X fragment_size
# Calculate the RMSD average for each fragments along all models
structure_rmsd_fragments = np.average(structure_rmsd_fragments, axis=0)  # no_fragments X fragment_size
# Pad with right zeros to reach the sequence length (no_fragments + fragment_size)
structure_rmsd_fragments = np.pad(structure_rmsd_fragments, ((0, 0), (0, structure_rmsd_fragments.shape[0])))
print(structure_rmsd_fragments.shape, len(ref_model))

# Roll the fragments one by one (add heading zeros)
for i, row in enumerate(structure_rmsd_fragments):
    structure_rmsd_fragments[i] = np.roll(row, i)


# Calculate average along columns of overlapping fragments (average RMSD per residue)
structure_rmsd_average = np.average(structure_rmsd_fragments, axis=0)
print(structure_rmsd_average.shape)


#################### PYMOL #####################

pymol.finish_launching()  # Open Pymol

#cmd.load(pdb_file, pdb_id)  # Download the PDB
cmd.load("{}".format(ModelFeatures._path), ModelFeatures._path)  # Load from file

cmd.remove("resn hoh")  # Remove water molecules
cmd.hide("lines", "all")  # Hide lines


norm = colors.Normalize(vmin=min(structure_rmsd_average), vmax=max(structure_rmsd_average))
for i, residue in enumerate(Selection.unfold_entities(structure[0], "R")):
    rgb = cm.bwr(norm(structure_rmsd_average[i]))
    # print(i, residue.id, structure_rmsd_average[i], rgb)
    cmd.set_color("col_{}".format(i), list(rgb)[:3])
    cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))

cmd.png(out_file, width=2000, height=2000, ray=1)

###################################################################################
# pymol.finish_launching()  # Open Pymol
# # Load the structure conformations
# h = g.nodes().value().pop()
# for j in g.nodes().value():
#     structure = PDBParser(QUIET=True).get_structure(self._id)
#     # Superimpose all models to the first model, fragment-by-fragment (sliding window)
#     super_imposer = Superimposer()
#
#     structure_feature_fragments = []
#     window_size = 9
#     ref_model = [atom for atom in structure[h].get_atoms() if atom.get_name() == "CA"]  # CA of the first model
#     ref_features = self._features[h]
#     model_features = []
#     alt_model = [atom for atom in structure[j].get_atoms() if atom.get_name() == "CA"]  # coords of the model
#     alt_features = self._features[j]
#
#     # Iterate fragments
#     for start in range(len(ref_model) - window_size):
#         end = start + window_size
#         ref_fragment = ref_model[start:end]
#         alt_fragment = alt_model[start:end]
#
#         # Calculate rotation/translation matrices
#         super_imposer.set_atoms(ref_fragment, alt_fragment)
#
#         # Rotate-translate coordinates
#         # alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
#         # alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
#         # alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]
#
#         #features of structure
#         #ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
#         # dist = np.diff(ref_features,alt_features)
#         # feat_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment
#         # model_features.append(feat_res)
#     dist = np.diff(ref_features,alt_features)
#     feat_res = np.sqrt(np.sum(dist * dist, axis=1))
#     structure_feature_fragments.append(feat_res)
#
#
#     structure_feature_fragments = np.array(structure_feature_fragments)  # no_models X no_fragments X fragment_size
#
#     structure_feature_fragments = np.average(structure_feature_fragments, axis=0)  # no_fragments X fragment_size
#     # Pad with right zeros to reach the sequence length (no_fragments + fragment_size)
#     structure_feature_fragments = np.pad(structure_feature_fragments, ((0, 0), (0, structure_feature_fragments.shape[0])))
#
#     # Roll the fragments one by one (add heading zeros)
#     for i, row in enumerate(structure_feature_fragments):
#         structure_feature_fragments[i] = np.roll(row, i)
#
#     # Calculate average along columns of overlapping fragments (average RMSD per residue)
#     structure_feature_average = np.average(structure_feature_fragments, axis=0)
#
#
#     #PYMOL SCRIPT
#
#     cmd.load("data/pdb{}.pdb".format(j), j)  # Load from file
#     cmd.remove("resn hoh")  # Remove water molecules
#     cmd.hide("lines", "all")  # Hide lines
#     cmd.show("cartoon", j)  # Show cartoon
#     norm = colors.Normalize(vmin=min(structure_feature_average), vmax=max(structure_feature_average))
#     for i, residue in enumerate(Selection.unfold_entities(structure[0], "R")):
#         rgb = cm.bwr(norm(structure_feature_average[i]))
#         # print(i, residue.id, structure_rmsd_average[i], rgb)
#         cmd.set_color("col_{}".format(i), list(rgb)[:3])
#         cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))
#     cmd.png("data/pymol_image", width=2000, height=2000, ray=1)
#     h=j