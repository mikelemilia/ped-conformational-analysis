from ModelFeatures import *
import pymol
from pymol import cgo, cmd, util



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