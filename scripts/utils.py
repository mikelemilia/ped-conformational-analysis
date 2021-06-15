import argparse
import glob
import os

def parser():
    """
    This function allows you a parser for the input path file
    :return: folder
    """

    # Create the parser
    my_parser = argparse.ArgumentParser(description='PED tool')

    # Add the arguments
    my_parser.add_argument('-p', '--path', metavar='path', type=str, help='the path to input folder', default='data')

    # Execute the parse_args() method
    args = my_parser.parse_args()

    folder = args.path

    return folder


def extract_filenames(folder, name="", extensions=None):
    """

    :param folder:
    :param name:
    :param extensions:
    :return:
    """
    if extensions is None:
        extensions = ['pdb', 'ent']

    files = []

    for extension in extensions:

        paths = glob.glob(folder + "/" + name + "*." + extension)

        for path in paths:
            base = os.path.basename(path)
            file = os.path.splitext(base)[0]
            files.append(file)

    return sorted(files)


# def pincopalla():
#     import pymol
#     from pymol import cgo, cmd, util
#     from pymol2 import PyMOL
#
#     # pymol.finish_launching()  # Open Pymol
#
#     p = PyMOL()
#     p.start()
#     # Input 1jsu (3 chains, non-standard amino acids), 1az5 (disordered loops, chain breaks)
#     pdb_id = '1jsu'
#
#     cmd.fetch(pdb_id, pdb_id, path="data/")  # Download the PDB
#     # cmd.load("/data/PED00020e001.pdb")  # Load from file
#
#     cmd.remove("resn hoh")  # Remove water molecules
#     cmd.hide("lines", "all")  # Hide lines
#     cmd.show("cartoon", pdb_id)  # Show cartoon
#     cmd.show("sticks", "hetatm")  # Show hetero atoms as sticks
#     # # cmd.spectrum(selection="all")  # Rainbow color
#     util.cbc(selection="all")  # Color by chain
#     util.cnc(selection="all")  # Color by atom type, but not the C atoms
#
#     # Select and color two residues
#     sele_name = "nodes"
#     res1 = 'B/200/'
#     res2 = 'C/52/'
#     cmd.select(sele_name, '{} or {}'.format(res1, res2))
#     cmd.show("spheres", sele_name)
#     cmd.set('sphere_transparency', 0.5, sele_name)  # Set transparency
#
#     # Get coordinates of two atoms
#     atom1 = 'B/200/SD'
#     atom2 = 'C/52/CE'
#     coord1 = cmd.get_coords(atom1, 1)  # shape N * 3, where N is the number of atoms in the selection
#     coord2 = cmd.get_coords(atom2, 1)
#     # model = cmd.get_model(pdb_id, 1)  # Fastest way to get all atom coordinates
#
#     # Calculate center of mass between two residues and create a new "pseudo" atom
#     center_of_mass = (coord1[0] + coord2[0]) / 2
#     print(coord1, coord2, center_of_mass)
#
#     obj_name = "ps_atom"
#     cmd.pseudoatom(obj_name, pos=list(center_of_mass))
#     cmd.show("spheres", obj_name)
#     # cmd.extract(...  # Move selected atoms to a new object
#     # cmd.create(...  # Create a new object from selection
#
#     cr, cg, cb = (1.0, 0.0, 0.0)  # RGB red
#
#     # Create lines object
#     obj = [cgo.BEGIN, cgo.LINES, cgo.COLOR, cr, cg, cb]
#     obj.append(cgo.VERTEX)
#     obj.extend(list(coord1[0]))
#     obj.append(cgo.VERTEX)
#     obj.extend(list(coord2[0]))
#     obj.append(cgo.END)
#
#     # Set the object
#     obj_name = 'edges'
#     cmd.load_cgo(obj, obj_name)
#     cmd.set("cgo_line_width", float(3), obj_name)
#
#     cmd.orient(pdb_id)  # Set the origin to full protein
#
#     input("Stupido")
#     # pymol.cmd.quit()
#     cmd.save('output/stupido.png')
#     cmd.save('output/stupido.mol')
#     p.stop()
