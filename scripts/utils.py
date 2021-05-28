import argparse
import glob
import os
import re


def parser():
    # Create the parser
    my_parser = argparse.ArgumentParser(description='PED tool')

    # Add the arguments
    my_parser.add_argument('path', metavar='path', type=str, help='the path to input folder')
    my_parser.add_argument('ped_name', metavar='ped_name', type=str, help='ID of the PED to be used (e.g. PEDxxxxx)')

    # Execute the parse_args() method
    args = my_parser.parse_args()

    folder = args.path
    ped_name = args.ped_name

    if not re.match('^PED[0-9]{5}$', ped_name):
        print("Wrong PED name")
        exit(1)

    return folder, ped_name


def extract_filenames(folder, name="", extension="pdb"):
    files = []

    paths = glob.glob(folder + "/" + name + "*." + extension)
    for path in paths:
        base = os.path.basename(path)
        file = os.path.splitext(base)[0]
        files.append(file)

    return sorted(files)


# # HSE exposure
# # https://en.wikipedia.org/wiki/Half_sphere_exposure
# # Low HSE-up means high ASA (absolute solvent accessibility)
# hse = HSExposureCB(structure[0])
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             if residue.id[0] == " ":
#                 print(chain.id, residue.id, hse[(chain.id, residue.id)], dssp_dict.get((chain.id, residue.id)))  # HSE beta up (EXP_HSE_B_U), HSE beta down (EXP_HSE_B_D)

