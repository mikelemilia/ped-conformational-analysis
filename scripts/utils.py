import argparse
import glob
import os


def parser():
    # Create the parser
    my_parser = argparse.ArgumentParser(description='PED tool')

    # Add the arguments
    my_parser.add_argument('path', metavar='path', type=str, help='the path to input folder')

    # Execute the parse_args() method
    args = my_parser.parse_args()

    folder = args.path

    return folder


def extract_filenames(folder, name="", extensions=None):

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
