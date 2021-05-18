import argparse
import sys
from Features import Features
from Clustering import Clustering

from utils import *


if __name__ == "__main__":

    # TODO : move the parser inside an utility function

    # Create the parser
    my_parser = argparse.ArgumentParser(description='PED tool')

    # Add the arguments
    my_parser.add_argument('path',
                           metavar='path',
                           type=str,
                           help='the path to input folder')

    # Execute the parse_args() method
    args = my_parser.parse_args()

    folder = args.path

    if not os.path.isdir(folder):
        print('The path specified does not exist')
        sys.exit()
    else:
        os.makedirs("{}/features_files".format(folder), exist_ok=True)
        pdb_id_list = extract_filenames(folder, "pdb")

        for pdb_id in pdb_id_list:

            path = os.path.join(folder, pdb_id + '.pdb')

            features = Features(path, pdb_id)

            path_features = "{}/features_files/{}_features.csv".format(folder, pdb_id)

            if os.path.exists(path_features):
                feat = features.extract_features(path_features)
            else:
                feat = features.compute_features()
                features.save(path_features)

            print(feat)
            test = Clustering(feat)

            test.compute_clustering()
            test.generate_graph()

            print('------------------\n')

