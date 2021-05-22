import sys
from Features import Features
from Clustering import Clustering
from utils import *

# TODO: resolve problems of dataset PED00142

if __name__ == "__main__":

    folder = parser()

    if not os.path.isdir(folder):
        print('The path specified does not exist')
        sys.exit()
    else:
        os.makedirs("{}/features_files".format(folder), exist_ok=True)
        pdb_id_list = extract_filenames(folder, "pdb")

        for pdb_id in pdb_id_list:

            print("Analyzing {}...".format(pdb_id))
            path = os.path.join(folder, pdb_id + '.pdb')

            features = Features(path, pdb_id)

            path_features = "{}/features_files/{}_features.csv".format(folder, pdb_id)

            if os.path.exists(path_features):
                print('\nLoading features...')
                feat = features.extract_features(path_features)
            else:
                print('\nComputing features...')
                feat = features.compute_features()
                features.save(path_features)

            print('\nClustering...')
            test = Clustering(pdb_id, feat)

            test.compute_clustering()
            test.generate_graph()

            print('\n------------------\n')

