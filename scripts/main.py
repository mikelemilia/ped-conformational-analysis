import sys
from Model_features import Model_Features
from Clustering import Clustering
from utils import *
from PED_features import *

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

            model_features = Model_Features(path, pdb_id)

            path_features = "{}/features_files/{}_features.csv".format(folder, pdb_id)

            if os.path.exists(path_features):
                print('\nLoading features...')
                feat = model_features.extract(path_features)
            else:
                print('\nComputing features...')
                feat = model_features.compute()
                model_features.save(path_features)

            print('\nClustering...')
            test = Clustering(pdb_id, feat)

            test.compute_clustering()
            test.generate_graph()

            print('\n------------------\n')

        print('Task2...')
        ped_features = PED_features(folder+'/features_files')

        path_PED_features = "{}/features_files/{}_features.csv".format(folder, ped_features.get_ped_name())

        if os.path.exists(path_PED_features):
            print('\nLoading features comparison...')
            ped_feat = ped_features.extract(path_PED_features)
            print(ped_feat)
        else:
            print('\nComparing features...')
            ped_feat = ped_features.compare()
            ped_features.save(path_PED_features)

        # ATTENZIONE: alla fine di ped_feat, ultima riga, ci sono dei Nan!!
