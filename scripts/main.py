import sys
from Model_features import Model_Features
from Clustering import Clustering
from utils import *
from PED_features import *

# TODO: resolve problems of dataset PED00142

if __name__ == "__main__":

    folder, ped_name = parser()

    if not os.path.isdir(folder):
        print('The path specified does not exist')
        sys.exit()
    else:

        model_features_folder = folder + '/model_features/'
        os.makedirs(model_features_folder, exist_ok=True)
        pdb_id_list = extract_filenames(folder, ped_name, "pdb")

        for pdb_id in pdb_id_list:

            print("Analyzing {}...".format(pdb_id))
            path = os.path.join(folder, pdb_id + '.pdb')

            model_features = Model_Features(path, pdb_id)

            path_features = "{}/{}_features.csv".format(model_features_folder, pdb_id)

            if os.path.exists(path_features):
                print('\nLoading features...')
                feat = model_features.extract(path_features)
            else:
                print('\nComputing features...')
                feat = model_features.compute()
                model_features.save(path_features)

            # print('\nClustering...')
            # test = Clustering(pdb_id, feat)
            #
            # test.compute_clustering()
            # test.generate_graph()

            print('\n------------------\n')

        print('Task2...')

        ped_features = PED_features(model_features_folder, ped_name)
        ped_features_folder = folder + '/ped_features/'

        path_ped_features = "{}/{}_features.csv".format(ped_features_folder, ped_name)

        if os.path.exists(path_ped_features):
            print('\nLoading features comparison...')
            ped_feat = ped_features.extract(path_ped_features)
            print(ped_feat)
        else:
            print('\nComparing features...')
            ped_feat = ped_features.compare()
            ped_features.save(path_ped_features)

        # ATTENZIONE: alla fine di ped_feat, ultima riga, ci sono dei Nan!!
