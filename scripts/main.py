import math
import sys
from Model_features import Model_Features
from Clustering import Clustering
from utils import *
from PED_features import *

if __name__ == "__main__":

    folder, ped_name = parser()

    if not os.path.isdir(folder):
        print('The path specified does not exist')
        sys.exit()
    else:

        pdb_id_list = extract_filenames(folder, ped_name, "pdb")

        for pdb_id in pdb_id_list:

            print("Analyzing {}...".format(pdb_id))

            model_features = Model_Features(folder, pdb_id)
            feat = model_features.choice_maker()

            # print('\nClustering...')
            # test = Clustering(pdb_id, feat)
            #
            # test.compute_clustering()
            # test.generate_graph()

            print('\n------------------\n')

        print('Task2...')

        ped_features = PED_features(folder, ped_name)
        ped_features.choice_maker()

        # ATTENZIONE: alla fine di ped_feat, ultima riga, ci sono dei Nan!!
