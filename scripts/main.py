import math
import sys
from ModelFeatures import ModelFeatures
from ModelComparisons import *
from utils import *
from PedFeatures import *
from PedComparisons import *

if __name__ == "__main__":

    folder, ped_names = parser()

    if not os.path.isdir(folder):
        print('The path specified does not exist')
        sys.exit()
    else:

        n_residues = []
        n_conformations = []

        for ped_name in ped_names:

            pdb_id_list = extract_filenames(folder, ped_name, "pdb")

            print('Task1...')

            # for pdb_id in pdb_id_list:
            #
            #     print("Analyzing {}...".format(pdb_id))
            #
            #     model_features = ModelFeatures(folder, pdb_id)
            #     feat = model_features.choice_maker()
            #
            #     print("Clustering {}...".format(pdb_id))
            #
            #     model_features.compute_clustering()
            #     graph = model_features.generate_graph()
            #     # model_features.generate_pymol_img(graph)
            #
            #     print('\n------------------\n')

            print('Task2...')

            ped_features = PedFeatures(folder, ped_name)
            ped_feat = ped_features.choice_maker()

            n_residues = ped_features.num_residues
            n_conformations = ped_features.num_conformations

            # ATTENZIONE: alla fine di ped_feat, ultima riga, ci sono dei Nan!!

            print('\nComparison between PEDs')
            comparison = PedComparison(ped_feat, n_conformations, n_residues)
            comparison.global_dendrogram()
            comparison.global_heatmap()
            # comparison.distance_matrix_med_rmsd_peds()
            comparison.local_metric()
