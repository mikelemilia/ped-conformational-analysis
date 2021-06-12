import re
import sys

from ModelFeatures import ModelFeatures
from PedFeatures import PedFeatures
from utils import extract_filenames


class Menu:

    def __init__(self, folder):
        self._folder = folder

    def do_1(self):
        """Task 1 : Analyze features of models inside a PED"""

        exit_menu = False

        while not exit_menu:

            ped_name = input("Insert the name of the PED you want to analyze: ")

            if ped_name == 'Q':

                exit_menu = True

            else:

                if not re.match('^PED[0-9]{5}$', ped_name):

                    exit_menu = False

                    print("You provide a wrong PED name, it should match the format PEDxxxxx\n", file=sys.stderr)
                    sys.stderr.flush()

                else:

                    pdb_id_list = extract_filenames(self._folder, ped_name, ['pdb', 'ent'])

                    if len(pdb_id_list) == 0:

                        exit_menu = False

                        print('PED name not found. You should retry with another PED\n', file=sys.stderr)
                        sys.stderr.flush()

                    else:

                        exit_menu = True
                        self.first_task(pdb_id_list)

    def do_2(self):
        """Task 2 : Compare inter and intra features of a PED"""

        exit_menu = False

        while not exit_menu:

            ped_name = input("Insert the name of the PED you want to analyze: ")

            if ped_name == 'Q':

                exit_menu = True

            else:

                if not re.match('^PED[0-9]{5}$', ped_name):

                    exit_menu = False

                    print("You provide a wrong PED name, it should match the format PEDxxxxx\n", file=sys.stderr)
                    sys.stderr.flush()

                else:

                    exit_menu = True

                    self.second_task(ped_name)

    @staticmethod
    def do_3():
        """Exit"""
        print('Exiting...')

    def first_task(self, pdb_id_list):

        for pdb_id in pdb_id_list:
            print("\nAnalyzing {}...".format(pdb_id))

            model_features = ModelFeatures(self._folder, pdb_id)
            model_features.choice_maker()

            # model_features.compute_clustering()
            # graph = model_features.generate_graph()
            # model_features.generate_pymol_img(graph)

    def second_task(self, ped_name):

        ped_features = PedFeatures(self._folder, ped_name)
        x = ped_features.choice_maker()

        if x == -1:

            print('Model features not found. You must run Task 1 in advance for this specific PED\n', file=sys.stderr)
            sys.stderr.flush()

        else:

            print('\n\t- Comparing {}'.format(ped_name))
            ped_features.global_dendrogram()
            ped_features.global_heatmap()
            ped_features.distance_matrix_med_rmsd_peds()
            ped_features.local_metric()

    def run(self):
        user_input = 0
        self.generate_menu()
        while user_input != 3:
            user_input = int(input("\nSelect what do you want to do: "))
            self.execute(user_input)
        print('Program stopped.')

    def execute(self, user_input):
        controller_name = 'do_{}'.format(user_input)
        try:
            controller = getattr(self, controller_name)
        except AttributeError:
            print('Method not found')
        else:
            controller()

    def generate_menu(self):
        print('== == == == == == == == == == == == == == == == == == ==')
        do_methods = [m for m in dir(self) if m.startswith('do_')]
        menu_string = "\n".join(
            ['{}.    {}'.format(method[-1], getattr(self, method).__doc__) for method in do_methods])
        print(menu_string)
        print('== == == == == == == == == == == == == == == == == == ==')
        print('Undo any selection with Q')
        print('== == == == == == == == == == == == == == == == == == ==')
