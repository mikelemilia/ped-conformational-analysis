import re
import sys

import numpy as np

from ModelFeatures import ModelFeatures
from PedFeatures import PedFeatures
from utils import *


class Menu:

    def do_1(self):
        """Task 1 : Analyze features of models inside a PED"""

        inputs = ask_input()

        if inputs is not None:
            folder = inputs[0]
            ped_name = inputs[1]
            self.first_task(folder, ped_name)

    def do_2(self):
        """Task 2 : Compare inter and intra features of a PED"""

        inputs = ask_input()

        if inputs is not None:
            folder = inputs[0]
            ped_name = inputs[1]
            self.second_task(folder, ped_name)

    @staticmethod
    def do_3():
        """Exit"""
        print('Exiting...')

    def first_task(self, folder, ped_name):
        """
        First task of the project
        :param pdb_id_list: list of the pdb id of one ped
        :return: null
        """

        pdb_id_list = extract_filenames(folder, ped_name, ['pdb', 'ent'])
        for pdb_id in pdb_id_list:
            print("\nAnalyzing {}...\n".format(pdb_id))

            model_features = ModelFeatures(folder, pdb_id)
            model_features.choice_maker()

            model_features.compute_clustering()
            graph = model_features.generate_graph()
            model_features.generate_pymol_image(graph)

    def second_task(self, folder, ped_name):
        """
        Second task of the project
        :param ped_name: name of ped to analyse
        :return: null
        """
        print('\nAnalyzing {}'.format(ped_name))
        ped_features = PedFeatures(folder, ped_name)
        ped_features.choice_maker()
        ped_features.global_dendrogram()
        ped_features.global_heatmap()
        ped_features.local_metric()

    def run(self):

        user_input = 0

        self.generate_menu()

        while user_input != 3:
            user_input = input("\nSelect what do you want to do: ")

            if user_input.isdigit():
                user_input = int(user_input)
                self.execute(user_input)
            else:
                print('WARNING : An invalid choice has been identified. Please enter a valid one')

        print('\nProgram stopped.')

    def execute(self, user_input):

        controller_name = 'do_{}'.format(user_input)

        try:

            controller = getattr(self, controller_name)

        except AttributeError:

            print('WARNING : Selected option not found')

        else:

            controller()

    def generate_menu(self):
        print('== == == == == == == == == == == == == == == == == == ==')
        do_methods = [m for m in dir(self) if m.startswith('do_')]
        menu_string = "\n".join(
            ['{}.    {}'.format(method[-1], getattr(self, method).__doc__) for method in do_methods])
        print(menu_string)
        print('== == == == == == == == == == == == == == == == == == ==')
        print('You can reverse any wrong selection with \n\t- Q (QUIT)\n\t- B (BACK)\n\t- U (UNDO)')
        print('== == == == == == == == == == == == == == == == == == ==')
