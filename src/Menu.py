from ModelFeatures import ModelFeatures
from PedFeatures import PedFeatures
from utils import *


class Menu:

    def do_1(self):
        """Task 1 : Analyze features of models inside a PED"""

        # Input management
        inputs = ask_input()

        # If provided input is not Q, B, U, extracts folder and ped_name and performs first task
        if inputs is not None:
            folder = inputs[0]
            ped_name = inputs[1]
            self.first_task(folder, ped_name)

    def do_2(self):
        """Task 2 : Compare inter and intra features of a PED"""

        # Input management
        inputs = ask_input()

        # If provided input is not Q, B, U, extracts folder and ped_name and performs second task
        if inputs is not None:
            folder = inputs[0]
            ped_name = inputs[1]
            self.second_task(folder, ped_name)

    @staticmethod
    def do_3():
        """Exit"""
        # Method for exiting the program
        print('Exiting...')

    def first_task(self, folder, ped_name):
        """
        First task of the project
        :param folder: path of the folder containing PDB files of the ensembles with name ped_name
        :param ped_name: name of the PED under analysis
        """

        # Extract PED ID (PEDxxxxxexx) list
        pdb_id_list = extract_filenames(folder, ped_name, 'pdb')

        # For each PED ID, perform the first task
        for pdb_id in pdb_id_list:

            print("\nAnalyzing {}...\n".format(pdb_id))
            model_features = ModelFeatures(folder, pdb_id)

            # Compute or load the features
            model_features.choice_maker()

            # Requested data analysis
            model_features.compute_clustering()
            graph = model_features.generate_graph()
            model_features.generate_pymol_image(graph)

    def second_task(self, folder, ped_name):
        """
        Second task of the project
        :param folder: path of the folder containing PDB files of the ensembles with name ped_name
        :param ped_name: name of the PED under analysis
        """

        print('\nAnalyzing {}'.format(ped_name))
        ped_features = PedFeatures(folder, ped_name)

        # Compute or load the features
        ped_features.choice_maker()

        # Requested data analysis
        ped_features.global_dendrogram()
        ped_features.global_heatmap()
        ped_features.local_metric()

    def run(self):
        """
        Function to generate the menu, allow to user to make his choice and generate the action relative to that choice.
        """

        user_input = 0
        self.generate_menu()

        # Repeat until the user doesn't provide value 3 (to exit)
        while user_input != 3:
            user_input = input("\nSelect what do you want to do: ")

            # Parsing of the input provided
            if user_input.isdigit():
                user_input = int(user_input)
                self.execute(user_input)
            else:
                print('WARNING : An invalid choice has been identified. Please enter a valid one')

        print('\nProgram stopped.')

    def execute(self, user_input):
        """
        Function to execute the do_{} function chosen by the user with the corresponding check of existance.
        :param user_input: user input of do_{} function
        """

        controller_name = 'do_{}'.format(user_input)

        try:
            controller = getattr(self, controller_name)
        except AttributeError:
            print('WARNING : Selected option not found')
        else:
            controller()

    def generate_menu(self):
        """
        Print of the corresponding menu based on the do_{} function provided in the current class.
        """

        print('== == == == == == == == == == == == == == == == == == ==')
        do_methods = [m for m in dir(self) if m.startswith('do_')]
        menu_string = "\n".join(
            ['{}.    {}'.format(method[-1], getattr(self, method).__doc__) for method in do_methods])
        print(menu_string)
        print('== == == == == == == == == == == == == == == == == == ==')
        print('You can reverse any wrong selection with \n\t- Q (QUIT)\n\t- B (BACK)\n\t- U (UNDO)')
        print('== == == == == == == == == == == == == == == == == == ==')
