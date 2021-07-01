import glob
import os
import re

from numpy import unique


def extract_names(folder):
    """
    Function to extract the PED ID (format PEDxxxxx) of the files present in the folder passed as input and return
    a list of unique ID founds
    :param folder: folder path where to find PED IDs
    :return: list of unique PED ID
    """

    # Extract all the filenames
    files = sorted([s for s in os.listdir(folder) if os.path.isfile(os.path.join(folder, s))])
    # Maintain only the files names of format PEDxxxxxexxx.pdb
    names = [s for s in files if re.match('^PED[0-9]{5}e[0-9]{3}.pdb$', s)]
    # Split the filename in order to maintain only the PED ID
    ped_names = [s.split('e')[0] for s in names]

    return list(unique(ped_names))


def check_path(folder, name=""):
    """
    Function to check that the provided data folder is correct: it must be an existing folder and contain at least one
    pdb file for the analysis
    :param folder: data folder path to be checked
    :param name: eventual name of PED ID to be checked inside the folder
    :return: True if the folder is valid, False otherwise
    """

    # Check existence
    if not os.path.isdir(folder):
        print('You must provide an existing data location')
        return False

    # Check that it contains at least one pdb file
    pdb_id_list = extract_filenames(folder, name, 'pdb')
    if len(pdb_id_list) == 0:
        print('WARNING : PED name not found. You should retry with another PED\n')
        return False

    return True


def check_input(choice: str, choices: list):
    """
    Function to check that the provided input is correct, so if it belongs to the list of possible choices. Note that
    the user could digit the index relative to the choice or the choice itself.
    :param choice: provided input to be checked
    :param choices: list of possible choices that the user could do
    :return: True/False depending on the validity of the input and corresponding choice inside the list
    """
    valid = False

    # Check user selection
    if choice.isdigit():
        # If the user selects an index, check that is a valid one (in the list of choices index)
        try:
            choice = choices[int(choice)]
            print('Your choice is {}'.format(choice))
            valid = True
        except IndexError:
            print("WARNING : Wrong choice, unable to find it.")
    else:
        # If the user selects a choice, check that is a valid one (in the list of choices)
        if choice.upper() not in choices:
            print("WARNING : Wrong choice, unable to find it.")
        else:
            print('Your choice is {}'.format(choice.upper()))
            valid = True

    return valid, choice


def ask_input():
    """
    Function for the implementation of the input providing with the corresponding checks. For the two tasks, it is
    necessary to have the data folder (to be checked), to list the PED ID of the pdb files found inside it and to
    allow the user to (correctly) select one of them to be analyzed. The possibility to have undo operations from the
    user must be taken into consideration
    :return: tuple containing folder path and PED name, None if the user decides to undo
    """

    folder_found = False
    folder = ""
    ped_name = ""
    # Accepted undo operations
    _undo_operation = [
        'Q', 'QUIT',
        'B', 'BACK',
        'U', 'UNDO'
    ]

    # Request for the data folder and check of its existence (if undo operation, return to the menu)
    while not folder_found:
        folder = input("Insert the folder path containing the PED you want to analyze: ")
        if folder in _undo_operation:
            return None
        if not os.path.isdir(folder):
            print('You must provide an existing data location')
        else:
            folder_found = True

    exit_menu = False
    while not exit_menu:

        # Extract PED ID of files contained in the folder: if no one is present, return to the menu
        ped_names = extract_names(folder)
        if len(ped_names) == 0:
            print('WARNING : Folder not containing PED files!')
            return None

        # Request for the PED to be analyzed showing the ones present in the folder
        print('\nWhich PED do you want to analyze: ')
        for i, names in enumerate(ped_names):
            print('\t{} - {}'.format(i, names))

        ped_name = input('\nYour choice: ')

        # Check of the choice
        if ped_name in _undo_operation:
            exit_menu = True
        else:
            exit_menu, ped_name = check_input(ped_name, ped_names)

    return [folder, ped_name]


def extract_filenames(folder, name="", extension='pdb'):
    """
    Extract the list of filenames inside a folder (facultative the matching with the name provided) of a certain extension.
    :param folder: folder path in which look for the files
    :param name: file name to be matched
    :param extension: extension of the files to be found
    :return: list of filenames found
    """

    files = []
    # Extract all the paths corresponding to the request
    paths = glob.glob(folder + "/" + name + "*." + extension)

    # Extract only the filenames (without parent folder path and extension) from each filename
    for path in paths:
        base = os.path.basename(path)
        file = os.path.splitext(base)[0]
        files.append(file)

    return sorted(files)
