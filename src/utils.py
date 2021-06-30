import glob
import os
import re

from numpy import unique


def extract_names(folder):
    files = sorted([s for s in os.listdir(folder) if os.path.isfile(os.path.join(folder, s))])
    names = [s for s in files if re.match('^PED[0-9]{5}e[0-9]{3}.pdb$', s)]
    ped_names = [s.split('e')[0] for s in names]
    return list(unique(ped_names))


def check_path(folder, name=""):

    if not os.path.isdir(folder):
        print('You must provide an existing data location')
        return False

    pdb_id_list = extract_filenames(folder, name, ['pdb', 'ent'])

    if len(pdb_id_list) == 0:
        print('WARNING : PED name not found. You should retry with another PED\n')
        return False

    return True


def check_input(choice: str, choices: list):
    valid = False

    # Check user selection
    if choice.isdigit():
        try:
            choice = choices[int(choice)]
            print('Your choice is {}'.format(choice))
            valid = True
        except IndexError:
            print("WARNING : Wrong choice, unable to find it.")
    else:
        if choice.upper() not in choices:
            print("WARNING : Wrong choice, unable to find it.")
        else:
            print('Your choice is {}'.format(choice.upper()))
            valid = True

    return valid, choice


def ask_input():
    folder_found = False
    folder = ""
    ped_name = ""
    _undo_operation = [
        'Q', 'QUIT',
        'B', 'BACK',
        'U', 'UNDO'
    ]

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

        ped_names = extract_names(folder)
        if len(ped_names) == 0:
            print('WARNING : Folder not containing PED files!')
            return None

        print('\nWhich PED do you want to analyze: ')
        for i, names in enumerate(ped_names):
            print('\t{} - {}'.format(i, names))

        ped_name = input('\nYour choice: ')

        if ped_name in _undo_operation:
            exit_menu = True
        else:
            exit_menu, ped_name = check_input(ped_name, ped_names)

    return [folder, ped_name]


def extract_filenames(folder, name="", extensions='pdb'):

    files = []

    for extension in extensions:

        paths = glob.glob(folder + "/" + name + "*." + extension)

        for path in paths:
            base = os.path.basename(path)
            file = os.path.splitext(base)[0]
            files.append(file)

    return sorted(files)
