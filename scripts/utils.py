import glob
import os


def extract_filenames(folder, extension):
    files = []

    paths = glob.glob(folder + "/*." + extension)
    for path in paths:
        base = os.path.basename(path)
        file = os.path.splitext(base)[0]
        files.append(file)

    return sorted(files)