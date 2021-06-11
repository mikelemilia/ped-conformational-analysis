import os
import sys

from utils import parser
from Menu import Menu

if __name__ == "__main__":

    folder = parser()

    if not os.path.isdir(folder):
        print('You must provide an existing data location')
        sys.exit()

    menu = Menu(folder)
    menu.run()
