import os
import sys
from Menu import Menu

if __name__ == "__main__":

    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    menu = Menu()
    menu.run()
