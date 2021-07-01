import os
from Menu import Menu

if __name__ == "__main__":

    # Set of the working directory to the project folder
    os.chdir(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    # Menu execution
    menu = Menu()
    menu.run()
