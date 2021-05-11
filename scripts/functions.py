import numpy as np
import math

def radius_gyration(chain):
    '''
    Calculates the Radius of Gyration (Rg) of a protein in Angstroms.
    Does not use mass and assume heavy atoms have the same mass.

    https://en.wikipedia.org/wiki/Radius_of_gyration  (formula considering mass)
    https://link.springer.com/article/10.1134/S0026893308040195  (formula without mass)
    '''

    # Heavy atoms coordinates
    coord = list()
    for atom in chain.get_atoms():
        if atom.get_name()[0] in ['C', 'O', 'N', 'S']:
            coord.append(atom.get_coord())
    coord = np.array(coord)  # N X 3

    barycenter = np.sum(coord, axis=0) / coord.shape[0]  # center of mass is more correct

    # Calculate distance of each atom from the barycenter
    dist = coord - barycenter
    dist = dist * dist
    dist = np.sqrt(np.sum(dist, axis=1))
    print(dist)

    return round(math.sqrt(np.sum(dist * dist) / len(coord)), 3)
