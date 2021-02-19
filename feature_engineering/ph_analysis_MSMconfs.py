
'''
IMPORTS
'''
import pdbutils
import structuralPH
import scipy
import numpy as np
from ripser import ripser
from persim import plot_diagrams
import matplotlib.pyplot as plt
import Bio

'''
FUNCTIONS
'''

# Setup structure
def getStrucAtomPositions(filename:str) -> np.array:
    '''
    Extract the 3d positions of all atoms in a structure.
    '''
    ### SET UP THE STRUCTURE OBJECT -- A NUMPY ARRAY WITH ATOM POSITIONS
    struc = pdbutils.createStruct("6VXX.pdb")               # get the structure
    coords = elemSpecificPH.getPositionArray(struc)         # separate it by its elements
    strucarray = elemSpecificPH.struc2elemPosArrays(struc)  # construct the structure array
    a = strucarray[1]                                       # take the first element from the structure array
    #a = scipy.spatial.distance.pdist(strucarray)
    return a

def runPH(positions:np.array) -> list:
    '''
    Run persistent homology calculations and return the objects from Ripser here.
    '''
    diagrams = ripser(a)['dgms']                            # compute the persistence diagram
    plot_diagrams(diagrams, show=False)                     # plot the diagrams
    plt.savefig("persistent_homology.png")                  # save the figure
    return diagrams


