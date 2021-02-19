'''
IMPORTS
'''
import Bio.PDB
import numpy as np
import pdbutils

'''
FUNCTIONS
'''
def sepEN(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                for atom in res:    # for each atom,
                    if atom.get_id() == "N": # if the atom is an carbon,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

carbontypes = {"C","CA","CB","CG","CE"}

def sepEC(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                for atom in res:    # for each atom,
                    if atom.get_id() in carbontypes: # if the atom is a carbon,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

def sepEO(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                for atom in res:    # for each atom,
                    if atom.get_id() == "O": # if the atom is an oxygen,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

def sepES(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                for atom in res:    # for each atom,
                    if atom.get_id() == "S": # if the atom is a sulfur,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

###### AMINO ACID TYPES ######

hydrophobics = {"PHE","ALA","LEU","MET","ILE","TRP","PRO","VAL"}

def sepRHψ(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in hydrophobics: # if residue is in the class of hydrophobic aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

polarunchargeds = {"CYS","GLY","GLN","ASN","SER","TYR","THR"}

def sepRδ0(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in polarunchargeds:  # if residue is in the class of polar aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray


acidics = {"ASP","GLU"}

def sepRδa(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in acidics:  # if residue is in the class of polar aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

basics = {"HIS","ARG","LYS"}

def sepRδb(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in basics:   # if residue is in the class of polar aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

larges = {"ARG","LYS","TRP","MET","GLU","GLN","TYR","PHE"}

def sepRl(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in larges:   # if residue is in the class of large aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

smalls = {"HIS","ILE","LEU","ALA","PRO","VAL","CYS","GLY","SER","ASN","THR","ASP"}

def sepRs(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in smalls: # if residue is in the class of small aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray


###### TYPE ######

specials = {"GLY", "CYS", "PRO"}

def sepTs(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    '''
    Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    '''
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res.get_resname() in specials: # if residue is in the class of small aa's,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

def runAllSeps(struc: Bio.PDB.Structure.Structure) -> list:
    '''
    Runs all previously-defined seps and returns a list of them in order. Sep types are separated by sublist.
    '''
    elemtype = [sepEC(struc), sepEN(struc), sepEO(struc), sepES(struc)]
    aatype = [sepRHψ(struc),sepRl(struc),sepRs(struc),sepRδ0(struc),sepRδa(struc),sepRδb(struc)]
    sstype = [sepTs(struc)]

    return [elemtype,aatype,sstype]

###### SECONDARY STRUCTURE ######
'''
def sssphβ(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    #Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res # if residue is in a beta-sheet,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray

# modify the sep function

def sssphα(struc: Bio.PDB.Structure.Structure) -> np.ndarray:
    #Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res # if residue is in an alpha-helix,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray


def sssphLoop(struc: Bio.PDB.Structure.Structure) -> np.ndarray:

    #Obtains a position array for the 3D-coordinates of each atom in the provided structure.
    coords = []   # construct a (3 x natoms) array
    for model in struc: # for each model,
        for chain in model:    # for each chain,
            for res in chain:   # for each residue,
                if res # if residue is in a loop,
                    for atom in res:    # for each atom,
                        coords.append(atom.get_coord()) # get the coordinates, append to a new array
    return np.array(coords)    # convert to np.ndarray
    '''
