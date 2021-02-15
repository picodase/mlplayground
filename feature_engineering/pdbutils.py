

'''
IMPORTS
'''

import Bio.PDB


'''
FUNCTIONS
'''


def createStruct(filename: str) -> Bio.PDB.Structure.Structure:
    '''
    Creates and returns a structure object from a PDB file, given as the filename/filepath (should work for both!).
    '''
    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser(PERMISSIVE=1)

    structure_id = filename.replace(".pdb", "")     # replace the file extension with nothing, this is the structure ID!
    return parser.get_structure(structure_id, filename)     # return the structure
