# Import packages
from rdkit import Chem
from rdkit.Chem.Draw import MolsToGridImage
from rdkit.Chem.Draw import MolToFile
from rdkit.Chem.rdmolfiles import SmilesMolSupplier
from pathlib import Path
import os

# Convert dataset

## Convert smiles to images
# SOURCE
paths = Path('/home/jacobnorth/Documents/GitHub/mdcode/ml_models/data/ml_2020_09_03_09_50_00').glob('**/*.smi')

attempt = 1

# DESTINATION
destdirname = ("drw_attempt"+str(attempt))      # name the dir
os.mkdir(destdirname)       # make a dir to hold each object

# make image files at the specified resolution

idx = 0         # start at 0 for each image name

for path in paths:
    path_in_str = str(path)
    # CREATE ALL IMAGES FOR THE GIVEN FILESET
    suppl = SmilesMolSupplier(path_in_str)      # create a SMILES supplier from the current image
    for mol in suppl:       # for each molecule specified in the supplier,
        filename = ("drw_"+str(attempt)+"_"+str(idx)+".png")     # name the new file
        MolToFile(mol, filename, size=(128,128))    # create the new file
        idx=idx+1       # increase the idx
        os.rename(filename, destdirname+"/"+filename)   # move the file into the dest directory
attempt=attempt+1