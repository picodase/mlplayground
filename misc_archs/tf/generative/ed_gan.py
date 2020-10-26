## ed_gan
### A GAN for generation of protein structures inside an electron density map

# %% 
# import packages
import tensorflow as tf
import Bio.PDB 
import scipy

# %% 
# take in a protein structure

# convert it to an electron density map with randomized uncertainty
# https://docs.scipy.org/doc/scipy/reference/tutorial/fft.html

# 