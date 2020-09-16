# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %% [markdown]
# # tf_DCGAN_2Dchemim
# -------
# ### Tensorflow Deep convolutional GAN using 2D images of chemical structures
# %% [markdown]
# ## Setup
# ### Packages

# %%
# import needed packages
import glob
import imageio
import matplotlib.pyplot as plt
import numpy as np
import os
import PIL
from tensorflow.keras import layers
import time
import pathlib

from IPython import get_ipython
from IPython import display


# %%
import tensorflow as tf

# ## Use `imageio` to create an animated gif using the images saved during training.

# %%
anim_file = 'dcgan.gif'

with imageio.get_writer(anim_file, mode='I') as writer:
  filenames = glob.glob('image*.png')
  filenames = sorted(filenames)
  for filename in filenames:
    image = imageio.imread(filename)
    writer.append_data(image)
  image = imageio.imread(filename)
  writer.append_data(image)


# %%
import tensorflow_docs.vis.embed as embed
embed.embed_file(anim_file)
