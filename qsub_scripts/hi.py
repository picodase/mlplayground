import os
import wget

os.mkdir("hi")

os.chdir("hi")

wget.download("https://files.rcsb.org/view/1QO1.pdb")