## Testing out psi4!
# %%
try:
    import os, sys
    sys.path.insert(1, os.path.abspath('/scratch/psilocaluser/conda-builds/psi4-docs-multiout_1557966099526/work/build/stage//scratch/psilocaluser/conda-builds/psi4-docs-multiout_1557966099526/_h_env_placehold_placehold_place/lib/python3.6/site-packages'))
except ImportError:
    pass

import psi4
# %% 
psi4.core.set_output_file('output.dat', False)
# %%
#! Sample HF/cc-pVDZ H2O Computation

psi4.set_memory('4000 MB')

h2o = psi4.geometry("""
O
H 1 0.96
H 1 0.96 2 104.5
""")

psi4.energy('scf/cc-pvdz')
# %%


#! Sample UHF/6-31G** CH2 Computation

R = 1.075
A = 133.93

ch2 = psi4.geometry("""
0 3
C
H 1 {0}
H 1 {0} 2 {1}
""".format(R, A)
)

psi4.set_options({'reference': 'uhf'})
psi4.energy('scf/6-31g**')


# %%
# Example SAPT computation for ethene*ethyne (*i.e.*, ethylene*acetylene).
# Test case 16 from S22 Database

dimer = psi4.geometry("""
0 1
C   0.000000  -0.667578  -2.124659
C   0.000000   0.667578  -2.124659
H   0.923621  -1.232253  -2.126185
H  -0.923621  -1.232253  -2.126185
H  -0.923621   1.232253  -2.126185
H   0.923621   1.232253  -2.126185
--
0 1
C   0.000000   0.000000   2.900503
C   0.000000   0.000000   1.693240
H   0.000000   0.000000   0.627352
H   0.000000   0.000000   3.963929
units angstrom
""")
psi4.set_options({'scf_type': 'df',
                  'freeze_core': 'true'})

psi4.energy('sapt0/jun-cc-pvdz', molecule=dimer)