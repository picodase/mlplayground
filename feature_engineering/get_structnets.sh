# setup conda environment
#eval "$(command conda 'shell.bash' 'hook' 2> /dev/null)"
#conda activate getcontacts

# set variables
expt_name="mdsim_07_10_10_13_49_07"    # experiment name
filepath="data/JLN_1_21_1A/"    # filepath
suffix=".tsv"
filename=${filepath}${expt_name}${suffix}
cores=2

# run getcontacts.py scripts

################################################### <get_stat> ##############################################

#for fn in */*/*/*.pdb

for fn in */*.pdb
do
    # all
    # select only protein-ligand (incl. everything that's not a biomolecule) interactions
    #get_static_contacts.py --structure ${fn} --output ${fn/.pdb/}_statcont_all_protlig.tsv --itypes all --sele2 'not (protein or nucleic or solv or lipid)'

    # select every atom in the network
    #get_static_contacts.py --structure ${fn} --output ${fn/.pdb/}_statcont_all.tsv --itypes all --sele 'all' --sele2 'all'

    # cluster version
    /home/other/northj/packages/getcontacts/get_static_contacts.py --structure ${fn} --output ${fn/.pdb/}_statcont_all.tsv --itypes all --sele 'all' --sele2 'all'

    #   --solv HOH 
    #   --sele "chain A"
    #   --sele2 "chain B"
    #   --itypes sb hb
done
