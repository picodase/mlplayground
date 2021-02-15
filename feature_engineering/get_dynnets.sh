# setup conda environment
conda activate getcontacts

# set variables
expt_name="mdsim_07_10_10_13_49_07"    # experiment name
filepath="data/JLN_1_21_1A/"    # filepath
suffix=".tsv"
filename=${filepath}${expt_name}${suffix}
cores=2

# run getcontacts.py scripts

################################################### <get_dyn> ##############################################

# salt bridge
get_dynamic_contacts.py --topology ${filepath}md_0_1.gro --trajectory ${filepath}md_0_1.xtc --output dyncontacts_sb_${expt_name}.tsv --cores ${cores} --itypes sb

# pi-cation
get_dynamic_contacts.py --topology ${filepath}md_0_1.gro --trajectory ${filepath}md_0_1.xtc --output dyncontacts_pc_${expt_name}.tsv --cores ${cores} --itypes pc

# pi stacking
get_dynamic_contacts.py --topology ${filepath}md_0_1.gro --trajectory ${filepath}md_0_1.xtc --output dyncontacts_ps_${expt_name}.tsv --cores ${cores} --itypes ps

# t-stacking
get_dynamic_contacts.py --topology ${filepath}md_0_1.gro --trajectory ${filepath}md_0_1.xtc --output dyncontacts_ts_${expt_name}.tsv --cores ${cores} --itypes ts

# hydrophobic
get_dynamic_contacts.py --topology ${filepath}md_0_1.gro --trajectory ${filepath}md_0_1.xtc --output dyncontacts_hp_${expt_name}.tsv --cores ${cores} --itypes hp

# van der Waals
get_dynamic_contacts.py --topology ${filepath}md_0_1.gro --trajectory ${filepath}md_0_1.xtc --output dyncontacts_vdw_${expt_name}.tsv --cores ${cores} --itypes vdw

#   --solv IP3
#   --sele "(chain A and resid 100 to 160)" \
#   --sele2 "resname EJ4" \
#   --itypes sb hb

################################################### </get_dyn> #############################################

# get_contact_frequencies.py
#get_contact_frequencies.py --input_files dyncontacts_* --output_file contfreqs_${filename}

