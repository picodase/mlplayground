#!/bin/sh

# rename_traj_as_foldername.sh

direct=${PWD}


# make a directory in the base experiment folder, to hold the processed, named trajectory files by experiment characteristics
mkdir to_mdanaly

for dir in ${PWD}/*/    # for each subfolder,
do 
    cd "$dir"      # move into that folder

    #direct=$(basename ${PWD})
    #echo ${direct}
    filename=$(basename ${PWD})

    mv md_0_1.xtc ${filename}.xtc    # rename the trajectory to the name of the folder + .xtc
    mv md_0_1.gro ${filename}.gro    # also rename the topology file

    cp ${filename}.xtc ${filename}.gro ../to_mdanaly  # copy the file to an upper folder for analysis

    cd ..   # cd out of that folder to the base
done