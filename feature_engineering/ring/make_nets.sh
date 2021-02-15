# export VICTOR_ROOT
export VICTOR_ROOT=/home/jnorth/Documents/programs/ring_exe/dist/

# run loop
for fn in */*.pdb; do "Ring -i ${fn} -E ${fn/.pdb/}_edges.txt"; done
