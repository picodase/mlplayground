# ring_analyzer.py

# %% 
# Imports
import pandas as pd
import numpy as np
import igraph
# %%
# write a function to take in data from ring-format txt-files
import os

def importRingE(fpref:str):
    edges = pd.read_csv(fpref+"_edges.txt", sep="\t")    # read in the data to a pandas dataframe (TSV)
    #nodes = pd.read_csv(fpref+"_nodes.txt", sep="\t")
    return edges

# %%
# import the network of interest

path = "/home/jnorth/Nextcloud/extracurriculars/publishing/research/hsu/experiments/JLN_1_21b"
filename = "5x8x"
fpref = path+"/"+filename

df = importRingE(fpref)
# %%
# extract NodeId1, interaction type, and NodeId2 to a new df
df_pruned = df[["NodeId1","Interaction","NodeId2"]]

df_pruned

# %% 
# partitioning scheme for each network type
hbond = {"HBOND:MC_MC", "HBOND:SC_MC", "HBOND:MC_SC", "HBOND:SC_SC"}
vdw = {"VDW:SC_SC", "VDW:MC_SC", "VDW:SC_MC", "VDW:LIG_SC"}
lig = {"IAC:LIG_SC","IAC:LIG_MC","VDW:LIG_SC"}
ππ = {"PIPISTACK:SC_SC"}
# %%
# get protein structures
# Libraries
import pypdb       # For searching and downloading PDB files
import datetime     # for organizing file output
import Bio      # Bio for downloading PDBs, printing protein lengths in residues
from Bio.PDB import PDBList
from Bio.PDB import PDBParser

# %%
query = input("Supply a query (term, accession number, etc.) :")

query = "nuclear receptor"
# Tag the time right when the query is entered

now = datetime.datetime.now()
def now_dir_ts():
    now_ts = str(now.year)+"_"+str(now.month)+"_"+str(now.day)+"_"+str(now.hour)+"_"+str(now.minute)+"_"+str(now.second)
    return now_ts

now = now_dir_ts()
PDB_dl_dir = "ds_"+now

search_dict = pypdb.Query(query)     # create a dictionary containing search information
found = search_dict.search(search_dict)[:500]      # create a list of these PDBs by searching RCSB

# create a list with the information and the metadata
metadata = []

for proteins in found:  # for items in # for the items in the list,
    metadata.append(pypdb.describe_pdb(proteins))  # append the dictionary 

# Save the metadata list as a CSV file
dfm = pd.DataFrame(metadata) # convert to a Pandas DF
dfm.to_csv('metadata_'+now+'.csv')      # save as a CSV file

# %%
parser = PDBParser()       # create a parser
pdbs = list()
pdbl = PDBList()

# Download all PDB structures in the previous list if they aren't there
for id in found:
    pdbl.retrieve_pdb_file(pdb_code=id, file_format='pdb', pdir=PDB_dl_dir)   # Retrieve in PDB format, put in directory 'PDB'

# Finished, print "Downloading ... finished!"
print('\n#############~DOWNLOADING COMPLETE~#############\n')
# %%
# convert pdb*.ent to *.pdb
for file in os.scandir(PDB_dl_dir):
    if (file.path.endswith(".ent") and file.is_file()):
        newfn = file.name.replace("pdb","").replace(".ent",".pdb")
        os.rename(file, PDB_dl_dir+"/"+newfn)

# %%
# run ring2.0 on each pdb file
!for fn in */*.pdb; do "./dist/bin/Ring -i ${fn} -E ${fn/.pdb/_edges.txt}"; done

# !for fn in */*.pdb; do "bash /home/jnorth/Documents/programs/ring_exe/dist/bin/Ring -i ${fn/.pdb/}.pdb -E ${fn/.pdb/_edges.txt}"; done

#%%
import os

rin_dfs = []
fileNames = []
for file in os.scandir(path):
    if (file.path.endswith("_edges.txt") and file.is_file()):
      #print(file.name.replace(".sif",""))
      fileNames.append(file.name.split(".")[0])
      df = pd.read_csv(file, sep='\t')[["NodeId1","Interaction","NodeId2"]].rename(columns = {'NodeId1':'resA', 'Interaction':'interacType', 'NodeId2':'resB'}, inplace = True) 
      rin_dfs.append(df)

# %%
# create lists to hold each set of values for the specified interaction types
interac_cts = []
edgelists = []
graphs = []
clusters = []
cluster_sizes = []
diams = []
#d_norm = []
thr_motifs = []
four_motifs = []

for df in rin_dfs:
  # interaction counts
  interac_cts.append(df.interacType.value_counts())

  # interaction types
  i = {}
  i["tot"] = df[['resA','resB']].values.tolist()
  i["hbd"] = df[df.interacType.isin(hbond)][['resA','resB']].values.tolist()
  i["vdw"] = df[df.interacType.isin(vdw)][['resA','resB']].values.tolist()
  i["lig"] = df[df.interacType.isin(lig)][['resA','resB']].values.tolist()
  i["ππ"] = df[df.interacType.isin(ππ)][['resA','resB']].values.tolist()
  edgelists.append(i)

  # graphs
  g = {}
  g["tot"] = igraph.Graph.TupleList(i["tot"])
  g["hbd"] = igraph.Graph.TupleList(i["hbd"], directed=True)
  g["vdw"] = igraph.Graph.TupleList(i["vdw"])
  g["lig"] = igraph.Graph.TupleList(i["lig"])
  g["ππ"] = igraph.Graph.TupleList(i["ππ"])
  graphs.append(g)

  # clusters
  c = {}
  c["tot"] = g["tot"].clusters()
  c["hbd"] = g["hbd"].clusters()    
  c["vdw"] = g["vdw"].clusters()
  c["lig"] = g["lig"].clusters()
  c["ππ"] = g["ππ"].clusters()
  clusters.append(c)

  # cluster sizes
  c_s = {}
  c_s["tot"] = c["tot"].sizes()
  c_s["hbd"] = c["hbd"].sizes()
  c_s["vdw"] = c["vdw"].sizes()
  c_s["lig"] = c["lig"].sizes()
  c_s["ππ"] = c["ππ"].sizes()
  cluster_sizes.append(c_s)

  # diameters
  d = {}
  d["tot"] = g["tot"].diameter()
  d["hbd"] = g["hbd"].diameter()
  d["vdw"] = g["vdw"].diameter()
  d["lig"] = g["lig"].diameter()
  d["ππ"] = g["ππ"].diameter()
  diams.append(d)

  # diameter normalized to number of vertices in the network
  #d_n = {}

  #d_n["hbd"] = (d["hbd"]/len(g["hbd"].vs) if (len(g["hbd"].vs) != 0) else np.nan)
  #d_n["vdW"] = (d["vdw"]/len(g["vdw"].vs) if (len(g["vdw"].vs) != 0) else np.nan)
  #d_n["lig"] = (d["lig"]/len(g["lig"].vs) if (len(g["lig"].vs) != 0) else np.nan)
  #d_n["ππ"] = (d["ππ"]/len(g["ππ"].vs) if (len(g["ππ"].vs) != 0) else np.nan)
  #d_norm.append(d)

  # motifs
  # calculate motifs for each graph

  # three-motifs
  t_m = {}
  t_m["tot"] = g["tot"].motifs_randesu()
  t_m["hbd"] = g["hbd"].motifs_randesu()
  t_m["vdw"] = g["vdw"].motifs_randesu()
  t_m["lig"] = g["lig"].motifs_randesu()
  t_m["ππ"] = g["ππ"].motifs_randesu()
  thr_motifs.append(t_m)

  # four-motifs
  f_m = {}
  f_m["tot"] = g["tot"].motifs_randesu(size=4)
  f_m["hbd"] = g["hbd"].motifs_randesu(size=4)
  f_m["vdw"] = g["vdw"].motifs_randesu(size=4)
  f_m["lig"] = g["lig"].motifs_randesu(size=4)
  f_m["ππ"] = g["ππ"].motifs_randesu(size=4)
  four_motifs.append(f_m)