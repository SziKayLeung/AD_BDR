#!/usr/bin/python
# Szi Kay Leung
# 21/06/2022: Download RNA-Seq data from ROSMAP dataset
# python 1_download_dataset.py <output_dir>

# Modules and arguments
import synapseclient
import sys
import pandas as pd

syn = synapseclient.Synapse()
output_dir = sys.argv[1] 

# login to Synapse
print("Logging into Synapse")
syn.login("e.pishva@exeter.ac.uk", "Ellecuylgaard105*")

# synpase ID ordered from merged_metadata.R
# RNASeq dataset from https://www.synapse.org/#!Synapse:syn8612203
synapseid = pd.read_table("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/Scripts/AD_BDR/2_Differential_Analysis/Mayo/Metadata/mayo_ordered_synapseid.txt")

# Transfer all files using synapseID
for counter, id in enumerate(synapseid["SynID"].values):
  print("##################", counter,"####################")
  print("Transferring", id)
  entity = syn.get(id, downloadLocation=output_dir)
