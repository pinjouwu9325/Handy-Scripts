# This is a python script for systemical renaming files
# It is for the files downloaded from GDC(tcga-controlled files)

# Author: PJ Wu
# Last update: 2020-01-03

import os 
import pandas as pd


os.chdir("/NAS_lab/tcga/hnsc")
print("Current directory: " + os.getcwd())

meta = pd.read_table("gdc_sample_sheet.2019-12-26_paired.tsv")
print(meta.head())
# testID = "2d1ec8c9-ae37-456c-b1a3-595d72a5a2c1_gdc_realn_rehead.bam"
# meta["File Name" == testID, "Sample ID"]

# print(os.scandir("."))

for entry in os.scandir("."):
    if entry.is_dir():
        print("Directory: " + entry.name)
        for i in os.scandir(entry):
            if i.name.endswith(".bam"):
                 print("File to rename: " + i.name)
                 new = str(meta.loc[meta["File Name"] == i.name].iloc[0,6] + ".bam")
                 print("File's new name: " + new)
                 os.rename(i, new)
                 print("Renaming done")
