#! /home/siu/anaconda/envs/chem/bin/python

import argparse
import os

parser = argparse.ArgumentParser(description = 'convert maegz to sdf and csv')
parser.add_argument('file', metavar='file', help='input .maegz files')
args = parser.parse_args()

file = args.file
name = file.split('.maegz')[0]

print(name)

if '_pv' in name:
    os.system(f"$SCHRODINGER/utilities/glide_sort {name}.maegz -o {name}.mae -norecep -best_by_title")
else:
    os.system(f"$SCHRODINGER/utilities/structconvert {name}.maegz {name}.mae")

os.system(f"$SCHRODINGER/utilities/structconvert {name}.mae {name}.sdf")
os.system(f"$SCHRODINGER/utilities/structconvert {name}.sdf {name}.csv")
