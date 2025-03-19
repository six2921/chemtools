import argparse

parser = argparse.ArgumentParser(description = 'split into N files from maegz (max=100T cps)')
parser.add_argument('file', metavar='file', help='input .maegz files')
args = parser.parse_args()

file = args.file
name = file.split('.maegz')[0]

import schrodinger
from schrodinger.structutils import sort

schrodinger.structutils.sort.split_file(file, max_count=100000, dir=None)