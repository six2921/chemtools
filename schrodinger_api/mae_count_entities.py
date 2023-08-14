import argparse
import os

parser = argparse.ArgumentParser(description = 'count cpds in maegz')
parser.add_argument('file', metavar='file', help='input .maegz files')
args = parser.parse_args()

file = args.file
name = file.split('.maegz')[0]

import gzip
from schrodinger import structure

# maegz 파일 읽기
suppl = structure.StructureReader(file)

# 물질 갯수 세기
count = 0
for st in suppl:
    count += 1

print(count)