#! /home/siu/anaconda/bin/python

import argparse
import wget

parser = argparse.ArgumentParser(description = 'download pdb file from your code (1crn)')
parser.add_argument('code', metavar='code', help='pdb code')
args = parser.parse_args()

code = args.code

wget.download(f"https://files.rcsb.org/download/{code}.pdb")

print('')
print(f"You got {code}.pdb")