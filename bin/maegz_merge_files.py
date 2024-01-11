#! /home/siu/anaconda/envs/chem/bin/python

import argparse
import glob

parser = argparse.ArgumentParser(description = 'Merge name_prefix_*.maegz files into one file')
parser.add_argument('prefix', metavar='prefix', help='prefix like file_name_')
args = parser.parse_args()

prefix = args.prefix
file_list = glob.glob(prefix + '*.maegz')

import schrodinger
from schrodinger.structutils import sort

schrodinger.structutils.sort.merge_files(file_list, sort_criteria=[('s_m_title', sort.ASCENDING)], out_file_name='sort_tmp_MERGED.maegz', remove_file_list=False)