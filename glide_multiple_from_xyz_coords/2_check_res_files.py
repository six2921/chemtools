import argparse
import os
import pandas as pd
import glob
import subprocess

parser = argparse.ArgumentParser(description='Check number of generated files in your dir')
parser.add_argument('dir_path', type=str, help='The dir is what you ran generate_grids_from_coordinates.py or run_glide_multiple_grids.py')

# 경로 변경
os.chdir(parser.parse_args().dir_path)

# CSV 파일에서 데이터 읽기
df = pd.read_csv('coordinates.csv')
print(f'Number of coordinates: {len(df)}')

# Check number of .zip files
grids = glob.glob('grid_*.zip')
print(f'Number of .zip files: {len(grids)}')

# Check number of .zip files
poses = glob.glob('glide_*.sdf')
print(f'Number of glide_*.sdf files: {len(poses)}')