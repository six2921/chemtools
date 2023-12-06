import argparse
import os
import sys
import pandas as pd
import glob
import subprocess

parser = argparse.ArgumentParser(description='Generate grids from coordinates in csv file')
parser.add_argument('dir_path', type=str, help='coordinates.csv and receptor.mae must be in this directory')

# 슈뢰딩거 파일 경로 읽기
with open('path_of_schrodinger.txt', 'r') as file:
    SRDG = file.readline().strip('PATH: ').strip()
    HOST = file.readline().strip('HOST: ').strip()

# 경로 변경
os.chdir(parser.parse_args().dir_path)

if not os.path.exists('coordinates.csv') or not os.path.exists('receptor.mae'):
    print('Error: 경로 내에 coordinates.csv 또는 receptor.mae 파일이 존재하지 않습니다.')
    sys.exit()

# CSV 파일에서 데이터 읽기
df = pd.read_csv('coordinates.csv')

# 각 행에 대해 파일 생성
for i, row in df.iterrows():
    with open(f'grid_{i}.in', 'w') as f:
        f.write(f'GRID_CENTER {row["x"]}, {row["y"]}, {row["z"]}\n')
        f.write('RECEP_FILE receptor.mae\n')
        f.write(f'GRIDFILE grid_{i}.zip\n')

configs = glob.glob('grid_*.in')

for conf in configs:
    subprocess.run([f'{SRDG}/glide', conf, '-LOCAL', '-HOST', f'{HOST}'])