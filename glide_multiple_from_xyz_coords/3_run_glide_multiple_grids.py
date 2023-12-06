import argparse
import os
import sys
import glob
import subprocess

parser = argparse.ArgumentParser(description='run glide for multiple grids in your dir')
parser.add_argument('dir_path', type=str, help='grid_*.in files and ligands_prep.sdf must be in this directory')

# 슈뢰딩거 파일 경로 읽기
with open('path_of_schrodinger.txt', 'r') as file:
    SRDG = file.readline().strip('PATH: ').strip()
    HOST = file.readline().strip('HOST: ').strip()

# 경로 변경
os.chdir(parser.parse_args().dir_path)

if not os.path.exists('ligands_prep.sdf'):
    print('Error: 경로 내에 ligands_prep.sdf 파일이 존재하지 않습니다.')
    sys.exit()

grids = glob.glob('grid_*.zip')

# 각 행에 대해 파일 생성
num = 0
for gr in grids:
    with open(f'glide_{num}.in', 'w') as f:
        f.write(f'GRIDFILE {gr}\n')
        f.write('LIGANDFILE ligands_prep.sdf\n')
        f.write(f'PRECISION SP\n')
        f.write(f'POSE_OUTTYPE ligandlib_sd\n')
        f.write(f'COMPRESS_POSES False\n')
    num += 1

configs = glob.glob('glide_*.in')

for conf in configs:
    subprocess.run([f'{SRDG}/glide', conf, '-LOCAL', '-HOST', f'{HOST}'])