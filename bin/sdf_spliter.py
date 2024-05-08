#! /home/siu/anaconda/bin/python

import os
import argparse

parser = argparse.ArgumentParser(description = 'Organize your sdf file')
parser.add_argument('sdf', help='your sdf file') # 필요한 인수를 추가
parser.add_argument('bunch', help='how many molecules per a file') # 필요한 인수를 추가
args = parser.parse_args()

sdf = args.sdf
bunch = int(args.bunch)

fn = os.path.splitext(sdf)[0]
counter = 0
file_counter = 0
output_file = open(f'{fn}_{file_counter}.sdf', 'w')

with open(sdf, 'r') as f:
    for line in f:
        output_file.write(line)
        if line.strip() == "$$$$":
            counter += 1
            if counter % bunch == 0:
                output_file.close()
                file_counter += 1
                output_file = open(f'{fn}_{file_counter}.sdf', 'w')
                print(f'{fn}_{file_counter}.sdf is generated')

output_file.close()

print("File_N.sdf are generated")