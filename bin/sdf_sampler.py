#! /home/siu/anaconda3/envs/rdkit/bin/python
import argparse
import os

# 인자 파싱
parser = argparse.ArgumentParser()
parser.add_argument('filename', help='input file name')
parser.add_argument('-s', '--size', type=int, default=10, help='number of $$$$ to find')
args = parser.parse_args()

# 파일 읽기
with open(args.filename, 'r') as f:
    new_lines = []
    count = 0
    for line in f:
        if line.startswith('$$$$'):
            count += 1
            if count == args.size:
                new_lines.append(line)
                break
        new_lines.append(line)

# 파일 저장
filename, ext = os.path.splitext(args.filename)
new_filename = f'{filename}_sample{ext}'
with open(new_filename, 'w') as f:
    f.writelines(new_lines)

# 결과 출력
print(f'File saved to {new_filename}')