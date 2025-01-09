import os
import numpy as np
import argparse

def calculate_centroid(file_path):
    coordinates = []

    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) >= 3:
                x, y, z = map(float, parts[:3])
                coordinates.append([x, y, z])

    coordinates_array = np.array(coordinates)
    centroid = np.mean(coordinates_array, axis=0)
    return centroid

def process_xyzrg_files(directory):
    for filename in os.listdir(directory):
        if filename.endswith('.xyzrg'):
            file_path = os.path.join(directory, filename)
            centroid = calculate_centroid(file_path)
            
            # Extract the identifier from the filename, assuming it's in the form "p0", "p1", etc.
            identifier = filename.split('_')[-1].split('.')[0]

            # Print the identifier and centroid in the desired format
            print(f"{identifier}")
            print(f"center_x = {centroid[0]:.3f}")
            print(f"center_y = {centroid[1]:.3f}")
            print(f"center_z = {centroid[2]:.3f}\n")

# argparse로 폴더 경로를 입력받음
parser = argparse.ArgumentParser(description="Process .xyzrg files and calculate centroids.")
parser.add_argument("directory", type=str, help="Directory containing .xyzrg files")
args = parser.parse_args()

# 지정된 폴더에서 .xyzrg 파일 처리
process_xyzrg_files(args.directory)

