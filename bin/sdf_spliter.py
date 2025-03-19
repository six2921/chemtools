import argparse
from sdf_modules import split_sdf_file

# Set up argument parsing
parser = argparse.ArgumentParser(description='Split your SDF file')
parser.add_argument('sdf', help='Your SDF file')
parser.add_argument('bunch', type=int, help='How many molecules per file')
args = parser.parse_args()

# Split the SDF file
generated_files = split_sdf_file(args.sdf, chunk_size=args.bunch)

# Print the generated files
for file in generated_files:
    print(f'{file} is generated')

print("All files are generated.")
