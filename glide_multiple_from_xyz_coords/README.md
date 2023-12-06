# Requirement
1. Make folder and put below files there.
2. receptor.mae: Prepared in Schroidnger 
3. coordinates.csv: x,y,z coordinates which is centroid you want to generate grid. Multiple rows are avaiable. x,y,z columns are must be exist. 
4. ligands_prep.sdf: Prepared in Schrodinger

# How to run
1. vi path_of_schroinger.txt.   # Modify the installation path of schrodinger and name of host for -HOST parameter. 
2. python 1_generate_grids_from_coordinates.py <your-folder>   # Read x,y,z from csv file and generate in for grid generation and run it. 
3. python 2_check_res_files.py   # You can check number of generated result files from sciprt 1 and 3.
4. python 3_run_glide_multiple_grids.py   # List up .zip(grid) files and generate in file for glide docking and run it. 
