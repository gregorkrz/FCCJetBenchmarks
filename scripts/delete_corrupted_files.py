
import os

filelist = "/fs/ddn/sdf/group/atlas/d/gregork/fastsim/jetbenchmarks/IDEA_20260120/corrupted_2601_2.txt"

with open(filelist, "r") as f:
    corrupted_files = f.readlines() # each line is a file path + \t + reason
    corrupted_files = [line.split("\t")[0] for line in corrupted_files]

print(corrupted_files)

for file in corrupted_files:
    if os.path.exists(file):
        print(f"Deleting corrupted file: {file}")
        os.remove(file)
    else:
        print(f"File not found, skipping deletion: {file}")