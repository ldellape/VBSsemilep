import os
import glob
import re

# Base search pattern
base_pattern = "/eos/user/l/ldellape/VBS/VBS_cards/production/jobs_output/ssWW_LL/job_*"
relative_path = "ssWW_LL_Run3Summer23wmLHEGS/root/Nanov12/NanoAODv12_Run3Summer23wmLHEGS.root"

# Find all job_* directories
job_dirs = sorted(glob.glob(base_pattern))

# Keep only those with job number > 200
filtered_dirs = []
for job_dir in job_dirs:
    match = re.search(r'job_(\d+)', job_dir)
    if match and int(match.group(1)) > 200:
        filtered_dirs.append(job_dir)

# Check for existing .root files
valid_paths = []
for job_dir in filtered_dirs:
    full_path = os.path.join(job_dir, relative_path)
    if os.path.isfile(full_path):
        valid_paths.append(full_path)

# Format list as a Python-style list of strings
formatted_list = ',\n'.join(f'"{p}"' for p in valid_paths)
output_content = f"[{formatted_list}]\n"

# Write to text file
output_file = "valid_paths.txt"
with open(output_file, "w") as f:
    f.write(output_content)

print(f"✅ Found {len(valid_paths)} valid ROOT files with job number > 200.")
print(f"List written to {output_file}")
