import awkward as ak
import glob

# Path to your parquet files (adjust as needed)
files = glob.glob("/eos/user/l/ldellape/VBS/parquet/ssWW_TT_mg5_madspin/Single*_AK8/*.parquet")

total_entries = 0

for f in files:
    arr = ak.from_parquet(f)
    n_entries = len(arr)
    print(f"{f}: {n_entries} entries")
    total_entries += n_entries

print(f"\nTotal entries across all files: {total_entries}")
