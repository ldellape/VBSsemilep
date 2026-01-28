import json
import ast
import glob

INPUT_PATTERN = "./*/2023/*.json"
OUTPUT_FILE = "das_names.py"

das_set = set()

for json_file in glob.glob(INPUT_PATTERN):
    with open(json_file, "r") as f:
        data = json.load(f)

    for sample_name, sample_data in data.items():
        metadata = sample_data.get("metadata", {})
        das_names_str = metadata.get("das_names")

        if das_names_str:
            try:
                das_list = ast.literal_eval(das_names_str)
                for das in das_list:
                    das_set.add(das)
            except Exception as e:
                print(f"⚠ Errore in {json_file} ({sample_name}): {e}")

# scrittura come lista Python
with open(OUTPUT_FILE, "w") as out:
    out.write("das_names = [\n")
    for das in sorted(das_set):
        out.write(f"    '{das}',\n")
    out.write("]\n")

print(f"✔ Creata lista Python con {len(das_set)} entries")
print(f"📄 Salvata in {OUTPUT_FILE}")
