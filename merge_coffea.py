#!/usr/bin/env python3

import coffea.util as util

output_file = "2023_preBPix_isoMu27.coffea"

input_files = [
#    "topCR_VVV/output_all.coffea",
#    "topCR_Wjets/output_all.coffea",
  #  "topCR_data_2023preBPix/output_all.coffea",
#    "topCR_data_2/output_all.coffea",
    "data2023C/output_all.coffea",
    "2023C_top/output_all.coffea",
]

def merge_objects(a, b):
    """
    Recursively merge coffea outputs.
    - dict -> recurse
    - set  -> union
    - everything else -> use +
    """
    if isinstance(a, dict) and isinstance(b, dict):
        out = a.copy()
        for k, v in b.items():
            if k in out:
                out[k] = merge_objects(out[k], v)
            else:
                out[k] = v
        return out

    elif isinstance(a, set) and isinstance(b, set):
        return a | b

    else:
        return a + b


def merge_coffea(output_file, input_files):
    merged = None

    for fname in input_files:
        print(f"Loading {fname}")
        acc = util.load(fname)

        if merged is None:
            merged = acc
        else:
            merged = merge_objects(merged, acc)

    print(f"Saving merged output to {output_file}")
    util.save(merged, output_file)


if __name__ == "__main__":
    merge_coffea(output_file, input_files)
