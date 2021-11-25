import os
import pandas as pd
import pyarrow as pa

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("lattice_str", type=str)
parser.add_argument("type", type=str, choices=["small", "dense", "sparse"])
args = parser.parse_args()

df = pa.ipc.open_file(
    pa.OSFile(
        os.path.join(os.path.dirname(__file__), "..", "data", args.lattice_str, "schedule.arrow"),
        "rb"
    ),
).read_pandas()
idx_list = df.idx[df.type==args.type].tolist()
print("\n".join(str(x) for x in idx_list))
