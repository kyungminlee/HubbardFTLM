#!/usr/bin/env python3
import os
import pandas as pd
import pyarrow as pa
import numpy as np

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("lattice_str", type=str)
parser.add_argument("nup", type=int)
parser.add_argument("ndn", type=int)
args = parser.parse_args()

df = pa.ipc.open_file(
    pa.OSFile(
        os.path.join(os.path.dirname(__file__), "..", "data", f"sectors-{args.lattice_str}.arrow"),
        "rb"
    ),
).read_pandas()
idx_list = df.idx[(df.nup == args.nup) & (df.ndn == args.ndn) & (df.idx == df.root_idx)].to_list()

print("\n".join(str(x) for x in idx_list))
