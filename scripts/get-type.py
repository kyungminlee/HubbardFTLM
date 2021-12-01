#!/usr/bin/env python3
import os
import pandas as pd
import pyarrow.feather
import argparse
from pathlib import Path

project_dir = (Path(os.path.dirname(__file__)) / "..").resolve()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    parser.add_argument("type", type=str, choices=["small", "dense", "sparse"])
    args = parser.parse_args()

    filepath = project_dir / "data" / args.lattice_str / "schedule.arrow"
    df = pyarrow.feather.read_feather(str(filepath))

    idx_list = df.idx[df.type==args.type].tolist()
    print("\n".join(str(x) for x in idx_list))

if __name__=='__main__':
    main()