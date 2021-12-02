#!/usr/bin/env python3

import sys, os
import pandas as pd
import argparse
import pyarrow.feather
import glob
import json
import fastavro
from pathlib import Path

# project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
project_dir = (Path(os.path.dirname(__file__)) / "..").resolve()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    parser.add_argument("scheduletype", type=str)
    parser.add_argument("--cutoff", type=int, default=100)
    args = parser.parse_args()
    lattice_str = args.lattice_str

    print("Reading sectors")
    # sectors_df = pyarrow.feather.read_feather(str(project_dir / "data" / f"sectors-{lattice_str}.arrow"))
    schedule_df = pyarrow.feather.read_feather(str(project_dir / "data" / lattice_str /"schedule.arrow"))

    dense_filepath = project_dir / "data" / lattice_str / "eigen" / "dense-results.avro"
    sparse_filepath = project_dir / "data" / lattice_str / "eigen" / "sparsedata"

    if args.scheduletype == "dense":
        indices = schedule_df.idx[(schedule_df.type == "dense") & (schedule_df.type == "small")]
        with dense_filepath.open("rb") as io:
            reader = fastavro.reader(io)
            dense_table = pd.DataFrame([r for r in reader])
        print(find_missing_dense(dense_table, indices))
    elif args.scheduletype == "sparse":
        indices = schedule_df.idx[(schedule_df.type == "sparse")]
        print(find_missing_sparse(sparse_filepath, indices, args.cutoff))
    else:
        print(f"unsupported type {args.scheduletype}")


def find_missing_dense(table, indices: list):
    all_parameters = set()
    collection = set()
    for _, row in table.iterrows():
        t = row.hopping
        U = row.interaction
        idx = row.idx
        all_parameters.add((t, U))
        collection.add((t, U, idx))

    missing_values = []
    for (t, U) in all_parameters:
        for idx in indices:
            if (t, U, idx) not in collection:
                missing_values.append((idx, t, U))
    return missing_values

def find_missing_sparse(sparse_filepath: Path, indices: list, cutoff: int):
    all_parameters = set()
    sparse_group_count = {}
    for qed_path in sparse_filepath.glob("*.qed"):
        # parameter.json
        # {"idx":10000,"hopping":1.0,"interaction":20.0,"krylovdim":100,"seed":1}
        with (qed_path / "parameter.json").open("r") as io:
            parameter = json.load(io)
        t = parameter["hopping"]
        U = parameter["interaction"]
        idx = parameter["idx"]
        count = 0
        with (qed_path / "eigenvalues.jsonl").open("r") as io:
            for line in io:
                if line.strip():
                    count += 1
        sparse_group_count[(t, U, idx)] = count
        all_parameters.add((t, U))

    # print(sparse_group_count)

    missing_values = []
    for (t, U) in all_parameters:
        for idx in indices:
            if (t, U, idx) not in sparse_group_count:
                missing_values.append((t, U, idx, 0))
            elif sparse_group_count[t, U, idx] < cutoff:
                missing_values.append((t, U, idx, sparse_group_count[t, U, idx]))
                
    return missing_values



if __name__=='__main__':
    main()



