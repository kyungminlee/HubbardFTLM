#!/usr/bin/env python3

import sys, os
import h5py
import argparse
import pyarrow.feather
import glob
import json

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    parser.add_argument("scheduletype", type=str)
    parser.add_argument("--cutoff", type=int, default=100)
    args = parser.parse_args()
    lattice_str = args.lattice_str

    print("Reading sectors")
    sectors_df = pyarrow.feather.read_table(os.path.join(project_dir, "data", f"sectors-{lattice_str}.arrow")).to_pandas()
    schedule_df = pyarrow.feather.read_table(os.path.join(project_dir, "data", lattice_str, "schedule.arrow")).to_pandas()

    dense_filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigen-dense-results.hdf5")
    sparse_filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "sparsedata")

    if args.scheduletype == "dense":
        indices = schedule_df.idx[(schedule_df.type == "dense") & (schedule_df.type == "small")]
        print(find_missing_dense(dense_filepath, indices))
    elif args.scheduletype == "sparse":
        indices = schedule_df.idx[(schedule_df.type == "sparse")]
        print(find_missing_sparse(sparse_filepath, indices, args.cutoff))


def find_missing_dense(dense_filepath: str, indices: list):

    h5file = h5py.File(dense_filepath, "r")
    h5group = h5file["eigen-dense"]

    all_parameters = set()
    dense_group_lookup = {}
    for child_name in h5group:
        g = h5group[child_name]
        idx = g.attrs["idx"]
        t = g.attrs["hopping"]
        U = g.attrs["interaction"]
        dense_group_lookup[(idx, t, U)] = child_name
        all_parameters.add((t, U))
    h5file.close()

    missing_values = []
    for (t, U) in all_parameters:
        for idx in indices:
            if (idx, t, U) not in dense_group_lookup:
                missing_values.append((idx, t, U))
    return missing_values

def find_missing_sparse(sparse_filepath: str, indices: list, cutoff: int):
    all_parameters = set()
    sparse_group_count = {}
    for qed_path in glob.glob(os.path.join(sparse_filepath, "*.qed")):
        # parameter.json
        # {"idx":10000,"hopping":1.0,"interaction":20.0,"krylovdim":100,"seed":1}
        with open(os.path.join(qed_path, "parameter.json"), "r") as io:
            parameter = json.load(io)
        idx = parameter["idx"]
        t = parameter["hopping"]
        U = parameter["interaction"]
        count = len(glob.glob(os.path.join(qed_path, "eigenvalue_*.arrow")))
        sparse_group_count[(idx, t, U)] = count
        all_parameters.add((t, U))

    # print(sparse_group_count)

    missing_values = []
    for (t, U) in all_parameters:
        for idx in indices:
            if (idx, t, U) not in sparse_group_count or sparse_group_count[(idx, t, U)] < cutoff:
                missing_values.append((idx, t, U))
    return missing_values



if __name__=='__main__':
    main()



