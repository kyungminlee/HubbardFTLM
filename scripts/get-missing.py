import sys, os
import h5py
import argparse
import pyarrow.feather

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    args = parser.parse_args()

    print("Reading sectors")
    sectors_df = pyarrow.feather.read_table(os.path.join(project_dir, "data", f"sectors-{lattice_str}.arrow")).to_pandas()

    lattice_str = args.lattice_str
    dense_filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigen-dense-results.hdf5")

    h5file = h5py.File(filepath, "r")
    h5group = h5file["eigen-dense"]

    all_parameters = set()
    group_lookup = {}
    for child_name in h5group:
        g = h5group[child_name]
        idx = g.attrs["idx"]
        t = g.attrs["hopping"]
        U = g.attrs["interaction"]
        group_lookup[(idx, t, U)] = child_name
        all_parameters.add((t, U))


