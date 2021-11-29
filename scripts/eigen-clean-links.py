import h5py
import json
import pyarrow.feather

import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("h5file", type=str)
    args = parser.parse_args()

    h5f = h5py.File(args.h5file, "a")
    h5g = h5f["eigen-dense"]
    for child_name in h5g:
        d = h5g[child_name]
        l = d.get("eigenvalue", getlink=True)
        if isinstance(l, h5py.SoftLink):
            del h5g[child_name]

    h5f.close()

if __name__=='__main__':
    main()
