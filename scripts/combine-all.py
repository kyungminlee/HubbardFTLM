import sys, os
import h5py
import numpy as np
import pandas as pd
import pyarrow.feather
import argparse

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
lattice_str = "square-(4,0)x(0,4)"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    args = parser.parse_args()
    print("Reading sectors")
    sectors_df = pyarrow.feather.read_table(os.path.join(project_dir, "data", f"sectors-{lattice_str}.arrow")).to_pandas()

    print("Reading schedule")
    schedule_df = pyarrow.feather.read_table(os.path.join(project_dir, "data", lattice_str, f"schedule.arrow")).to_pandas()

    df = pd.merge(sectors_df, schedule_df, how="left")

    out_filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigenvalues.hdf5")
    out_h5file = h5py.File(out_filepath, "w")

    dense_filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigen-dense-results.hdf5")
    dense_h5file = h5py.File(dense_filepath, "r")
    complete_dense(df, out_h5file, dense_h5file)
    dense_h5file.close()

    sparse_filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigen-sparse-results.hdf5")
    sparse_h5file = h5py.File(sparse_filepath, "r")
    complete_sparse(df, out_h5file, sparse_h5file)
    sparse_h5file.close()

    out_h5file.flush()

    check_complete(df, out_h5file)

    out_h5file.close()


def complete_dense(df: pd.DataFrame, out_h5file, dense_h5file):
    in_h5 = dense_h5file["eigen-dense"]
    out_h5 = out_h5file

    print("Reading result file")
    all_parameters = set()
    for (child_name, in_group) in in_h5.items():
        scheduletype = in_group.attrs["type"]
        if isinstance(scheduletype, str):
            pass
        elif isinstance(scheduletype, np.bytes_):
            scheduletype = scheduletype.decode()
        else:
            assert(False)
        if scheduletype not in ["small", "dense"]: continue
        t = in_group.attrs["hopping"]
        U = in_group.attrs["interaction"]
        all_parameters.add((t, U))

        print(f"dense copying {child_name}")
        out_group = out_h5.create_group(child_name)
        for k, v in in_group.attrs.items():
            out_group.attrs[k] = v
        out_group["eigenvalue"] = np.array(in_group["eigenvalue"])

    #for irow, (idx, nup, ndn, tii, pii, pic, dim, root_idx, matrixtype) in sectors_df.iterrows():
    for irow, row in df.iterrows():
        assert(irow+1 == row.idx)
        if row.idx == row.root_idx: continue
        if row.dim == 0: continue

        assert(pd.isna(row.type))

        root_row = df.loc[row.root_idx-1]
        assert(root_row.dim != 0)

        if root_row.type not in ["small", "dense"]: continue

        for (t, U) in all_parameters:
            src_name = f"hopping={t}_idx={row.root_idx}_interaction={U}"
            dst_name = f"hopping={t}_idx={row.idx}_interaction={U}"
            if src_name not in in_h5:
                print(f"missing: {src_name}")
                continue
            #print(f"linking: idx={row.idx}, t={t}, U={U}, root_idx={row.root_idx}")
            #print(f"creating link: {dst_name}")
            if dst_name in out_h5:
                print(f"{dst_name} exists in output")
                print("old data:")
                for k, v in out_h5[dst_name].attrs.items():
                    print(f"  {k}: {v}")
                print("new data:")
                print(f"  idx: {row.idx}")
                print(f"  hopping: {t}")
                print(f"  interaction: {U}")
                print(f"  type: dense-link")
                print(f"  root_idx: {row.root_idx}")
            dst_group = out_h5.create_group(dst_name)
            dst_group.attrs["idx"] = row.idx
            dst_group.attrs["hopping"] = t
            dst_group.attrs["interaction"] = U
            dst_group.attrs["root_idx"] = row.root_idx
            dst_group.attrs["type"] = "dense-link"
            dst_group["eigenvalue"] = h5py.SoftLink(f"/{src_name}/eigenvalue")
        

def complete_sparse(df: pd.DataFrame, out_h5file, in_h5file):
    in_h5 = in_h5file["eigen-sparse"]
    out_h5 = out_h5file

    print("Reading sparse result file")
    all_parameters = set()
    for (child_name, in_group) in in_h5.items():
        scheduletype = in_group.attrs["type"]
        if isinstance(scheduletype, str):
            pass
        elif isinstance(scheduletype, np.bytes_):
            scheduletype = scheduletype.decode()
        else:
            assert(False)
        if scheduletype != "sparse": continue
        t = in_group.attrs["hopping"]
        U = in_group.attrs["interaction"]
        all_parameters.add((t, U))

        print(f"sparse copying {child_name}")
        out_group = out_h5.create_group(child_name)
        for k, v in in_group.attrs.items():
            out_group.attrs[k] = v
        out_group["eigenvalue"] = np.array(in_group["eigenvalue"])

    #for irow, (idx, nup, ndn, tii, pii, pic, dim, root_idx, matrixtype) in sectors_df.iterrows():
    for irow, row in df.iterrows():
        if row.idx == row.root_idx: continue
        if row.dim == 0: continue

        assert(pd.isna(row.type))

        root_row = df.loc[row.root_idx-1]
        assert(root_row.dim != 0)

        if root_row.type != "sparse": continue

        for (t, U) in all_parameters:
            src_name = f"hopping={t}_idx={row.root_idx}_interaction={U}"
            dst_name = f"hopping={t}_idx={row.idx}_interaction={U}"
            if src_name not in in_h5:
                print(f"missing: {src_name}")
                continue
            #print(f"linking: idx={row.idx}, t={t}, U={U}, root_idx={row.root_idx}")
            #print(f"creating link: {dst_name}")
            if dst_name in out_h5:
                print(f"{dst_name} exists in output")
                print("old data:")
                for k, v in out_h5[dst_name].attrs.items():
                    print(f"  {k}: {v}")
                print("new data:")
                print(f"  idx: {row.idx}")
                print(f"  hopping: {t}")
                print(f"  interaction: {U}")
                print(f"  type: sparse-link")
                print(f"  root_idx: {row.root_idx}")
                continue

            dst_group = out_h5.create_group(dst_name)
            dst_group.attrs["idx"] = row.idx
            dst_group.attrs["hopping"] = t
            dst_group.attrs["interaction"] = U
            dst_group.attrs["root_idx"] = row.root_idx
            dst_group.attrs["type"] = "sparse-link"
            dst_group["eigenvalue"] = h5py.SoftLink(f"/{src_name}/eigenvalue")


def check_complete(df, h5file):
    all_parameters = set()
    for (child_name, in_group) in h5file.items():
        t = in_group.attrs["hopping"]
        U = in_group.attrs["interaction"]
        all_parameters.add((t, U))

    missing_count = 0
    exist_count = 0
    for _, row in df[df.dim != 0].iterrows():
        for t, U in all_parameters:
            name = f"hopping={t}_idx={row.idx}_interaction={U}"
            if name not in h5file:
                print(f"missing: {name}")
                print(row)
                missing_count += 1
            else:
                exist_count += 1
    print(f"missing: {missing_count}")
    print(f"exist: {exist_count}")


if __name__=='__main__':
    main()
