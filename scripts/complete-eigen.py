import sys, os
import h5py
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

    complete_dense(df, args.lattice_str)
    complete_sparse(df, args.lattice_str)


def complete_dense(df: pd.DataFrame, lattice_str: str):
    filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigen-dense-results.hdf5")

    h5file = h5py.File(filepath, "a")
    h5group = h5file["eigen-dense"]

    print("Reading result file")
    all_parameters = set()
    group_lookup = {}
    for child_name in h5group:
        g = h5group[child_name]
        idx = g.attrs["idx"]
        t = g.attrs["hopping"]
        U = g.attrs["interaction"]
        group_lookup[(idx, t, U)] = child_name
        all_parameters.add((t, U))

    #for irow, (idx, nup, ndn, tii, pii, pic, dim, root_idx, matrixtype) in sectors_df.iterrows():
    for irow, row in df.iterrows():
        if row.idx == row.root_idx: continue
        if row.dim == 0: continue

        assert(pd.isna(row.type))
        root_row = df.loc[row.root_idx-1]
        assert(root_row.dim != 0)

        if root_row.type not in ["small", "dense"]: continue

        for (t, U) in all_parameters:
            src_name = f"hopping={t}_idx={row.root_idx}_interaction={U}"
            dst_name = f"hopping={t}_idx={row.idx}_interaction={U}"
            if src_name not in h5group:
                print(f"missing: {src_name}")
                continue
            if dst_name in h5group:
            #if (idx, t, U) in group_lookup:
                print(f"exists: idx={row.idx}, t={t}, U={U}, root_idx={row.root_idx}")
                #continue
                del h5group[dst_name]
            assert(dst_name not in h5group)
            print(f"linking: idx={row.idx}, t={t}, U={U}, root_idx={row.root_idx}")
            dst_group = h5group.create_group(dst_name)
            dst_group.attrs["idx"] = idx
            dst_group.attrs["hopping"] = t
            dst_group.attrs["interaction"] = U
            dst_group.attrs["type"] = "dense-link"
            dst_group["eigenvalue"] = h5py.SoftLink(f"/eigen-dense/{src_name}/eigenvalue")
        
    h5file.flush()
    h5file.close()

def complete_sparse(df: pd.DataFrame, lattice_str: str):
    filepath = os.path.join(project_dir, "data", lattice_str, "eigen", "eigen-sparse-results.hdf5")

    h5file = h5py.File(filepath, "a")
    h5group = h5file["eigen-sparse"]

    print("Reading result file")
    all_parameters = set()
    group_lookup = {}
    for child_name in h5group:
        g = h5group[child_name]
        idx = g.attrs["idx"]
        t = g.attrs["hopping"]
        U = g.attrs["interaction"]
        group_lookup[(idx, t, U)] = child_name
        all_parameters.add((t, U))

    #for irow, (idx, nup, ndn, tii, pii, pic, dim, root_idx, matrixtype) in sectors_df.iterrows():
    for irow, row in df.iterrows():
        if row.idx == row.root_idx: continue
        if row.dim == 0: continue

        assert(pd.isna(row.type))
        root_row = df.loc[row.root_idx-1]
        assert(root_row.dim != 0)

        if root_row.type not in ["sparse"]: continue

        for (t, U) in all_parameters:
            src_name = f"hopping={t}_idx={row.root_idx}_interaction={U}"
            dst_name = f"hopping={t}_idx={row.idx}_interaction={U}"
            if src_name not in h5group:
                print(f"missing: {src_name}")
                continue
            if dst_name in h5group:
            #if (idx, t, U) in group_lookup:
                print(f"exists: idx={row.idx}, t={t}, U={U}, root_idx={row.root_idx}")
                #continue
                del h5group[dst_name]
            assert(dst_name not in h5group)
            print(f"linking: idx={row.idx}, t={t}, U={U}, root_idx={row.root_idx}")
            dst_group = h5group.create_group(dst_name)
            dst_group.attrs["idx"] = idx
            dst_group.attrs["hopping"] = t
            dst_group.attrs["interaction"] = U
            dst_group.attrs["type"] = "sparse-link"
            dst_group["eigenvalue"] = h5py.SoftLink(f"/eigen-sparse/{src_name}/eigenvalue")
        
    h5file.flush()
    h5file.close()

if __name__=='__main__':
    main()
