import sys, os
import h5py
import pyarrow.feather

project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

lattice_str = "square-(4,0)x(0,4)"
filepath = "/gpfs/home/kyungminlee_42/Projects/HubbardFTLM/data/square-(4,0)x(0,4)/eigen/eigen-dense-results.hdf5"

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

print("Reading sectors")
sectors_df = pyarrow.feather.read_table(os.path.join(project_dir, "data", f"sectors-{lattice_str}.arrow")).to_pandas()

for irow, (idx, nup, ndn, tii, pii, pic, dim, root_idx) in sectors_df.iterrows():
    if idx != root_idx:
        for (t, U) in all_parameters:
            src_name = f"hopping={t}_idx={root_idx}_interaction={U}"
            dst_name = f"hopping={t}_idx={idx}_interaction={U}"
            if (idx, t, U) in group_lookup:
                print(f"exists: idx={idx}, t={t}, U={U}, root_idx={root_idx}")
                #continue
                del h5group[dst_name]
            #if (root_idx, t, U) in group_lookup:
            #assert(src_name == group_lookup[root_idx, t, U])
            assert(dst_name not in h5group)
            #src_group = h5group[src_name]
            dst_group = h5group.create_group(dst_name)
            dst_group.attrs["idx"] = idx
            dst_group.attrs["hopping"] = t
            dst_group.attrs["interaction"] = U
            #dst_group.attrs["description"] = f"copied from {root_idx}"
            dst_group["eigenvalue"] = h5py.SoftLink(f"/eigen-dense/{src_name}/eigenvalue")
            #else:
            #    print(f"missing: idx={idx}, t={t}, U={U}, root_idx={root_idx}")
        
h5file.flush()
h5file.close()
#print(group_lookup)
