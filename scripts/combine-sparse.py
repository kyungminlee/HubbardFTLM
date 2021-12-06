import os, sys
import json
import numpy as np
import scipy as sp
import fastavro
import pandas as pd
import argparse
from pathlib import Path
import progressbar

project_dir = (Path(os.path.dirname(__file__)) / "..").resolve()

with (project_dir / "src" / "eigenvalues-sparse.avsc").open("r") as f:
    schema = json.load(f)
parsed_schema = fastavro.parse_schema(schema)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    args = parser.parse_args()
    input_paths = (project_dir / "data" / args.lattice_str / "eigen" / "sparsedata").glob("*.qed")

    result_filepath = (project_dir / "data" / args.lattice_str / "eigen" / "sparse-results.avro")

    input_paths = list(input_paths)
    
    data = {}
    for input_path in progressbar.progressbar(input_paths):
        df = pd.read_json(str(input_path / "eigenvalues.jsonl"), lines=True)
        for _, row in df.iterrows():
            iden = (row.hopping, row.interaction, row.idx, row.seed, row.krylovdim)
            if iden not in data:
                data[iden] = {
                    'hopping': row.hopping,
                    'interaction': row.interaction,
                    'idx': row.idx,
                    'type': 'sparse',
                    'seed': row.seed,
                    'krylovdim': row.krylovdim,
                    'coefficients': [],
                    'eigenvalues': [],
                    'timestamp': 0,
                    'description': '',
                    # 'sampledim': 0,
                }
            data[iden]['coefficients'].extend(row.coefficients)
            data[iden]['eigenvalues'].extend(row.eigenvalues)
            t = (row.timestamp - pd.Timestamp(1970, 1,1)) // pd.Timedelta("1ms")
            data[iden]['timestamp'] = max(data[iden]['timestamp'], t)
            data[iden]['description'] += row.description
            # data[iden]['sampledim'] += row.krylovdim

    with result_filepath.open("wb+") as fi:
        fastavro.writer(fi, parsed_schema, data.values())  

if __name__=='__main__':
    main()

