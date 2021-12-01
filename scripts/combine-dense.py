import os, sys
import json
import numpy as np
import scipy as sp
import fastavro
import pandas as pd
import argparse
from pathlib import Path

project_dir = (Path(os.path.dirname(__file__)) / "..").resolve()

with (project_dir / "src" / "eigenvalues.avsc").open("r") as f:
    schema = json.load(f)
parsed_schema = fastavro.parse_schema(schema)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("lattice_str", type=str)
    args = parser.parse_args()
    filenames = (project_dir / "data" / args.lattice_str / "eigen").glob("*.jsonl")

    existing_parameters = set()
    result_filepath = (project_dir / "data" / args.lattice_str / "eigen" / "dense-results.avro")
    readflag = result_filepath.is_file()
    with result_filepath.open("ab+") as fi:
        fi.seek(0)
        if readflag:
            avro_reader = fastavro.reader(fi)
            for record in avro_reader:
                t = record["hopping"]
                U = record["interaction"]
                idx = record["idx"]
                existing_parameters.add((t, U, idx))
        else:
            fastavro.writer(fi, parsed_schema, [])
        
        for filename in filenames:
            df = pd.read_json(filename, lines=True)
            nr = len(df)
            df = df[df.apply(lambda row: (row.hopping, row.interaction, row.idx) not in existing_parameters, axis=1)]
            nr2 = len(df)
            print(f"{nr} -> keep {nr2} / throw {nr-nr2}")
            df["timestamp"] = (df.timestamp - pd.Timestamp(1970, 1,1)) // pd.Timedelta("1ms")
            df.apply(lambda row: existing_parameters.add((row.hopping, row.interaction, row.idx)), axis=1)
            fastavro.writer(fi, parsed_schema, df.to_dict('records'))
            os.rename(filename, f"{filename}.processed")
            

if __name__=='__main__':
    main()

