#!/usr/bin/env python

import argparse
import glob
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-ext')
parser.add_argument('-out')
args = parser.parse_args()

files = glob.glob("*{}".format(args.ext))

df_from_each_file = (pd.read_csv(f) for f in files)

concatenated_df   = pd.concat(df_from_each_file, ignore_index=True)

concatenated_df.to_csv(args.out, sep="\t", index=False)
