#!/usr/bin/env python

import os
import sys
import re
import argparse

from zopy.utils import iter_rowdicts, tprint
from zopy.table2 import table

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "species", help="" )
args = parser.parse_args( )

target = 12 * (365/12.0)
window = 60

T = table( args.species )
M, T = T.metasplit( "Day" )
T.float( )
T.groupby( lambda x: x.replace( "s__", "" ).split( "_" )[0], sum )

ss = M.row( "Subject" )
dd = M.row( "Day" )
vv = T.row( "Bacteroides" )

diffs = {}
values = {}
for s, d, v in zip( ss, dd, vv ):
    d = int( d )
    diff = abs( d - target )
    if diff < window:
        if s not in diffs or diffs[s] > diff:
            diffs[s] = diff
            values[s] = v
        
for k, v in values.items( ):
    tprint( k, v )
