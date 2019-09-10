#!/usr/bin/env python

import os
import sys
import argparse

from zopy.utils import iter_rows, tprint

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "data" )
parser.add_argument( "h", type=int )
parser.add_argument( "r", type=int )
args = parser.parse_args( )

name = os.path.split( args.data )[1]
h = args.h
r = args.r
f = 0

headers = []
rows = None

for row in iter_rows( args.data ):
    if len( headers ) < h:
        headers.append( row )
        continue
    elif len( headers ) == h and rows is None:
        rows = [row[:] for row in headers]
    if len( rows ) == r + h:
        f += 1
        fh = open( "{}.chunk-{:05d}".format( name, f ), "w" )
        for row2 in rows:
            tprint( *row2, file=fh )
        fh.close( )
        rows = [row[:] for row in headers]
    rows.append( row )

f += 1
fh = open( "{}.chunk-{:05d}".format( name, f ), "w" )
for row2 in rows:
    tprint( *row2, file=fh )
fh.close( )
