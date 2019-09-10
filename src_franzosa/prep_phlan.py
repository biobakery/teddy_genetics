#!/usr/bin/env python

import os
import sys
import argparse

from zopy.table2 import table

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "old" )
parser.add_argument( "new" )
parser.add_argument( "meta" )
args = parser.parse_args( )

# combine
P = table( args.old )
P.merge( table( args.new ) )
P.na2zero( )
P.rowsort( )
P.apply_colheads( lambda x: x.split( "_" )[0] )

# attach metadata / dump samples without metadata
M = table( args.meta )
M.transpose( )
P.select( "headers", set( M.colheads ), transposed=True )
P.metamerge( M )
P.dump( "all_phlan.tsv" )
