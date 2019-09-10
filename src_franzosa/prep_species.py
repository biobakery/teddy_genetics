#!/usr/bin/env python

import os
import sys
import argparse

from math import log

from zopy.table2 import table

# constants
c_min_peak = 0.01
c_min_frac = 0.10

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "phlan" )
parser.add_argument( "lastmeta" )
args = parser.parse_args( )

# combine
P = table( args.phlan )
M, P = P.metasplit( args.lastmeta )

# filter
P.apply_rowheads( lambda x: x.split( "|" )[-1] )
P.grep( "headers", "s__" )
P.float( )
P.apply_entries( lambda x: x / 100.0 )

# compute smoothing factors / remove empty samples
bad = set( )
eps = {}
for sample, values in P.iter_cols( ):
    if sum( values ) == 0:
        bad.add( sample )
    else:
        eps[sample] = 0.5 * min( [k for k in values if k > 1e-20] )
P.delete( "headers", bad, transposed=True )

# save all species
P.unfloat( )
P.rowsort( )
P.metamerge( M )
P.dump( "species.tsv" )

# isolate major species
M, P = P.metasplit( args.lastmeta )
P.float( )
kids = M.row( "Subject" )
bad = set( )
for s in P.rowheads:
    peaked = set( )
    for k, v in zip( kids, P.row( s ) ):
        if v >= c_min_peak:
            peaked.add( k )
    if float( len( peaked ) ) / len( set( kids ) ) < c_min_frac:
        bad.add( s )
P.unfloat( )
P.delete( "headers", bad )
P.metamerge( M )
P.dump( "major_species.tsv" )

# transform major species
M, P = P.metasplit( args.lastmeta )
P.float( )
for r, c, v in P.iter_entries( ):
    P.set( r, c, log( v + eps[c] ) / log( 10.0 ) )
P.unfloat( )
P.metamerge( M )
P.dump( "major_species_log10.tsv" )
