#!/usr/bin/env python

import os
import sys
import argparse

from zopy.utils import tprint, qw, say
from zopy.table2 import table, nesteddict2table

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "species" )
parser.add_argument( "mode", choices=["init", "prev"] )
args = parser.parse_args( )

T = table( args.species )
M, F = T.metasplit( "Day" )
F.float( )

def bc( da, db ):
    ret = 0
    bugs = set( da )
    bugs.update( set( db ) )
    for bug in bugs:
        aval = da.get( bug, 0 )
        bval = db.get( bug, 0 )
        ret += min( aval, bval )
    return 1 - ret

kids = {}
for sample in M.colheads:
    kid = M.entry( "Subject", sample )
    day = int( M.entry( "Day", sample ) )
    kids.setdefault( kid, [] ).append( [day, sample] )

fname = "bc_dists_{}.details.tsv".format( args.mode )
fh = open( fname, "w" )

headers = qw( """
Subject
IndexA
DayA
SampleA
IndexB
DayB
SampleB
Delta
BCDist
""" )
tprint( *headers, file=fh )

T = {}
inner = T.setdefault( fname.split( "." )[0], {} )

kid_count = 0
for kid in sorted( kids ):
    kid_count += 1
    say( "Working on kid", kid_count, "of", len( kids ) )
    pairs = kids[kid]
    pairs.sort( )
    if len( pairs ) < 2:
        continue
    i = 1
    while i < len( pairs ):
        a = {"init":0, "prev":i-1}[args.mode]
        b = i
        da = F.coldict( pairs[a][1] )
        db = F.coldict( pairs[b][1] )
        dist = bc( da, db )       
        row = [kid, a] + pairs[a] + [b] + pairs[b] + [pairs[b][0] - pairs[a][0]] + [dist]
        tprint( *row, file=fh )
        inner[pairs[b][1]] = dist
        # critical
        i += 1

fh.close( )

T = nesteddict2table( T )
T.metamerge( M )
T.data[0][0] = "#"
T.dump( fname.replace( ".details", "" ) )
