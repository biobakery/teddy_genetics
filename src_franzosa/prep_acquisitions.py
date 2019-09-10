#!/usr/bin/env python

import os
import sys
import argparse

from zopy.utils import tprint, qw, say
from zopy.table2 import table, nesteddict2table

# constants
c_detect = 0.01
c_early = 365
c_later = 1000000

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "species" )
args = parser.parse_args( )

T = table( args.species )
M, F = T.metasplit( "Day" )
F.float( )

# reorganize
pairs = {}
for sample in M.colheads:
    subject, day = M.col( sample )
    day = int( day )
    if day > c_later:
        continue
    pairs.setdefault( subject, [] ).append( [day, sample] )

# subset subjects
ok = set( )
for subject in pairs:
    test1 = any( [d[0] <= c_early for d in pairs[subject]] )
    test2 = any( [d[0] <= c_later for d in pairs[subject]] )
    if test1 and test2:
        ok.add( subject )
pairs = {k:v for k, v in pairs.items( ) if k in ok}

# filter bugs
data = {}
for bug in F.rowheads:
    inner = data.setdefault( bug, {} )
    for subject in pairs:
        phenos = ["2_never"]
        for [day, sample] in sorted( pairs[subject] ):
            if F.entry( bug, sample ) >= c_detect:
                phenos.append( "0_early" if day < c_early else "1_later" )
        phenos.sort( )
        inner[subject] = phenos[0]

# output
T = nesteddict2table( data )
T.data[0][0] = "#"
T.dump( "acquisitions.tsv" )
