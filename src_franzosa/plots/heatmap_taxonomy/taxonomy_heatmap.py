#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import argparse
import random
from math import log10

from numpy import mean

from zopy.utils import iter_rows, qw
from zopy.table2 import table, nesteddict2table

parser = argparse.ArgumentParser( )
parser.add_argument( "species", help="" )
args = parser.parse_args( )

T = table( args.species )

# simplfy for testing
"""
ok = {c for c in T.colheads if random.random( ) < 0.1}
T.select( "headers", ok, transposed=True )
"""

# adjust Day to ordinal
T.apply_rowheads( lambda x: x if x != "Day" else "Time" )
for c in T.colheads:
    day = T.entry( "Time", c )
    day = int( day )
    mos = 6
    ival = 6 * (365/12.0)
    while day > ival and mos <= 36:
        mos += 6
        day -= ival
    label = "{:02d} mos".format( mos ) if mos <= 36 else "99 mos"
    T.set( "Time", c, label )

# remove rare bugs
M, T = T.metasplit( "Time" )
T.apply_rowheads( lambda x: x.replace( "s__", "" ).replace( "_", " " ) )
T.grep( "headers", "unclassified", invert=True )
T.float( )
means = {}
for r, row in T.iter_rows( ):
    means[r] = mean( row )
ok = {r for r in sorted( means, key=lambda x: means[x] )[-30:]}
T.select( "headers", ok )
T.rowsort( order=sorted( T.rowheads, key=lambda x: means[x] ) )
T.metamerge( M )

T.dump( "heatmap.tsv" )

command = qw( """
python /home/efranzosa/hg/zopy/scripts/heatmap.py heatmap.tsv
--overrides rownames_c:10
--title "TEDDY Taxonomy"
--lastrow Time
--colmeta Time
--colmeta-colors Blues
--cmap bbcry
--units "Relative abundance"
--transform log10
--rowsort mean
--colsort metadata
--limits -6 0
--vscale
--output "taxonomy-overview.pdf"
--rowlabel "top 30 named species, mean abund."
--collabel "metagenomes"
""" )
command = " ".join( command )
os.system( command )
