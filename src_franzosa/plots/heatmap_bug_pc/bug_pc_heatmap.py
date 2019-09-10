#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import argparse
from math import log10

from zopy.utils import iter_rows
from zopy.table2 import nesteddict2table

parser = argparse.ArgumentParser()
parser.add_argument( "results", help="" )
parser.add_argument( "mode", choices=["add", "int"], help="" )
parser.add_argument( "--simplify", action="store_true", help="" )
args = parser.parse_args()

"""
0 s__Akkermansia_muciniphila
1 PC1
2 Day|1.642854810976314e-06
3 SNP|0.3726319584622303
4 Day:SNP|0.31764727039982055
5 (Intercept)|0.01735867420137683|2.831566616544884e-27
6 Day|1.1021459007523974e-05|1.642854810976314e-06
7 SNP|-0.0652547950687024|0.3726319584622303
8 Day:SNP|0.00011036889339222005|0.31764727039982055
"""

coefs = {}
pvals = {}
for row in iter_rows( args.results ):
    bug = row[0].replace( "s__", "" ).replace( "_", " " )
    pc = row[1]
    pc = "PC{:02d}".format( int( pc.replace( "PC", "" ) ) )
    index = 7 if args.mode == "add" else 8
    coef = float( row[index].split( "|" )[1] )
    pval = float( row[index].split( "|" )[2] )
    coefs.setdefault( bug, {} )[pc] = (1 if coef > 0 else -1) * abs( log10( pval ) ) 
    pvals.setdefault( bug, {} )[pc] = pval

coefs = nesteddict2table( coefs )
coefs.data[0][0] = "#"

ok_bugs = set( )
ok_pcs = set( )

tests = len( coefs.rowheads ) * len( coefs.colheads )
with open( "BonferroniSig.txt", "w" ) as fh:
    print( "#kwargs\tcolor:white\tedgecolor:none", file=fh )
    for bug in pvals:
        for pc, p in pvals[bug].items( ):
            if p <= 0.05 / tests:
                print( "\t".join( [bug, pc] ), file=fh )
                ok_bugs.add( bug )
                ok_pcs.add( pc )

if args.simplify:
    coefs.select( "headers", ok_bugs )
    coefs.select( "headers", ok_pcs, transposed=True )
coefs.dump( "coefs.tsv" )
                
command = """
python /home/efranzosa/hg/zopy/scripts/heatmap.py coefs.tsv
--overrides colnames_r:3 rownames_c:10
--title "Major species vs. SNP PCs: {0}"
--cmap RdBu_r
--units "sign( coef ) * -log10( pval )"
--limits -5 5
--dots BonferroniSig.txt
--vscale
--rowsort euclidean
--colsort euclidean
--output "heatmap_bugs-vs-pcs_{1}{2}.pdf"
"""
command = command.replace( "\n", " " )
command = command.format(
    "Additive Effect" if args.mode == "add" else "Time Interaction",
    "additive" if args.mode == "add" else "time-interaction",
    "_simplfied" if args.simplify else "",
)
os.system( command )
