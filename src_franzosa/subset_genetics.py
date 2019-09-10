#!/usr/bin/env python

import os
import sys
import argparse

from zopy.utils import iter_rows
from zopy.table2 import table

"""
bc_dists_init
rs1320571
Day|1.8673952473457896e-10
SNP|0.0031287399459056643
Day:SNP|8.822016330162701e-05
(Intercept)|0.0503814762174525|0.76022577679799
Day|0.0010832352925434629|3.7988004757133905e-05
SNPAG|0.5398957290054202|0.00117960317375227
SNPGG|0.5138193831272909|0.001904932200939625
Day:SNPAG|-0.000814692846950005|0.001989626397885419
Day:SNPGG|-0.000754194045396658|0.004133830028756279
"""

# constant
crit_p = 5e-8

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "genetics" )
parser.add_argument( "--results", nargs="+", required=True )
parser.add_argument( "--mode", choices=["strict", "relaxed"], required=True )
args = parser.parse_args( )

# parse results
ok = set( )
for p in args.results:
    for row in iter_rows( p ):
        snp = row[1]
        snp_p = float( row[3].split( "|" )[1] )
        int_p = float( row[4].split( "|" )[1] )
        if int_p <= crit_p:
            if snp_p <= crit_p or args.mode == "relaxed":
                ok.add( snp )

# combine
T = table( args.genetics )
T.select( "headers", ok )
T.dump( "genetics_subset-{}.tsv".format( args.mode ) )
