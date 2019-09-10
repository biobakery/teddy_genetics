#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import argparse
from math import log10

from zopy.utils import iter_rows
from zopy.table2 import nesteddict2table

parser = argparse.ArgumentParser( )
parser.add_argument( "snpmap" )
parser.add_argument( "results" )
parser.add_argument( "mode", choices=["lenient", "strict"] )
args = parser.parse_args( )

#-------------------------------------------------------------------------------
# data munging
#-------------------------------------------------------------------------------

"""
1
rs61733845
0
1108138
"""

chroms = {}
snpmap = {}

for row in iter_rows( args.snpmap ):
    chrom, snp, other, coord = row
    chrom = int( chrom )
    if chrom > 22:
        continue
    coord = int( coord )
    snpmap[snp] = [chrom, coord]
    chroms.setdefault( chrom, [] ).append( coord )

shifts = {}
for i in sorted( chroms ):
    shifts[i] = -min( chroms[i] )
for i in sorted( chroms ):
    if i > 1:
        shifts[i] = max( chroms[i-1] ) + shifts[i-1]

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

pvals = {}
for row in iter_rows( args.results ):
    bug = row[0].replace( "s__", "" )
    snp = row[1]
    add_pval = float( row[7].split( "|" )[2] )
    int_pval = float( row[8].split( "|" )[2] )
    choices = [add_pval, int_pval]
    pval = min( choices ) if args.mode == "lenient" else max( choices )
    pvals.setdefault( bug, {} )[snp] = pval

#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------
    
import matplotlib as mpl
mpl.use( "Agg" )
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import matplotlib.pyplot as plt

import zopy.mplutils2 as mu

def plot( bug ):

    fig = plt.figure( )
    fig.set_size_inches( 12, 5 )

    ax = plt.subplot( 111 )

    xticks = []
    for i in sorted( chroms ):
        pos = 0
        pos += min( chroms[i] )
        pos += max( chroms[i] )
        pos /= 2.0
        pos += shifts[i]
        xticks.append( pos )
        
    ax.set_xticks( xticks )
    ax.set_xticklabels( sorted( chroms ) )

    to_label = []
    xmax = 0
    ymax = 0
    for m in [0, 1]:
        xx = []
        yy = []
        for snp, pval in pvals[bug].items( ):
            if snp not in snpmap:
                continue
            chrom, coord = snpmap[snp]
            if chrom % 2 != m:
                continue            
            pval = -log10( pval )
            xx.append( shifts[chrom] + coord )
            yy.append( pval )
            to_label.append( [yy[-1], xx[-1], snp] )
        ax.scatter( xx, yy, color="blue" if m == 1 else "cornflowerblue", edgecolor="none", alpha=0.5 )
        xmax = max( xmax, max( xx ) )
        ymax = max( ymax, max( yy ) )

    to_label.sort( )
    for y, x, snp in to_label[-10:]:
        ax.text( x, y, "  " + snp, ha="left", va="center", size=8 )
        
    ax.set_xlim( 0, xmax )
    ax.set_ylim( 0, ymax )
    ax.set_xlabel( "chromosome / coordinate" )
    ax.set_ylabel( "-log10( p )" )
    
    mu.tick_params( ax )
    mu.margin( ax )
    mu.hline( ax, -log10( 5e-8 ), color="red", linestyle="--" )

    title_bug = "log10( {} )".format( bug ) if "_log10" in args.results else bug
    title = """{} ~ time + snp + time * snp + (1|subject)
    {} of additive and time-linked SNP effect p-values (i.e. {} mode)"""
    ax.set_title( title.format( title_bug, "MINIMUM" if args.mode == "lenient" else "MAXIMUM", args.mode ) )
    
    plt.tight_layout( )

    path_bug = "{}_log10".format( bug ) if "_log10" in args.results else bug                                    
    plt.savefig( "manhattan.{}.{}.png".format( path_bug, args.mode ) )
    plt.close( )

for bug in pvals:
    plot( bug )
