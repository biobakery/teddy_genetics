#!/usr/bin/env python

import os
import sys
import re
import argparse
from collections import Counter

from numpy import polyfit, array

import matplotlib as mpl
mpl.use( "Agg" )
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import matplotlib.pyplot as plt

from zopy.table2 import table
import zopy.mplutils2 as mu

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "mbfs", help="" )
parser.add_argument( "snps", help="" )
parser.add_argument( "limit", help="" )
args = parser.parse_args( )

L = table( args.limit )
B = table( args.mbfs )
S = table( args.snps )
S.select( "headers", L.rowheads )

subs = set( S.colheads ).__and__( set( B.row( "Subject" ) ) )
S.select( "headers", subs, transposed=True )
B.select( "Subject", subs, transposed=True )

M, F = B.metasplit( "Day" )
F.float( )

def plot( mbf, snp ):
    snpmap = S.rowdict( snp )
    kidmap = M.rowdict( "Subject" )
    snpmap = {c:snpmap.get( v, "??" ) for c, v in kidmap.items( )}
    vv = [F.entry( mbf, c ) for c in F.colheads]
    dd = [int( M.entry( "Day", c ) ) for c in F.colheads]
    ll = [snpmap[c] for c in F.colheads]
    # reshape
    dv = {}
    dd = {}
    for v, d, l in zip( vv, dd, ll ):
        dv.setdefault( l, [] ).append( v )
        dd.setdefault( l, [] ).append( d )
    fig = plt.figure()
    fig.set_size_inches( 3, 4 )
    ax = plt.subplot( 111 )
    gens = Counter( ll )
    order = sorted( gens, key=lambda x: -gens[x] )
    colors = mu.ncolors( len( order ) + 1, "Blues_r" )[:-1]
    colors = {o:c for o, c in zip( order, colors)}

    # plotting
    for l in dv:
        m, b = polyfit( dd[l], dv[l], 1 )
        xx = array( dd[l] )
        yy = m * xx + b
        ax.scatter( xx, yy, edgecolor="none", color=colors[l], label=l )

    ax.set_ylim( 0, 1 )
    ax.set_title( snp )
    ax.set_ylabel( "Divergence from first visit\n(species-level Bray-Curtis distance)" )
    ax.set_xlabel( "Age (days)" )
    ax.legend( )

    mu.tick_params( ax )
    mu.margin( ax )
    plt.tight_layout( )
    plt.savefig( "boxplot.{}.{}.pdf".format( mbf, snp ) )
    plt.close( )

c = 0
for snp in S.rowheads[0:1]:
    c += 1
    print c, "of", len( S.rowheads )
    plot( "bc_dists_init", snp )
