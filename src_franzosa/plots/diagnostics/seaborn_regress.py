#!/usr/bin/env python

import os
import sys
import re
import argparse
from collections import Counter

import matplotlib as mpl
mpl.use( "Agg" )
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns

from zopy.table2 import table

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

def clean( df ):
    gg = df["genotype"]
    ok = [i for i, g in enumerate( gg ) if g != "00"]
    for k, v in df.items( ):
        df[k] = [v[i] for i in ok]                                   

def plot( mbf, snp ):
    snpmap = S.rowdict( snp )
    kidmap = M.rowdict( "Subject" )
    snpmap = {c:snpmap.get( v, "??" ) for c, v in kidmap.items( )}
    df = {}
    df["divergence" ] = [F.entry( mbf, c ) for c in F.colheads]
    df["day"] = [int( M.entry( "Day", c ) ) for c in F.colheads]
    df["genotype"] = [snpmap[c] for c in F.colheads]
    clean( df )
    order = Counter( df["genotype"] )
    order = sorted( order, key=lambda x: -order[x] )
    df = pd.DataFrame.from_dict( df )
    fig = plt.figure()
    fig.set_size_inches( 6, 4 )
    sns.lmplot( x="day", y="divergence", hue="genotype", data=df,
                scatter=False,
                palette="inferno",
                hue_order=order,
                ci=95,
    )
    plt.title( snp )
    plt.ylim( 0, 1 )
    plt.xlim( 0, 1000 )
    plt.savefig( "regress.{}.{}.pdf".format( mbf, snp ) )
    plt.close( )
    
c = 0
for snp in S.rowheads:
    c += 1
    print c, "of", len( S.rowheads )
    plot( "bc_dists_init", snp )
