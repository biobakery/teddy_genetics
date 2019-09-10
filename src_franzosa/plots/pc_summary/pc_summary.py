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
parser.add_argument( "eigen" )
parser.add_argument( "results" )
args = parser.parse_args( )

#-------------------------------------------------------------------------------
# data munging
#-------------------------------------------------------------------------------

eigen = {}
for row in iter_rows( args.eigen ):
    pc = "PC{:02d}".format( len( eigen ) + 1 )
    eigen[pc] = float( row[0] )
eigen = {k:v/max( eigen.values( ) ) for k, v in eigen.items( )}   

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

adds = {}
ints = {}
for row in iter_rows( args.results ):
    pc = row[1]
    pc = "PC{:02d}".format( int( pc.replace( "PC", "" ) ) )   
    adds[pc] = float( row[7].split( "|" )[2] )
    ints[pc] = float( row[8].split( "|" )[2] )

#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------
    
import matplotlib as mpl
mpl.use( "Agg" )
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import matplotlib.pyplot as plt

import zopy.mplutils2 as mu

fig = plt.figure( )
fig.set_size_inches( 8, 5 )

axes = []
ax = plt.subplot( 311 )
axes.append( ax )

kwargs_ylab = {"rotation":0, "rotation_mode":"anchor", "ha":"right", "va":"center"}

colors = ["blue" if i % 2 == 0 else "cornflowerblue" for i in range( len( eigen ) )]

mu.barplot( ax, [eigen[pc] for pc in sorted( eigen )], labels=None, colors=colors )
ax.set_ylabel( "Relative\nimportance\n(scaled eigenvalue)", **kwargs_ylab )

ax = plt.subplot( 312 )
axes.append( ax )

mu.barplot( ax, [-log10( adds[pc] ) for pc in sorted( eigen )], labels=None, colors=colors )
ax.set_ylabel( "Additive PC effect\non divergence\n-log10( P-value )", **kwargs_ylab )

ax = plt.subplot( 313 )
axes.append( ax )

mu.barplot( ax, [-log10( ints[pc] ) for pc in sorted( eigen )], labels=None, colors=colors )
ax.set_ylabel( "Time-linked PC effect\non divergence\n-log10( P-value )", **kwargs_ylab )

for i, ax in enumerate( axes ):
    mu.tick_params( ax )
    if i < 2:
        mu.hide_xaxis( ax )
    else:
        ax.set_xticks( range( len( eigen ) ) )
        ax.set_xticklabels( sorted( eigen ), rotation=35, rotation_mode="anchor", ha="right", va="center" )
    
plt.tight_layout( )
plt.savefig( "pc_summary.pdf" )
plt.close( )
