#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import argparse
from math import log10

from scipy.stats import spearmanr
from scipy.stats.mstats import mquantiles

import matplotlib as mpl
mpl.use( "Agg" )
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.sans-serif"] = "Arial"
import matplotlib.pyplot as plt

from zopy.utils import iter_rows, path2name
from zopy.table2 import table, nesteddict2table
from zopy.dictation import col2dict
import zopy.mplutils2 as mu

#-------------------------------------------------------------------------------
# cli
#-------------------------------------------------------------------------------

parser = argparse.ArgumentParser( )
parser.add_argument( "pca" )
parser.add_argument( "colorby" )
args = parser.parse_args( )

#-------------------------------------------------------------------------------
# munging
#-------------------------------------------------------------------------------

T = table( args.pca )
T.float( )
T.apply_rowheads( lambda x: x if len( x ) == 4 else "PC0" + x[-1] )

colorby_raw = col2dict( args.colorby, value=1, func=float )
T.select( "headers", colorby_raw, transposed=True )

def qform( d ):
    ret = {}
    qq = list( mquantiles( list( d.values( ) ) ) )
    qq.append( max( d.values( ) ) )
    for k, v in d.items( ):
        for i, q in enumerate( qq ):
            if v <= q:
                v = "Q{}: <{:.2g}".format( i+1, q )
                ret[k] = v
                break
    return ret
colorby = qform( colorby_raw )

qq = sorted( set( colorby.values( ) ) )
colors = mu.ncolors( 1 + len( qq ), "Blues" )[1:]
cmap = {q: c for q, c in zip( qq, colors )}

for c in T.colheads:
    if c not in colorby:
        colorby[c] = "Missing"
                                      
#-------------------------------------------------------------------------------
# plotting
#-------------------------------------------------------------------------------
    
def plot( a, b, da, db ):

    fig = plt.figure( )
    fig.set_size_inches( 8, 4 )

    P = mu.Plotrix( [3], [1] )
    ax = P.axes[0][0]
    xx = [da[k] for k in sorted( da )]
    yy = [db[k] for k in sorted( da )]
    cc = [cmap[colorby[k]] for k in sorted( da )]
    for x, y, c in zip( xx, yy, cc ):
        ax.scatter( x, y, color=c, edgecolor="none" )

    mu.dummy( P.legend )
    L = mu.Legendizer( P.legend, vscale=1.0 )
    L.color_guide( cmap, title=path2name( args.colorby ) )
    L.draw( )
    
    ax.set_xlabel( a )
    ax.set_ylabel( b )
    mu.tick_params( ax )
    mu.tick_grid( ax )

    ax.set_title( "Spearman = {:.3f}".format( spearmanr( xx, yy )[0] ), loc="left" )
    
    plt.tight_layout( )
    plt.savefig( "pc_scatters.{}.{}_vs_{}.pdf".format( path2name( args.colorby ), a, b ) )
    plt.close( )

for i in range( 3 ):
    for j in range( 3 ):
        if i < j:
            a = T.rowheads[i]
            b = T.rowheads[j]
            plot( a, b, T.rowdict( a ), T.rowdict( b ) )

"""
for i in range( 5 ):
    a = T.rowheads[i]
    plot( a, path2name( args.colorby ), T.rowdict( a ), colorby_raw )
"""
