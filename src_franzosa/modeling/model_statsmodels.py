#!/usr/bin/env python

import os
import sys
import argparse

import numpy as np
import statsmodels.api as sm

from zopy.table2 import table

parser = argparse.ArgumentParser( )
parser.add_argument( "outcomes" )
parser.add_argument( "genetics" )
args = parser.parse_args( )

O = table( args.outcomes )
G = table( args.genetics )
O.select( "Subject", set( G.colheads ), transposed=True )

M, F = O.metasplit( "Day" )
F.float( )

df = {}
df["Subject"] = M.row( "Subject" )
df["Day"] = [float( k ) for k in M.row( "Day" )]

for bug in F.rowheads:
    df["Bug"] = F.row( bug )
    for snp in G.rowheads:
        d = G.rowdict( snp )
        df["SNP"] = [d[k] for k in M.row( "Subject" )]
        model = sm.MixedLM.from_formula( "Bug ~ Day + SNP + Day * SNP", df, groups=df["Subject"] )
        result = model.fit( )
        print( result.summary( ) )
