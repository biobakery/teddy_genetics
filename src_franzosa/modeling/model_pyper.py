#!/usr/bin/env python

import os
import sys
import argparse

import pandas as pd
import pyper

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

r = pyper.R( use_pandas=True )
r( "library(lmerTest)" )

for bug in F.rowheads[0:1]:
    df["Bug"] = F.row( bug )
    for snp in G.rowheads[0:1]:
        d = G.rowdict( snp )
        df["SNP"] = [d[k] for k in M.row( "Subject" )]
        r.assign( "df", pd.DataFrame( df ) )
        result = r( "anova( lmer( 'Bug ~ Day + SNP + Day * SNP + (1|Subject)', data=df ) )" )
        print( result )
