#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse

from rpy2 import robjects
from rpy2.robjects import FloatVector, StrVector
from rpy2.robjects.packages import importr

from zopy.table2 import table

parser = argparse.ArgumentParser( )
parser.add_argument( "outcomes" )
parser.add_argument( "genetics" )
parser.add_argument( "results" )
args = parser.parse_args( )

fh = open( args.results, "w" )

O = table( args.outcomes )
G = table( args.genetics )
O.select( "Subject", set( G.colheads ), transposed=True )

M, F = O.metasplit( "Day" )
F.float( )

base = importr( "base" )
stats = importr( "stats" )
lmer = importr( "lmerTest" )

robjects.globalenv["Subject"] = StrVector( M.row( "Subject" ) )
robjects.globalenv["Day"] = FloatVector( [float( k ) for k in M.row( "Day" )] )

for bug in F.rowheads:
    robjects.globalenv["Bug"] = FloatVector( F.row( bug ) )
    for snp in G.rowheads:
        # general
        outline = [bug, snp]
        try:
            # per-snp prep
            d = G.rowdict( snp )
            snp_vals = [d.get( k, "??" ) for k in M.row( "Subject" )]
            robjects.globalenv["SNP"] = StrVector( snp_vals )
            # run the model
            model = lmer.lmer( "Bug ~ Day + SNP + Day * SNP + (1|Subject)" )
            # get overall effect from anova
            anova = stats.anova( model )
            coefs = anova.rownames
            pvals = anova.rx2( "Pr(>F)" )
            details = ["|".join( [str( k ) for k in items] ) for items in zip( coefs, pvals )]
            outline += details
            # extract levelwise details from summary
            summary = base.summary( model )
            a = summary.rx2( "coefficients" )
            levels = a.rownames
            coefs = a.rx( True, 1 )
            pvals = a.rx( True, 5 )
            details = ["|".join( [str( k ) for k in items] ) for items in zip( levels, coefs, pvals )]
            outline += details
        except:
            outline += ["# MODEL FAILED"]
        print( "\t".join( outline ), file=fh )

fh.close( )

