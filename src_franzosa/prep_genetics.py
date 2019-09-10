#!/usr/bin/env python

import os
import sys
import argparse
from collections import Counter

from zopy.utils import iter_rows, say
from zopy.table2 import table, nesteddict2table

# constants
c_minor_alleles = 0.05

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "ped" )
parser.add_argument( "map" )
parser.add_argument( "species" )
parser.add_argument( "--test", action="store_true" )
args = parser.parse_args( )

# restrict to samples with mgx data
phlan = set( table( args.species ).row( "Subject" ) )


"""
1  rs61733845  0   1108138
"""

snps = []
for R in iter_rows( args.map ):
    #snps.append( ":".join( [R[1], R[0], R[3]] ) )
    snps.append( R[1] )

"""
0       906136  0       0       0       -9      G G     G G     C A
"""

data = {}
counter = 0
for R in iter_rows( args.ped ):
    subject = R[1]
    if subject not in phlan:
        continue
    counter += 1
    say( "Processing subject #", counter )
    if args.test and counter > 100:
        break
    calls = R[6:]
    assert len( calls ) == len( snps ), exit( "SNP mismatch" )
    for s, c in zip( snps, calls ):
        inner = data.setdefault( s, {} )
        inner[subject] = c.replace( " ", "" )

# find "common" snps (minor alleles >=5%)
ok = set( )
for snp in data:
    C = Counter( "".join( data[snp].values( ) ) )
    if max( C.values( ) ) / float( sum( C.values( ) ) ) < 1 - c_minor_alleles:
        ok.add( snp )

data = {snp:d for snp, d in data.items( ) if snp in ok}
T = nesteddict2table( data, aRowheads=[snp for snp in snps if snp in data] )
T.data[0][0] = "#"
T.dump( "genetics.tsv" )
