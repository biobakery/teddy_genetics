#!/usr/bin/env python

import os
import sys
import re
import argparse

from zopy.utils import tprint
from zopy.table2 import table
from zopy.dictation import polymap
from zopy.enrichments import rank_enrich

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "loadings", help="" )
parser.add_argument( "snpmap", help="" )
args = parser.parse_args( )

T = table( args.loadings )
annotations = polymap( args.snpmap, key=1 )

headers = None
for pc in T.colheads:
    d = T.coldict( pc )
    d = {k:float( v ) for k, v in d.items( ) if v != "nan"}
    results = rank_enrich( d, annotations,
                           min_overlap=50,
                           fdr=0.1,
                           intersect_annotated=True,
                           verbose=True,
    )
    if headers is None:
        headers = ["#PC"] + results[0].fields
        tprint( *headers )
    for R in results:
        row = [pc] + R.row( )
        tprint( *row )
