#!/usr/bin/env python

import os
import sys
import re

"""
major_species_tsv_chunk-00089.genetics_subset-strict_tsv.results.txt
"""

m = {}
for p in sys.argv[1:]:
    name = os.path.split( p )[1]
    name = name.replace( "_tsv", "" )
    name = re.sub( "_chunk-[0-9]+", "", name )
    m.setdefault( name, [] ).append( p )

for n in sorted( m ):
    fh = open( n, "w" )
    for p in sorted( m[n] ):
        with open( p ) as fh2:
            for line in fh2:
                print( line.strip( ), file=fh )
    fh.close( )
