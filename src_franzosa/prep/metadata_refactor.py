#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import re
import argparse

import zopy.utils as zu
from zopy.table2 import nesteddict2table

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "config" )
parser.add_argument( "primary" )
parser.add_argument( "--recoding" )
args = parser.parse_args( )

# ---------------------------------------------------------------
# utils
# ---------------------------------------------------------------

def unformat( string ):
    if string.startswith( "b'" ):
        string = string[2:-1]
    if string.endswith( ".0" ):
        string = string[0:-2]
    return string

# ---------------------------------------------------------------
# load config
# ---------------------------------------------------------------

data = {}
for line in zu.iter_lines( args.config ):
    items = line.rstrip( ).split( "\t" )
    if len( items ) == 1:
        inner = data.setdefault( items[0], {} )
    elif len( items ) == 2:
        items = items[1].split( " >>> " )
        a = items[0]
        b = items[1] if len( items ) == 2 else a
        inner[a] = b

# ---------------------------------------------------------------
# load dictionary mappings
# ---------------------------------------------------------------

coding = {}
for p in data:
    inner = coding.setdefault( p, {} )
    for row in zu.iter_rows( p ):
        row = [unformat( k ) for k in row]
        if row[0] in data[p] and "=" in row[2]:
            temp = {}
            for match in re.finditer( "(\d+)\s*=\s*(.*?)(,|;|\"|$)", row[2] ):
                hits = match.groups( )
                temp[hits[0]] = hits[1]
            inner[row[0]] = temp

coding_path = args.config + ".coding"
with open( coding_path, "w" ) as fh:
    for p in sorted( coding ):
        for term in sorted( coding[p] ):
            for code in sorted( coding[p][term] ):
                recode = coding[p][term][code]
                zu.tprint( p, term, code, recode, file=fh )

# ---------------------------------------------------------------
# load remapping
# ---------------------------------------------------------------

if args.recoding is None:
    args.recoding = coding_path

recoding = {}
for row in zu.iter_rows( args.recoding ):
    p, term, code, recode = row[0:4]
    recode = recode.split( " >>> " )[-1]
    inner = recoding.setdefault( p, {} ).setdefault( term, {} )
    inner[code] = recode

# ---------------------------------------------------------------
# load values
# ---------------------------------------------------------------

"""
b'Mask Id'  b'Race'  b'Race Ethnicity'  b'Sex'     
925749.0    6.0      5.0                b'Male'    
244537.0    6.0      5.0                b'Male'    
733269.0    1.0      2.0                b'Female'  
965870.0    1.0      2.0                b'Male'
"""

output = {}
for p in data:
    headers = None
    p2 = p.replace( "_dictionary", "" )
    rows = []
    for row in zu.iter_rows( p2 ):
        row = [unformat( k ) for k in row]
        if headers is None:
            headers = [k.lower( ).replace( " ", "_" ) for k in row]
        else:
            row = {a:b for a, b in zip( headers, row )}
            rows.append( row )
    for row in rows:
        inner = output.setdefault( row[args.primary], {} )
        for term, value in row.items( ):
            if term in data[p] and term != args.primary:
                term2 = data[p][term]
                if p in recoding and term in recoding[p]:
                    value = recoding[p][term].get( value, "#N/A" )
                inner[term2] = value

# ---------------------------------------------------------------
# finalize
# ---------------------------------------------------------------
            
T = nesteddict2table( output )
T.data[0][0] = args.primary
T.transpose( )
T.dump( )
