#!/usr/bin/env python

import os
import sys
import re
import argparse

from zopy.utils import iter_rowdicts, tprint

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "bc_dist_init_details", help="" )
args = parser.parse_args( )

"""
Subject  202703
IndexA   0
DayA     126
SampleA  6382152
IndexB   1
DayB     154
SampleB  1107821
Delta    28
BCDist   0.6909557
"""

target = 12 * (365/12.0)
window = 60

diffs = {}
dists = {}
for R in iter_rowdicts( args.bc_dist_init_details ):
    s = R["Subject"]
    d = int( R["DayB"] )
    bc = float( R["BCDist"] )
    diff = abs( d - target )
    if diff < window:
        if s not in diffs or diffs[s] > diff:
            diffs[s] = diff
            dists[s] = bc
        
for k, v in dists.items( ):
    tprint( k, v )
    
