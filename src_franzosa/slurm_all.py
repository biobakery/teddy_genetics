#!/usr/bin/env python

import os
import sys
import re
import time
import argparse
from collections import Counter
import subprocess

# argument parsing (python argparse)
parser = argparse.ArgumentParser( )
parser.add_argument( "--features", nargs="+", required=True )
parser.add_argument( "--genetics", nargs="+", required=True )
parser.add_argument( "--execute", action="store_true" )
args = parser.parse_args( )

# ---------------------------------------------------------------
# read sacct data
# ---------------------------------------------------------------

output = subprocess.check_output( "sacct -u franzosa --format=jobid,jobname,state -P".split( ) )
output = str( output, "utf-8" ).split( "\n" )[1:-1]
jobstate = {}
for line in output:
    jobid, jobname, status = line.split( "|" )
    jobid = int( jobid.split( "." )[0] )
    status = status.split( )[0]
    jobstate.setdefault( jobname, [] ).append( [jobid, status] )
for j, pairs in jobstate.items( ):
    jobstate[j] = sorted( pairs )[-1][1]

# ---------------------------------------------------------------
# return to main flow 
# ---------------------------------------------------------------

cc = Counter( )

items = []
for f in args.features:
    for g in args.genetics:
        name = ".".join( [
                os.path.split( f )[1].replace( ".", "_" ),
                os.path.split( g )[1].replace( ".", "_" ),
                ] )
        items.append( [f, g, name] )

for f, g, name in items:
    print( f, g, name, end=" " )
    slurm = "slurm/{}.slurm".format( name )
    final = "slurm/{}.results.txt".format( name )
    state = jobstate.get( name, None )
    if state in ["RUNNING", "PENDING"]:
        label = state
    elif os.path.exists( final ) and state not in ["COMPLETED", None]:
        label = "FINISHED: {}".format( state )
    elif os.path.exists( final ):
        label = "FINISHED"
    else:
        os.system( "python src/slurm.py {} {} {}".format( f, g, name ) )
        label = "PRIMED"
        if args.execute:
            os.system( "sbatch slurm/{}.slurm".format( name ) )
            label = "DISPATCHED"
            time.sleep( 1 )
        if state is not None:
            label = "{} after {}".format( label, state )
    print( "<--", label )
    cc[label] += 1
    cc["[TOTAL]"] += 1

for label in sorted( cc ):
    print( label, cc[label] )
