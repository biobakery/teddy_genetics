#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import argparse

from zopy.utils import path2name

parser = argparse.ArgumentParser( )
parser.add_argument( "features" )
parser.add_argument( "genetics" )
parser.add_argument( "name" )
args = parser.parse_args( )

# ---------------------------------------------------------------
# SLURM
# ---------------------------------------------------------------

command = """#!/usr/bin/env bash

#SBATCH -J {NAME}
#SBATCH -o slurm/{NAME}.stdout
#SBATCH -e slurm/{NAME}.stderr
#SBATCH -p hii02
#SBATCH -t 0-03:00
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem=4000

# use <-p hii02>    for production
# use <-p hii-test> for testing

python src/model_rpy2.py \\
    {FEATURES} \\
    {GENETICS} \\
    slurm/{NAME}.results.txt \\

# END"""

# ---------------------------------------------------------------
# python
# ---------------------------------------------------------------

command = command.format( 
    FEATURES=args.features,
    GENETICS=args.genetics,
    NAME=args.name,
)

with open( os.path.join( "slurm", args.name + ".slurm" ), "w" ) as fh:
    print( command, file=fh )
