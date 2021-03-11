#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser( )
parser.add_argument( '-j' , metavar='job', help="job name" )
parser.add_argument( '-t' , metavar='time', help="runtime on cluster" )
parser.add_argument( '-m' , metavar='microbiome', help="microbiome data" )
parser.add_argument( '-f' , metavar='feature', help="microbiome feature" )
parser.add_argument( '-g' , metavar='snps', help="genetics data" )
parser.add_argument( '-d' , metavar='metadata', help="metadata" )
parser.add_argument( '-c' , metavar='cpus', help="number of cpus" )
parser.add_argument( '-o' , metavar='output', help="model output" )
args = parser.parse_args( )

shell = """#!/bin/sh
#SBATCH -J {job}
#SBATCH -o {job}.stdout
#SBATCH -e {job}.stderror
#SBATCH -p serial_requeue
#SBATCH -n {cpus}
#SBATCH --mem-per-cpu 4000
#SBATCH -t {time}

module load Anaconda2
source activate teddy_genetics 

run_models_snps.R -m {microbiome} -f {feature} -g {snps} -d {metadata} -c {cpus} -o {output}

gzip {output}

#END"""

shell = shell.format( job=args.j, time=args.t, microbiome=args.m, feature=args.f, snps=args.g, metadata=args.d, cpus=args.c, output=args.o )

print( shell )
