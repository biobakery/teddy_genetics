#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser( )
parser.add_argument( '-i' , metavar='INPUT', help="model input" )
parser.add_argument( '-o' , metavar='OUTPUT', help="model output" )
parser.add_argument( '-c' , metavar='CPUS', help="number of cpus" )
parser.add_argument( '-t' , metavar='TIME', help="runtime on cluster" )
args = parser.parse_args( )

shell = """#!/bin/sh
#SBATCH -J {JOB}
#SBATCH -o {JOB}.stdout
#SBATCH -e {JOB}.stderror
#SBATCH -p serial_requeue
#SBATCH -n {CPUS}
#SBATCH --mem-per-cpu 4000
#SBATCH -t {TIME}

source activate teddy_genetics 

run_models.R -i {INPUT}.model_input.tsv -c {CPUS} -o {OUTPUT}

#END"""

input_file = args.i.replace(".model_input.tsv", "")

output_file = args.o

cpus = args.c

job_name = input_file.replace("snps.", "")

job_time = args.t

shell = shell.format( INPUT=input_file, OUTPUT=output_file, CPUS=cpus, JOB=job_name, TIME=job_time )

print( shell )
