#! /usr/bin/env python

import os
import sys
import shutil
import re

for p in os.listdir( sys.argv[1] ):
    if p.endswith( ".tsv" ):
        name = os.path.split( p )[1]
        name = name.split( "_2019" )[0].lower( )
        shutil.copy( os.path.join( sys.argv[1], p ), name + ".tsv" )
