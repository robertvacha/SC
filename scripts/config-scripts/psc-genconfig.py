#!/usr/bin/env python

#this program generates random coordinations and orientations (for spherocylinders)

import os
import sys
import math
import optparse
import commands
import string
import random
import gzip
import usefulmath


def gennewcoord(pbc):
    #generate random vector on 1 sphere
    orient=usefulmath.vec_random()
    orient=usefulmath.vec_normalize(orient)

    #generate position
    pos=[random.random()*float(pbc[0]),random.random()*float(pbc[1]),random.random()*float(pbc[2])]
    
    #generate patch (perpendicular to orientation of cylinder)
    patch=usefulmath.perp_vec(orient)
    patch=usefulmath.vec_normalize(patch)
    
    #return [pos,orient]
    return [pos,orient,patch]

def make(num,pbc,outfilename):
    data=[]
    
    for i in range(int(num)):
        data.append(gennewcoord(pbc))
    
    outstring=""
    for line in data:
        for part in line:
            for piece in part:
                outstring=outstring+str(piece)+" "
            outstring=outstring+" 0\n"
    
    output = open(outfilename, 'wb')
    #print(pbc)
    boxsize = ""
    for i in pbc:
        boxsize = boxsize + i + " "
    boxsize = boxsize + "\n"
    output.write(boxsize)
    try:
        output.write(outstring)
    finally:
        output.close()
            
        
    return 0

parser=optparse.OptionParser()
help="""Usage:
%program [options] 
pbc is set as string ("100 15 24")
"""
parser.set_usage(help)
parser.add_option(
    "-n",
    "--number",
    help="Number of cylinders",
    dest="num"
    )
parser.add_option(
    "-o",
    "--output",
    help="Set to which file you want to save data",
    dest="outfilename",
    default="config.init"
    )
parser.add_option(
    "--pbc",
    help="Set size of box you want to use \"x y z\"]",
    dest="pbc"
    )

(options,arguments)=parser.parse_args()
if not options.num:
    sys.exit("Error: You have to set number of cylinders (-n)")
if not options.pbc:
    sys.exit("Error: You have to set size of box (--pbc)")
make(options.num,options.pbc.split(),options.outfilename)
