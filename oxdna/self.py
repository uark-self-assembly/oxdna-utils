#!/usr/bin/env python3

import base
import readers
try:
    import numpy as np
except:
    import mynumpy as np
import os.path
import sys
import external_forces as forces

if len(sys.argv) < 2:
    base.Logger.log("Usage is %s <file> [<block>]" % sys.argv[0], base.Logger.CRITICAL)
    sys.exit()

try:
    block = int(sys.argv[2])
except:
    print("# Using default block of 100", file=sys.stdout)
    block = 100

inp = open (sys.argv[1], 'r')

linea = inp.readline ()
nrows, avge, avgee = 0, 0., 0.
es = []
while linea:
    words = linea.split ()
    t, e = float (words[0]), float (words[1])
    es.append(e)
    avgee += e * e
    avge += e
    linea = inp.readline ()

inp.close ()

nrows = len (es)
avgee /= float (nrows)
avge /= float (nrows)
print("# nrows, <x>, <x^2>, sigma^2 :", nrows, avge, avgee, (avgee - avge * avge), file=sys.stderr)

nblocks = 0
vtv0 = [0. for i in range (block)]
for i in range (0, nrows - block):
    for j in range (block):
        vtv0[j] += (es[i + j] - avge)*(es[i] - avge)
    nblocks += 1

for i in range (len (vtv0)):
    print(i, vtv0[i] / float (nblocks), vtv0[i] / float (nblocks) / (avgee - avge * avge))
