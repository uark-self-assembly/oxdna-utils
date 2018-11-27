#!/usr/bin/env python3

import sys
try:
    import numpy as np
except:
    print("error: no numpy installed. See requirements.txt", file=sys.stderr)
import base, readers


def print_usage():
    print("USAGE:")
    print("\t%s trajectory topology skip" % sys.argv[0])
    sys.exit(1)

try:
    traj = sys.argv[1]
    print(traj)
    topo = sys.argv[2]
    print(topo)
    skip = int(sys.argv[3])
    print(skip)
except:
    print_usage()

l = readers.LorenzoReader(traj, topo)
s = l.get_system(N_skip=skip)

base.Logger.log ('Writing to decimated.dat')
base.Logger.log ('skipping %d confs' % (skip))
out = open ('decimated.dat', 'w')
niter = 1
while s:
    base.Logger.log("Working on conf %i..." % niter, base.Logger.INFO)
    s = l.get_system(N_skip=skip)
    if s:
        s.print_lorenzo_output ('tmpdec.dat', 'tmpdec.top')
        inp = open ('tmpdec.dat', 'r');
        for line in inp.readlines():
            out.write (line)
        inp.close()
    niter += 1

out.close()
