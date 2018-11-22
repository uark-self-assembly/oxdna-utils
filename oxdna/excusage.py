#!/usr/bin/env python3

# This is EXCUSAGE, as in EXclusively Complementary Useful StrAnd GEnerator


import sys
import random

LENGTH = 6
NSEQS = 12
MAXTRIALS = 5000
ALLOWED = 2

def get_nbp (u, v):
    n = 0
    #for off in range(-LENGTH + 1, LENGTH):
    #    for i in range(0, LENGTH):
    #        if i > 0 and i+off > 0 and i+off < LENGTH:
    #            if u[i] + v[i + off] == 3:
    #                n += 1
    for i in range(0, LENGTH):
        if u[i] + v[i] == 3:
            n += 1
    return n

seqs, cseqs = [], []
trials = 0
while len(seqs) < NSEQS:
    prova = [random.randint (0, 3) for i in range(LENGTH)]
    check = True
    for s in seqs + cseqs:
        if get_nbp(prova, s) > ALLOWED:
            check = False
    if check:
        seqs.append(prova)
        cseqs.append([3 - x for x in prova])
        print("ADDED sequence ", prova, len(seqs), file=sys.stderr)
    else:
        pass
        #print >> sys.stderr,"REFUSED sequence", len(seqs)
    trials = trials + 1
    if trials > MAXTRIALS * NSEQS:
        seqs, cseqs = [], []
        trials = 0
        print("TOO MANY TRIALS, RESETTING", file=sys.stderr)


print("DONE", file=sys.stderr)
import base

for i in range(NSEQS):
    str = ""
    for x in seqs[i]:
        str += base.number_to_base[x]
    print(str)
    str = ""
    for x in cseqs[i]:
        str += base.number_to_base[x]
    print(str)
