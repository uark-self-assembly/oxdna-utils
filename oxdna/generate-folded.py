#!/usr/bin/env python2
"""
Standalone version of generate.py
Created for stable release end-users and machines without numpy support.
This file is only supported for stable release, and does not include more recent functionality.

(05 December 2012)
"""

import sys
import os
try:
    import numpy as np
except:
    import mynumpy as np


# return parts of a string
def partition(string_value, delimiter):
    if delimiter in string_value:
        split = string_value.split(delimiter, 1)
        return split[0], delimiter, split[1]
    else:
        return string_value, '', ''


def is_float(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


# every defined macro in model.h must be imported in this module
def import_model_constants():
    PI = np.pi
    model_file_path = os.path.join(os.path.dirname(__file__), '..', 'src', 'model.h')
    model_file = open(model_file_path)
    for line in model_file.readlines():
        line = partition(line.strip(), '//')[0].strip()

        macro = partition(line, '#define ')[2].strip().split(' ', 1)
        if len(macro) > 1:
            key, value = [x.strip() for x in macro]
            value = value.replace('f', '')
            globals()[key] = eval(value)

    model_file.close()


def add_strands(new_positions, new_a1s, new_a3s):
    overlap = False

    for i in xrange(len(positions)):
        p = positions[i]
        pa1 = a1s[i]

        for j in xrange(len(new_positions)):
            q = new_positions[j]
            qa1 = new_a1s[j]

            p_pos_back = p + pa1 * POS_BACK
            p_pos_base = p + pa1 * POS_BASE
            q_pos_back = q + qa1 * POS_BACK
            q_pos_base = q + qa1 * POS_BASE

            dr = p_pos_back - q_pos_back
            dr -= box * np.rint(dr / box)

            if np.dot(dr, dr) < RC2_BACK:
                overlap = True

            dr = p_pos_base - q_pos_base
            dr -= box * np.rint(dr / box)
            if np.dot(dr, dr) < RC2_BASE:
                overlap = True

            dr = p_pos_back - q_pos_base
            dr -= box * np.rint(dr / box)
            if np.dot(dr, dr) < RC2_BACK_BASE:
                overlap = True

            dr = p_pos_base - q_pos_back
            dr -= box * np.rint(dr / box)
            if np.dot(dr, dr) < RC2_BACK_BASE:
                overlap = True

            if overlap:
                return False

    if not overlap:
        for p in new_positions:
            positions.append(p)
        for p in new_a1s:
            a1s.append(p)
        for p in new_a3s:
            a3s.append(p)

    return True


def get_rotation_matrix(axis, anglest):
    # the argument anglest can be either an angle in radiants
    # (accepted types are float, int or np.float64 or np.float64)
    # or a tuple [angle, units] where angle a number and
    # units is a string. It tells the routine whether to use degrees,
    # radiants (the default) or base pairs turns
    if not isinstance(anglest, (np.float64, np.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ['degrees', 'deg', 'o']:
                angle = (np.pi / 180.) * (anglest[0])
            elif anglest[1] in ['bp']:
                angle = int(anglest[0]) * (np.pi / 180.) * 35.9
            else:
                angle = float(anglest[0])
        else:
            angle = float(anglest[0])
    else:
        angle = float(anglest)  # in degrees, I think

    axis = np.array(axis)
    axis /= np.sqrt(np.dot(axis, axis))

    ct = np.cos(angle)
    st = np.sin(angle)
    olc = 1. - ct
    x, y, z = axis

    return np.array([[olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
                    [olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
                    [olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct]])


class StrandGenerator(object):
    def generate(self,
                 bp,
                 sequence=None,
                 start_pos=np.array([0, 0, 0]),
                 direction=np.array([0, 0, 1]),
                 perpendicular=False,
                 double_strand=True,
                 rotation=0.):

        new_positions = []
        new_a1s = []
        new_a3s = []

        # we need a numpy array for these
        start_pos = np.array(start_pos, dtype=float)
        direction = np.array(direction, dtype=float)
        if sequence is None:
            sequence = np.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += np.random.randint(0, 4, n)
            print >> sys.stderr, 'sequence is too short, adding %d random bases' % n

        # create the sequence of the second strand as made of complementary bases
        sequence2 = [(3 - s) for s in sequence]
        sequence2.reverse()

        # we need to find a vector orthogonal to dir
        dir_norm = np.sqrt(np.dot(direction, direction))
        if dir_norm < 1e-10:
            print >> sys.stderr, 'direction must be a valid vector, defaulting to (0, 0, 1)'
            direction = np.array([0, 0, 1])
        else:
            direction /= dir_norm

        if perpendicular is None or perpendicular is False:
            v1 = np.random.random_sample(3)
            v1 -= direction * (np.dot(direction, v1))
            v1 /= np.sqrt(sum(v1*v1))
        else:
            v1 = perpendicular

        # and we need to generate a rotational matrix
        R0 = get_rotation_matrix(direction, rotation)
        R = get_rotation_matrix(direction, [1, 'bp'])

        a1 = v1
        a1 = np.dot(R0, a1)
        rb = np.array(start_pos)
        a3 = direction
        for i in range(bp):
            rcdm = rb - CM_CENTER_DS * a1
            new_positions.append(rcdm)
            new_a1s.append(a1)
            new_a3s.append(a3)

            if i != bp-1:
                a1 = np.dot(R, a1)
                rb += a3 * BASE_BASE

        if double_strand:
            a1 = -a1
            a3 = -direction
            R = R.transpose()
            for i in range(bp):
                rcdm = rb - CM_CENTER_DS * a1
                new_positions.append(rcdm)
                new_a1s.append(a1)
                new_a3s.append(a3)

                a1 = np.dot(R, a1)
                rb += a3 * BASE_BASE

        assert len(new_positions) > 0

        return [new_positions, new_a1s, new_a3s]


def read_strands(filename):
    """
    The main() function for this script
    Reads a text file with the following format:
    - Each line contains the sequence for a single strand (A,C,T,G)
    - Lines begining in DOUBLE produce double stranded DNA

    Ex: Two ssDNA (single stranded DNA)
    ATATATA
    GCGCGCG

    Ex: Two strands, one double stranded, the other single stranded.
    DOUBLE AGGGCT
    CCTGTA

    """
    # we take sequences from a file; each line one strand if it's a
    # double strand, there must be a DOUBLE in front of it

    try:
        infile = open(filename)
    except:
        print >> sys.stderr, 'Could not open file ', filename, 'Aborting now'
        sys.exit(2)

    # get the number of strands and nucleotides
    stand_count = 0
    nucleotide_count = 0

    lines = infile.readlines()
    for line in lines:
        line = line.upper().strip()
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            print >> sys.stderr, '## Found duplex of %i bases' % len(line)
            nucleotide_count += 2 * len(line)
            stand_count += 2
        else:
            line = line.split()[0]
            print >> sys.stderr, '## Found single strand of %i bases' % (len(line))
            nucleotide_count += len(line)
            stand_count += 1

    infile.seek(0)

    print >> sys.stderr, '## nstrands, nnucl = ', stand_count, nucleotide_count

    # here we generate the topology file
    try:
        topology_file = open('generated.top', 'w')
    except:
        print >> sys.stderr, 'Could not open generated.top for writing. Aborting'
        sys.exit(4)

    print >> topology_file, nucleotide_count, stand_count

    myns = 1
    mynn = 0
    for line in lines:
        line = line.upper().strip()
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1].upper()
            print >> topology_file, myns, line[0], -1, mynn + 1
            mynn += 1
            for i in xrange(1, len(line) - 1):
                print >> topology_file, myns, line[i], mynn - 1, mynn + 1
                mynn += 1
            print >> topology_file, myns, line[-1], mynn - 1, -1
            mynn += 1
            myns += 1

            # get the compl sequence in numbers
            seq = [3 - base_to_number[x] for x in line]
            # put it back in letters
            line = [number_to_base[x] for x in seq[::-1]]

            print >> topology_file, myns, line[0], -1, mynn + 1
            mynn += 1
            for i in xrange(1, len(line) - 1):
                print >> topology_file, myns, line[i], mynn - 1, mynn + 1
                mynn += 1
            print >> topology_file, myns, line[-1], mynn - 1, -1
            mynn += 1
            myns += 1
        else:
            line = line.split()[0].upper()
            print >> topology_file, myns, line[0], -1, mynn + 1
            mynn += 1
            for i in xrange(1, len(line) - 1):
                print >> topology_file, myns, line[i], mynn - 1, mynn + 1
                mynn += 1
            print >> topology_file, myns, line[-1], mynn - 1, -1
            mynn += 1
            myns += 1
    topology_file.close()
    infile.seek(0)

    # generate the strands
    strand_generator = StrandGenerator()
    lines = infile.readlines()
    lines_count = len(lines)
    i = 1

    for line in lines:
        line = line.upper().strip()

        # skip empty lines
        if len(line) == 0:
            continue
        if line[:6] == 'DOUBLE':
            line = line.split()[1]
            seq = [base_to_number[x] for x in line]
            length = len(line)
            print >> sys.stderr, '## Adding duplex of %i bases' % length
            cdm = np.random.random_sample(3) * box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            new_positions, new_a1s, new_a3s = strand_generator.generate(len(line), sequence=seq, direction=axis, start_pos=cdm, double_strand=True)
            while not add_strands(new_positions, new_a1s, new_a3s):
                cdm = np.random.random_sample(3) * box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                new_positions, new_a1s, new_a3s = strand_generator.generate(len(line), sequence=seq, direction=axis, start_pos=cdm, double_strand=True)
            print >> sys.stderr, '##  done line %i / %i, now at %i/%i' % (i, lines_count, len(positions), nucleotide_count)
        else:
            seq = [base_to_number[x] for x in line]
            cdm = np.random.random_sample(3) * box
            axis = np.random.random_sample(3)
            axis /= np.sqrt(np.dot(axis, axis))
            print >> sys.stderr, '## Adding single strand of %i bases' % (len(line))
            new_positions, new_a1s, new_a3s = strand_generator.generate(len(line), sequence=seq, direction=axis, start_pos=cdm, double_strand=False)
            while not add_strands(new_positions, new_a1s, new_a3s):
                cdm = np.random.random_sample(3) * box
                axis = np.random.random_sample(3)
                axis /= np.sqrt(np.dot(axis, axis))
                new_positions, new_a1s, new_a3s = strand_generator.generate(len(line), sequence=seq, direction=axis, start_pos=cdm, double_strand=False)
            print >> sys.stderr, '##  done line %i / %i, now at %i/%i' % (i, lines_count, len(positions), nucleotide_count)

        i += 1

    if not len(positions) == nucleotide_count:
        print len(positions), nucleotide_count
        raise AssertionError

    # here we generate the configuration file (coordinates)
    try:
        outfile = open('generated.dat', 'w')
    except:
        print >> sys.stderr, 'Could not open generated.dat for writing.  Aborting'
        sys.exit(5)

    print >> outfile, 't = 0'
    print >> outfile, 'b = ', box_side, box_side, box_side
    print >> outfile, 'E = 0. 0. 0.'
    for i in xrange(nucleotide_count):
        print >> outfile, positions[i][0], positions[i][1], positions[i][2],
        print >> outfile, a1s[i][0], a1s[i][1], a1s[i][2],
        print >> outfile, a3s[i][0], a3s[i][1], a3s[i][2],
        print >> outfile, 0., 0., 0., 0., 0., 0.  # v and L

    outfile.close()
    print >> sys.stderr, '## ALL DONE. just generated "generated.dat" and "generated.top"'


##
#  MAIN
##
if __name__ == '__main__':
    try:
        box_side = float(sys.argv[1])
        input_file_path = sys.argv[2]
    except:
        print >> sys.stderr, 'Usage: %s <%s> <%s>' % (sys.argv[0], 'box size', 'file with sequences')
        sys.exit(1)

    try:
        input_file = open(input_file_path, 'r')
        input_file.close()
    except:
        print >> sys.stderr, 'Could not open file "%s" for reading. Aborting' % input_file_path
        sys.exit(2)

    box = np.array([box_side, box_side, box_side])

    import_model_constants()

    CM_CENTER_DS = POS_BASE + 0.2
    BASE_BASE = 0.3897628551303122

    RC2_BACK = EXCL_RC1**2
    RC2_BASE = EXCL_RC2**2
    RC2_BACK_BASE = EXCL_RC3**2

    number_to_base = {0: 'A', 1: 'G', 2: 'C', 3: 'T'}
    base_to_number = {'A': 0, 'a': 0, 'G': 1, 'g': 1,
                      'C': 2, 'c': 2, 'T': 3, 't': 3}

    positions = []
    a1s = []
    a3s = []
    newpositions = []
    newa1s = []
    newa3s = []

    read_strands(input_file_path)
