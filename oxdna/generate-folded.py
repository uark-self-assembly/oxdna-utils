#!/usr/bin/env python3
"""
Standalone version of generate.py
Created for stable release end-users and machines without numpy support.
This file is only supported for stable release, and does not include more recent functionality.
"""

import sys
import os
import numpy
import argparse

from typing import List, Dict

constants = {'PI': numpy.pi}


class Strand:
    _complement_table: Dict[int, int] = {
        ord('A'): ord('T'), ord('T'): ord('A'), ord('G'): ord('C'), ord('C'): ord('G'),
        ord('a'): ord('T'), ord('t'): ord('A'), ord('g'): ord('C'), ord('c'): ord('G'),
    }

    _base_to_number_table = {
        'A': 0, 'a': 0, 'G': 1, 'g': 1,
        'C': 2, 'c': 2, 'T': 3, 't': 3
    }

    def __init__(self, sequence: str='', double: bool=False, line: str=None):
        """
        Initializes a new Strand object. To properly initialize a Strand, either supply a sequence as a string,
        or supply a "line" from a file.

        Example sequences:

        ATCGCCA, actggcaa, GCtAccG

        Example lines:

        ATCGCCA
        DOUBLE ATCGCCA
        atcgACcctG
        douBLe ttcgAAc

        :param sequence: a `str` object, with characters containing only 'a', 'g', 'c', 't', either upper or lower case
        :param double: a `boolean`, whether or not this is a double strand
        :param line: a `str` to parse into a `sequence` and `double`. If supplied, `sequence` and `double` are ignored
        """
        if line:
            split = line.split()
            if split[0].upper() == 'DOUBLE':
                self.double: bool = True
                self.sequence: str = split[1].upper()
            else:
                self.double: bool = False
                self.sequence: str = split[0].upper()
        else:
            self.double: bool = double
            self.sequence: str = sequence.upper()

    def complement(self) -> str:
        return self.sequence.translate(Strand._complement_table)

    def nucleotide_count(self) -> int:
        return len(self.sequence) * self.strand_count()

    def strand_count(self) -> int:
        return 2 if self.double else 1

    def to_integer_array(self) -> List[int]:
        return [Strand._base_to_number_table[i] for i in self.sequence]


class System:
    def __init__(self, strands: List[Strand]=None):
        if strands is None:
            strands = []
        self.strands: List[Strand] = strands
        self.nucleotide_count: int = sum([strand.nucleotide_count() for strand in self.strands])
        self.strand_count: int = sum([strand.strand_count() for strand in self.strands])
        self.positions: List[numpy.ndarray] = []
        self.a1s: List[numpy.ndarray] = []
        self.a3s: List[numpy.ndarray] = []

    def add_strand(self,
                   positions: List[numpy.ndarray],
                   a1s: List[numpy.ndarray],
                   a3s: List[numpy.ndarray]) -> bool:
        overlap = False

        for i in range(len(self.positions)):
            p = self.positions[i]
            pa1 = self.a1s[i]

            for j in range(len(positions)):
                q = positions[j]
                qa1 = a1s[j]

                p_pos_back = p + pa1 * constants['POS_BACK']
                p_pos_base = p + pa1 * constants['POS_BASE']
                q_pos_back = q + qa1 * constants['POS_BACK']
                q_pos_base = q + qa1 * constants['POS_BASE']

                dr = p_pos_back - q_pos_back
                dr -= box * numpy.rint(dr / box)

                if numpy.dot(dr, dr) < constants['RC2_BACK']:
                    overlap = True

                dr = p_pos_base - q_pos_base
                dr -= box * numpy.rint(dr / box)
                if numpy.dot(dr, dr) < constants['RC2_BASE']:
                    overlap = True

                dr = p_pos_back - q_pos_base
                dr -= box * numpy.rint(dr / box)
                if numpy.dot(dr, dr) < constants['RC2_BACK_BASE']:
                    overlap = True

                dr = p_pos_base - q_pos_back
                dr -= box * numpy.rint(dr / box)
                if numpy.dot(dr, dr) < constants['RC2_BACK_BASE']:
                    overlap = True

                if overlap:
                    return False

        if not overlap:
            self.positions += positions
            self.a1s += a1s
            self.a3s += a3s

        return True


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

    model_file_path = os.path.join(os.path.dirname(__file__), '..', 'src', 'model.h')
    model_file = open(model_file_path)
    for line in model_file.readlines():
        line = partition(line.strip(), '//')[0].strip()

        macro = partition(line, '#define ')[2].strip().split(' ', 1)
        if len(macro) > 1:
            key, value = [x.strip() for x in macro]
            value = value.replace('f', '')
            constants[key] = eval(value, constants)

    model_file.close()


def get_rotation_matrix(axis, anglest):
    # the argument anglest can be either an angle in radiants
    # (accepted types are float, int or np.float64 or np.float64)
    # or a tuple [angle, units] where angle a number and
    # units is a string. It tells the routine whether to use degrees,
    # radiants (the default) or base pairs turns
    if not isinstance(anglest, (numpy.float64, numpy.float32, float, int)):
        if len(anglest) > 1:
            if anglest[1] in ['degrees', 'deg', 'o']:
                angle = (numpy.pi / 180.) * (anglest[0])
            elif anglest[1] in ['bp']:
                angle = int(anglest[0]) * (numpy.pi / 180.) * 35.9
            else:
                angle = float(anglest[0])
        else:
            angle = float(anglest[0])
    else:
        angle = float(anglest)  # in degrees, I think

    axis = numpy.array(axis)
    axis /= numpy.sqrt(numpy.dot(axis, axis))

    ct = numpy.cos(angle)
    st = numpy.sin(angle)
    olc = 1. - ct
    x, y, z = axis

    return numpy.array([[olc*x*x+ct, olc*x*y-st*z, olc*x*z+st*y],
                       [olc*x*y+st*z, olc*y*y+ct, olc*y*z-st*x],
                       [olc*x*z-st*y, olc*y*z+st*x, olc*z*z+ct]])


class StrandGenerator(object):
    def generate(self,
                 bp,
                 sequence=None,
                 start_pos=numpy.array([0, 0, 0]),
                 direction=numpy.array([0, 0, 1]),
                 perpendicular=False,
                 double_strand=True,
                 rotation=0.):

        new_positions = []
        new_a1s = []
        new_a3s = []

        # we need a numpy array for these
        start_pos = numpy.array(start_pos, dtype=float)
        direction = numpy.array(direction, dtype=float)
        if sequence is None:
            sequence = numpy.random.randint(0, 4, bp)
        elif len(sequence) != bp:
            n = bp - len(sequence)
            sequence += numpy.random.randint(0, 4, n)
            print('sequence is too short, adding {} random bases'.format(n), file=sys.stderr)

        # create the sequence of the second strand as made of complementary bases
        sequence2 = [(3 - s) for s in sequence]
        sequence2.reverse()

        # we need to find a vector orthogonal to dir
        dir_norm = numpy.sqrt(numpy.dot(direction, direction))
        if dir_norm < 1e-10:
            print('direction must be a valid vector, defaulting to (0, 0, 1)', file=sys.stderr)
            direction = numpy.array([0, 0, 1])
        else:
            direction /= dir_norm

        if perpendicular is None or perpendicular is False:
            v1 = numpy.random.random_sample(3)
            v1 -= direction * (numpy.dot(direction, v1))
            v1 /= numpy.sqrt(sum(v1*v1))
        else:
            v1 = perpendicular

        # and we need to generate a rotational matrix
        R0 = get_rotation_matrix(direction, rotation)
        R = get_rotation_matrix(direction, [1, 'bp'])

        a1 = v1
        a1 = numpy.dot(R0, a1)
        rb = numpy.array(start_pos)
        a3 = direction
        for i in range(bp):
            rcdm = rb - constants['CM_CENTER_DS'] * a1
            new_positions.append(rcdm)
            new_a1s.append(a1)
            new_a3s.append(a3)

            if i != bp-1:
                a1 = numpy.dot(R, a1)
                rb += a3 * constants['BASE_BASE']

        if double_strand:
            a1 = -a1
            a3 = -direction
            R = R.transpose()
            for i in range(bp):
                rcdm = rb - constants['CM_CENTER_DS'] * a1
                new_positions.append(rcdm)
                new_a1s.append(a1)
                new_a3s.append(a3)

                a1 = numpy.dot(R, a1)
                rb += a3 * constants['BASE_BASE']

        assert len(new_positions) > 0

        return [new_positions, new_a1s, new_a3s]


def parse_strands(file_path):
    sequence_file = open(file_path, 'r')
    lines = sequence_file.readlines()
    return [Strand(line=line) for line in lines]


def generate_topology(system: System):
    print('strand_count, nucleotide_count = {}, {}'.format(system.strand_count, system.nucleotide_count), file=sys.stderr)

    # here we generate the topology file
    try:
        topology_file = open('generated.top', 'w')
    except:
        print('Could not open "generated.top" for writing. Aborting', file=sys.stderr)
        sys.exit(4)

    print(system.nucleotide_count, system.strand_count, file=topology_file)

    strand_index = 1
    nucleotide_index = 0

    for strand in system.strands:
        if strand.nucleotide_count() == 0:
            continue

        double = strand.double
        sequence = strand.sequence

        print(strand_index, sequence[0], -1, nucleotide_index + 1, file=topology_file)

        nucleotide_index += 1

        for i in range(1, len(sequence) - 1):
            print(strand_index, sequence[i], nucleotide_index - 1, nucleotide_index + 1, file=topology_file)
            nucleotide_index += 1

        print(strand_index, sequence[-1], nucleotide_index - 1, -1, file=topology_file)

        nucleotide_index += 1
        strand_index += 1

        if double:
            sequence = strand.complement()

            print(strand_index, sequence[0], -1, nucleotide_index + 1, file=topology_file)

            nucleotide_index += 1

            for i in range(1, len(sequence) - 1):
                print(strand_index, sequence[i], nucleotide_index - 1, nucleotide_index + 1, file=topology_file)
                nucleotide_index += 1

            print(strand_index, sequence[-1], nucleotide_index - 1, -1, file=topology_file)

            nucleotide_index += 1
            strand_index += 1

    topology_file.close()


def generate_dat(system: System):
    # generate the strands
    strand_generator = StrandGenerator()
    lines_count = len(system.strands)
    i = 1

    for strand in system.strands:
        double = strand.double
        sequence = strand.to_integer_array()
        # skip empty lines
        if len(sequence) == 0:
            continue

        print(
            '## Adding {} of {} bases'.format('duplex' if double else 'single strand', len(sequence)),
            file=sys.stderr)

        strands_added = False
        while not strands_added:
            cdm = numpy.random.random_sample(3) * box
            axis = numpy.random.random_sample(3)
            axis /= numpy.sqrt(numpy.dot(axis, axis))
            positions, a1s, a3s = strand_generator.generate(len(sequence), sequence=sequence, direction=axis,
                                                            start_pos=cdm, double_strand=double)

            strands_added = system.add_strand(positions, a1s, a3s)

        print('##  done line {} / {}, now at {}/{}'.format(i, lines_count, len(system.positions), system.nucleotide_count),
              file=sys.stderr)

        i += 1

    if len(system.positions) != system.nucleotide_count:
        print(len(system.positions), system.nucleotide_count)
        raise AssertionError

    # here we write the configuration file (coordinates)
    try:
        dat_file = open('generated.dat', 'w')
    except:
        print('Could not open generated.dat for writing. Aborting', file=sys.stderr)
        sys.exit(5)

    print('t = 0', file=dat_file)
    print('b = {} {} {}'.format(box_size, box_size, box_size), file=dat_file)
    print('E = 0. 0. 0.', file=dat_file)

    for i in range(system.nucleotide_count):
        print(
            system.positions[i][0], system.positions[i][1], system.positions[i][2],
            system.a1s[i][0], system.a1s[i][1], system.a1s[i][2],
            system.a3s[i][0], system.a3s[i][1], system.a3s[i][2],
            0., 0., 0.,
            0., 0., 0.,
            file=dat_file)

    dat_file.close()


##
#  MAIN
##
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate an oxDNA topology using a sequence file.')
    parser.add_argument(
        'input_file_path', metavar='SEQUENCE_FILE', type=str, help='A file with DNA sequences on separate lines')
    parser.add_argument(
        '-b', '--box-size', type=int, dest='box_size', default=50, help='The box size to use. Default is 50')

    args = parser.parse_args()

    input_file_path = args.input_file_path
    box_size = args.box_size

    try:
        input_file = open(input_file_path, 'r')
        input_file.close()
    except:
        print('Could not open file {} for reading. Aborting.'.format(input_file_path), file=sys.stderr)
        sys.exit(2)

    box = numpy.array([box_size, box_size, box_size])

    import_model_constants()

    constants['CM_CENTER_DS'] = constants['POS_BASE'] + 0.2
    constants['BASE_BASE'] = 0.3897628551303122

    constants['RC2_BACK'] = constants['EXCL_RC1']**2
    constants['RC2_BASE'] = constants['EXCL_RC2']**2
    constants['RC2_BACK_BASE'] = constants['EXCL_RC3']**2

    parsed_strands = parse_strands(input_file_path)
    oxdna_system = System(strands=parsed_strands)
    generate_topology(oxdna_system)
    generate_dat(oxdna_system)

    print('## ALL DONE. Generated "generated.dat" and "generated.top"', file=sys.stderr)
