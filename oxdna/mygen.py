#!/usr/bin/env python2

import base
import generators as gen
import numpy as np

s = gen.StrandGenerator()
strand = s.generate(4, double=False, dir=np.array([0,1,0]), start_position=np.array([0, 0, 1]))
syst = base.System([3, 3, 3])
syst.add_strand(strand)
syst.print_lorenzo_output ("prova.conf", "prova.top")
