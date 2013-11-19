#!/usr/bin/env python2

import orca_parser
from copper_imidazole_analysis import CopperImidazoleAnalysis
import argparse

cia = CopperImidazoleAnalysis()

parser = argparse.ArgumentParser(description="Given pathnames of ORCA output files, pretty-print the electronic g-tensor and the copper/far nitrogen hyperfine tensors.")

parser.add_argument(dest="namelist", metavar="<orca filename>", nargs="+", type=str, default=None, help="ORCA input or output files.")

args = parser.parse_args()
namelist = args.namelist

s = "{:>34s} {:>31s} {:>3s} {:>28s} {:>3s} {:<s}"
print s.format("g-tensor", "a_copper", "id", "a_nitrogen_far", "id", "name")

for name in namelist:
    orcafile = orca_parser.ORCAOutputParser(name)
    gtensor, giso = orcafile.return_gtensor()
    id_cu = cia.copper_id(orcafile)
    id_far = cia.nitrogen_far_id(orcafile)
    atensor_cu = cia.hyperfine(orcafile, id_cu)
    atensor_far = cia.hyperfine(orcafile, id_far)

    s = "{:>34s} {:>28s} {:>3d} {:>25s} {:>3d} {:<s}"
    print s.format(gtensor, atensor_cu, id_cu, atensor_far, id_far, name)
