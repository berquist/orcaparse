#!/usr/bin/env python

import numpy as np
import orca_parser
from copper_imidazole_analysis import CopperImidazoleAnalysis
import argparse

cia = CopperImidazoleAnalysis()

parser = argparse.ArgumentParser()
parser.add_argument(dest="namelist", nargs="+")
args = parser.parse_args()
namelist = args.namelist

for filename in namelist:
    orcafile = orca_parser.ORCAOutputParser(filename)

    gtensor, giso = orcafile.return_gtensor()
    g_para = gtensor[2]

    id_copper = cia.copper_id(orcafile)
    atom_copper = orcafile.molecule.atoms[id_copper]

    hyp = atom_copper.hyperfine

    A_para = hyp.atensor[2]
    A_FC = np.sum(hyp.afc)/len(hyp.afc)
    A_SOC_x = hyp.aso[0]
    A_SOC_y = hyp.aso[1]
    A_SOC_z = hyp.aso[2]
    A_x = hyp.atensor[0]
    A_y = hyp.atensor[1]
    A_z = hyp.atensor[2]
    A_iso = (A_x + A_y + A_z) / 3
    A_D_x = A_x - A_iso
    A_D_y = A_y - A_iso
    A_D_z = A_z - A_iso

    print "   name:", filename
    print " g_para:", g_para
    print " A_para:", A_para
    print "   A_FC:", A_FC
    print "  A_D_x:", A_D_x
    print "  A_D_y:", A_D_y
    print "  A_D_z:", A_D_z
    print "A_SOC_x:", A_SOC_x
    print "A_SOC_y:", A_SOC_y
    print "A_SOC_z:", A_SOC_z
    print "    A_x:", A_x
    print "    A_y:", A_y
    print "    A_z:", A_z
    print "  A_iso:", A_iso
    print "\n"
