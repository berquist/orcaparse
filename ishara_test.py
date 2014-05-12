#!/usr/bin/env python

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
    g_perp = (gtensor[0] + gtensor[1])/2

    id_copper = cia.copper_id(orcafile)
    ids_nitrogen = []

    atom_copper = orcafile.molecule.atoms[id_copper]
    atoms_nitrogen = []

    # if we are within the distance cutoff, we must be directly bonded, and do not
    # contribute to the ESEEM results, so leave those out
    cutoff = 3.0
    for id_n in cia.nitrogen_ids(orcafile):
        if orcafile.pair_distance(id_copper, id_n) > cutoff:
            ids_nitrogen.append(id_n)
            atoms_nitrogen.append(orcafile.molecule.atoms[id_n])

    A_para = atom_copper.hyperfine.atensor[2]
    A_perp = (atom_copper.hyperfine.atensor[0] + atom_copper.hyperfine.atensor[1])/2

    # need to average A_iso, T_dip, kappa, and eta
    vals_A_iso = []
    vals_T_dip = []
    vals_kappa = []
    vals_eta = []
    for id_n in ids_nitrogen:
        hyp = orcafile.molecule.atoms[id_n].hyperfine
        efg = orcafile.molecule.atoms[id_n].efg

        vals_A_iso.append(hyp.aiso)
        vals_T_dip.append(hyp.tdip)
        vals_kappa.append(efg.nqcc)
        vals_eta.append(efg.eta)

    A_iso = sum(vals_A_iso)/len(vals_A_iso)
    T_dip = sum(vals_T_dip)/len(vals_T_dip)
    kappa = sum(vals_kappa)/len(vals_kappa)
    eta = sum(vals_eta)/len(vals_eta)

    print "  name:", filename
    print " atoms:", ids_nitrogen
    print "g_para:", g_para
    print "A_para:", A_para
    print "g_perp:", g_perp
    print "A_perp:", A_perp
    print " A_iso:", vals_A_iso, A_iso
    print " T_dip:", vals_T_dip, T_dip
    print " kappa:", vals_kappa, kappa
    print "   eta:", vals_eta, eta
    print "\n"
