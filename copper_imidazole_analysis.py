#!/usr/bin/env python2

import orca_parser

import numpy as np
from numpy import nan
from piratechem.utils import one_smallest, two_smallest

def g_to_mhz(gauss):
    """1 Gauss = 2.8025 MHz"""
    return gauss * 2.8025

def copper_id(orcafile):
    """Return the id of the copper atom. Assume there's only one."""
    return orcafile.find_element("Cu")[0]

def nitrogen_close_id(orcafile):
    """Return the id of the nitrogen closest to the copper."""
    if not orcafile._has_coords: return nan
    ids_n = orcafile.find_element("N")
    if (ids_n == []): return nan
    id_cu = copper_id(orcafile)
    distances = [orcafile.pair_distance(id_cu, n) for n in ids_n]
    idx = ids_n[distances.index(one_smallest(distances))]
    return idx

def nitrogen_far_id(orcafile):
    """Return the id of the nitrogen 2nd closest to the copper."""
    if not orcafile._has_coords: return nan
    ids_n = orcafile.find_element("N")
    if (ids_n == []): return nan
    id_cu = copper_id(orcafile)
    distances = [orcafile.pair_distance(id_cu, n) for n in ids_n]
    idx = ids_n[distances.index(two_smallest(distances)[1])]
    return idx

def hyperfine(orcafile, idx_nucleus):
    """Return the full hyperfine tensor for the given atom id."""
    if np.isnan(idx_nucleus): return np.array([nan, nan, nan])
    return orcafile.get_hyperfine(orcafile.molecule[idx_nucleus])[0]

def hyperfine_zz(orcafile, idx_nucleus):
    """Return the A_{zz} element for the given atom id."""
    if np.isnan(idx_nucleus): return nan
    orcafile.get_hyperfine(orcafile.molecule[idx_nucleus])
    return orcafile.molecule[idx_nucleus].amatrix[2,2]

def formatted_output_copper_imidazole(namelist):
    """
    Given pathnames of output files, pretty print the electronic g-tensor
    and the copper/far nitrogen hyperfine tensors to STDOUT.
    """
    s = "{:>34s} {:>28s} {:>3s} {:>25s} {:>3s} {:<s}"
    print s.format("g-tensor", "copper", "id", "nitrogen_far", "id", "name")

    for name in namelist:

        orcafile = orca_parser.ORCAOutputParser(name)
        gtensor, giso = orcafile.get_gtensor()
        id_cu, id_far = copper_id(orcafile), nitrogen_far_id(orcafile)
        atensor_cu, atensor_far = hyperfine(orcafile, id_cu), hyperfine(orcafile, id_far)

        a = abs(atensor_cu[-1])
        s = "{:>34s} {:>28s} {:>3d} {:>25s} {:>3d} {:<s}"
        print s.format(gtensor, atensor_cu, id_cu, atensor_far, id_far, name)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Given pathnames of ORCA output files, pretty-print the electronic g-tensor and the copper/far nitrogen hyperfine tensors.")
    parser.add_argument(dest="namelist", metavar="<orca filename>", nargs="+", type=str, default=None, help="ORCA input or output files.")
    args = parser.parse_args()
    namelist = args.namelist


    formatted_output_copper_imidazole(namelist)
