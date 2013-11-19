#!/usr/bin/env python2

import orca_parser
import numpy as np
from numpy import nan
from piratechem.utils import one_smallest, two_smallest

class CopperImidazoleAnalysis:
    def g_to_mhz(self, gauss):
        """1 Gauss = 2.8025 MHz"""
        return gauss * 2.8025

    def copper_id(self, orcafile):
        """Return the id of the copper atom. Assume there's only one."""
        return orcafile.find_element("Cu")[0]

    def nitrogen_close_id(self, orcafile):
        """Return the id of the nitrogen closest to the copper."""
        if not orcafile._has_coords: return nan
        ids_n = orcafile.find_element("N")
        if (ids_n == []): return nan
        id_cu = self.copper_id(orcafile)
        distances = [orcafile.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(one_smallest(distances))]
        return idx

    def nitrogen_far_id(self, orcafile):
        """Return the id of the nitrogen 2nd closest to the copper."""
        if not orcafile._has_coords: return nan
        ids_n = orcafile.find_element("N")
        if (ids_n == []): return nan
        id_cu = self.copper_id(orcafile)
        distances = [orcafile.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(two_smallest(distances)[1])]
        return idx

    def hyperfine(self, orcafile, idx_nucleus):
        """Return the full hyperfine tensor for the given atom id."""
        if np.isnan(idx_nucleus): return np.array([nan, nan, nan])
        return orcafile.return_atom_hyperfine(orcafile.molecule[idx_nucleus])[0]

    def hyperfine_zz(self, orcafile, idx_nucleus):
        """Return the A_{zz} element for the given atom id."""
        if np.isnan(idx_nucleus): return nan
        orcafile.get_hyperfine(orcafile.molecule[idx_nucleus])
        return orcafile.molecule[idx_nucleus].hyperfine.amatrix[2,2]
    
    def euler(self, orcafile, idx_nucleus):
        """Return the alpha, beta, gamma angles between the given hyperfine
        tensor and the g-tensor."""
        if np.isnan(idx_nucleus): return np.array([nan, nan, nan])
        return orcafile.molecule[idx_nucleus].euler.hyperfine.return_angles()

    def nqi(self, orcafile, idx_nucleus):
        """Return the nuclear quadrupolar interactions (NQI tensor, NQCC, eta)
        for the given atom id."""
        if np.isnan(idx_nucleus): return nan
        return orcafile.return_atom_nqi(orcafile.molecule[idx_nucleus])
