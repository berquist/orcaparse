#!/usr/bin/env python2



import numpy as np
from numpy import nan

from . import orca_parser
from piratechem.utils import one_smallest, two_smallest

class CopperImidazoleAnalysis:
    def g_to_mhz(self, gauss):
        """1 Gauss = 2.8025 MHz"""
        return gauss * 2.8025

    def copper_id(self, orcafile):
        """Return the id of the copper atom. Assume there's only one."""
        if not orcafile._has_coords: return nan
        return orcafile.find_element("Cu")[0]

    def nitrogen_ids(self, orcafile):
        """Return a list of ids of all nitrogen atoms."""
        if not orcafile._has_coords: return nan
        return orcafile.find_element("N")

    def nitrogen_close_id(self, orcafile):
        """Return the id of the nitrogen closest to the copper."""
        ids_n = self.nitrogen_ids(orcafile)
        if (ids_n == []): return nan
        id_cu = self.copper_id(orcafile)
        distances = [orcafile.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(one_smallest(distances))]
        return idx

    def nitrogen_far_id(self, orcafile):
        """Return the id of the nitrogen 2nd closest to the copper."""
        ids_n = self.nitrogen_ids(orcafile)
        if (ids_n == []): return nan
        id_cu = self.copper_id(orcafile)
        distances = [orcafile.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(two_smallest(distances)[1])]
        return idx

    def hyperfine(self, orcafile, idx_nucleus):
        """
        Return the diagonalized hyperfine tensor for the given atom id.
        """
        if np.isnan(idx_nucleus):
            return np.array([nan, nan, nan])
        return orcafile.return_atom_hyperfine(orcafile.molecule[idx_nucleus])[0]

    def euler(self, orcafile, idx_nucleus):
        """
        Return the alpha, beta, and gamma angles between the given hyperfine
        tensor and the g-tensor.
        """
        if np.isnan(idx_nucleus):
            return np.array([nan, nan, nan])
        return orcafile.molecule[idx_nucleus].euler.hyperfine.return_angles()

    def nqi(self, orcafile, idx_nucleus):
        """
        Return the nuclear quadrupolar interactions (NQI tensor, NQCC, eta)
        for the given atom id.
        """
        if np.isnan(idx_nucleus):
            return nan
        return orcafile.return_atom_nqi(orcafile.molecule[idx_nucleus])

    def determine_copper_covalency(self, orcafile):
        """
        Calculate the covalence metric \alpha^{2} from the the equation of
        Kivelson and Nieman:

        \alpha^{2} = -\frac{A_{\parallel}}{P} + (g_{\parallel} - g_{e}) +
          \frac{3}{7}(g_{\perp} - g_{e}) + 0.04

        where P is the interaction dipolar term of the free Cu(II) ion
        (= 0.036 cm**-1 == 1079.25 MHz).

        \alpha^{2} = 1.0 -> completely ionic
        \alpha^{2} = 0.5 -> completely covalent

        Reference:
        Daniel Kivelson and Robert Neiman, Journal of Chemical Physics,
        ESR Studies on the Bonding in Copper Complexes, 1961, 35, 149-155.
        DOI: 10.1063/1.1731880
        """

        A_para = abs(self.hyperfine(orcafile, self.copper_id(orcafile))[2])
        g_principal = orcafile.molecule.gtensor.gtensor
        g_para = g_principal[2]
        g_perp = (g_principal[0] + g_principal[1]) / 2
        g_e = orcafile.molecule.gtensor.gel

        c = 299792458
        P = 0.036 * 100 * c * 10**-6

        alpha_sq = -(A_para/P) + (g_para - g_e) + (3/7)*(g_perp - g_e) + 0.04

        return alpha_sq
