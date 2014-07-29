#!/usr/bin/env python2

import piratechem as pc
import numpy as np
from numpy import nan

class Molecule(pc.molecule.Molecule):
    """
    Allow each molecule to contain more specific quantum chemical properties
    than piratechem can currently handle.
    """
    def __init__(self, name):
        pc.molecule.Molecule.__init__(self, name)
        self.gtensor = GTensor()

        self.ssq_actual = nan
        self.ssq_ideal = nan
        self.ssq_deviation = nan

    def __getitem__(self, key):
        """
        """
        return self.atoms[key]

    def __str__(self):
        """
        Allow for pretty-printing a molecule.
        """
        s = ""
        for atom in self.atoms:
            s += atom.__str__() + "\n"
        return s

    def append(self, atom):
        """
        Append an Atom directly on to a Molecule.
        """
        self.atoms.append(atom)

class GTensor:
    """
    Hold all of the fields that may be present in the output file
    from an electronic g-tensor calculation.
    """
    gel = 2.0023193
    giso = nan
    dgiso = nan
    gtensor = np.array([nan, nan, nan])
    delgtensor = np.array([nan, nan, nan])
    gmatrix = np.array([[nan, nan, nan],
                        [nan, nan, nan],
                        [nan, nan, nan]])

    grmc = np.array([nan, nan, nan])
    gdso1el = np.array([nan, nan, nan])
    gdso2el = np.array([nan, nan, nan])
    gdsotot = np.array([nan, nan, nan])
    gpso1el = np.array([nan, nan, nan])
    gpso2el = np.array([nan, nan, nan])
    gpsotot = np.array([nan, nan, nan])
    gori = np.array([[nan, nan, nan],
                     [nan, nan, nan],
                     [nan, nan, nan]])

    def __str__(self):
        s = "GTensor([{0}, {1}, {2}]; {3})"
        return s.format(self.gtensor[0],
                        self.gtensor[1],
                        self.gtensor[2],
                        self.giso)

    def _diag(self):
        """
        Diagonalize the raw g-matrix to get the three principal g values
        values and an isotropic result.
        """
        self.eigvals = np.sqrt(spl.eigvals(np.dot(self.gmatrix.T, self.gmatrix)).real)
        self.iso = np.sum(self.eigvals) / 3.0

