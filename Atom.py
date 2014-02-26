#!/usr/bin/env python2

import piratechem as pc
import numpy as np
from numpy import nan
import scipy.linalg as spl

class Atom(pc.atom.Atom):
    """
    Allow each atom to contain more specific quantum chemical properties
    than piratechem can currently handle.
    """
    def __init__(self, index, name, r):
        pc.atom.Atom.__init__(self, name, r)

        self.index = index

        self.nmr = NMR()
        self.hyperfine = Hyperfine()
        self.efg = EFG()
        self.euler = Euler()

    def __str__(self):
        s = "Atom(%d, %s, [%6.3f, %6.3f, %6.3f])"
        return s % (self.index, self.name, self.posx, self.posy, self.posz)

class Euler:
    """
    Store all possible Euler angle information for a single atom.
    """
    def __init__(self):
        self.hyperfine = self.Hyperfine()
        self.efg = self.EFG()

    class Hyperfine:
        def __init__(self):
            self.alpha = self.beta = self.gamma = nan
            self.ax = self.ay = self.az = nan

        def __str__(self):
            s = "EulerHyperfine([{0}, {1}, {2}]; [{3} {4} {5}])"
            return s.format(self.alpha, self.beta, self.gamma,
                            self.ax, self.ay, self.az)

        def return_angles(self):
            """
            Return the three angles as a NumPy row vector.
            """
            return np.array([self.alpha, self.beta, self.gamma])

    class EFG:
        def __init__(self):
            self.alpha = self.beta = self.gamma = nan
            self.efgx = self.efgy = self.efgz = nan

        def __str__(self):
            s = "EulerEFG([{0}, {1}, {2}]; [{3} {4} {5}])"
            return s.format(self.alpha, self.beta, self.gamma,
                            self.efgx, self.efgy, self.efgz)

        def return_angles(self):
            """
            Return the three angles as a NumPy row vector.
            """
            return np.array([self.alpha, self.beta, self.gamma])

class NMR:
    """
    Hold all of the fields that may be present in the output file
    from an NMR shift calculation.
    """
    def __init__(self):
        self.shiftmat = np.array([[nan, nan, nan],
                                  [nan, nan, nan],
                                  [nan, nan, nan]])
        self.sdso = np.array([nan, nan, nan])
        self.spso = np.array([nan, nan, nan])
        self.shiftpri = np.array([nan, nan, nan])
        self.sdso_iso = nan
        self.spso_iso = nan
        self.shiftiso = nan
        self.shiftori = np.array([[nan, nan, nan],
                                  [nan, nan, nan],
                                  [nan, nan, nan]])

        self.eigvals = np.array([nan, nan, nan])
        self.iso = nan

    def __str__(self):
        s = "NMR([{0} {1} {2}]; {3})"
        return s.format(self.shiftpri[0],
                        self.shiftpri[1],
                        self.shiftpri[2],
                        self.shiftiso)

    def _scale(self):
        """
        Convert the absolute values given by ORCA to ppm.
        """
        abs_to_ppm = 1e6
        self.shiftmat *= abs_to_ppm
        self.sdso *= abs_to_ppm
        self.spso *= abs_to_ppm
        self.shiftpri *= abs_to_ppm
        self.sdso_iso *= abs_to_ppm
        self.spso_iso *= abs_to_ppm
        self.shiftiso *= abs_to_ppm

    def _diag(self):
        """
        Diagonalize the raw shift matrix to get the three principal shift
        values and an isotropic result.
        """
        self.eigvals = np.sqrt(spl.eigvals(np.dot(self.shiftmat.T, self.shiftmat)).real)
        self.iso = np.sum(self.eigvals) / 3.0

class Hyperfine:
    """
    Hold all of the fields that may be present in the output file
    from nuclear property calculations.
    """
    def __init__(self):
        self.aiso = nan
        self.atensor = np.array([nan, nan, nan])
        self.amatrix = np.array([[nan, nan, nan],
                                 [nan, nan, nan],
                                 [nan, nan, nan]])

        self.afc = np.array([nan, nan, nan])
        self.asd = np.array([nan, nan, nan])
        self.aso = np.array([nan, nan, nan])
        self.apc = nan
        self.aori = np.array([[nan, nan, nan],
                              [nan, nan, nan],
                              [nan, nan, nan]])

    def __str__(self):
        s = "Hyperfine([{0} {1} {2}]; {3})"
        return s.format(self.atensor[0],
                        self.atensor[1],
                        self.atensor[2],
                        self.aiso)

class EFG:
    """
    """
    def __init__(self):
        self.vmatrix = np.array([[nan, nan, nan],
                                 [nan, nan, nan],
                                 [nan, nan, nan]])
        self.vel = np.array([nan, nan, nan])
        self.vnuc = np.array([nan, nan, nan])
        self.vtot = np.array([nan, nan, nan])
        self.vori = np.array([[nan, nan, nan],
                              [nan, nan, nan],
                              [nan, nan, nan]])

        self.nqcc = nan
        self.k = nan
        self.eta = nan

        self.px = nan
        self.py = nan
        self.pz = nan
        self.p = np.array([nan, nan, nan])

    def __str__(self):
        s = "EFG([{0} {1} {2}]; {3})"
        return s.format(self.vtot[0],
                        self.vtot[1],
                        self.vtot[2],
                        self.nqcc)

    def _calc_nqi_tensor(self):
        """
        """
        self.px = self.k * (-(1-self.eta))
        self.py = self.k * (-(1+self.eta))
        self.pz = self.k * 2
        self.p = np.array([self.px, self.py, self.pz])
