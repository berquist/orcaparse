#!/usr/bin/env python2

import piratechem as pc
import numpy as np
from numpy import nan

class Atom(pc.atom.Atom):
    """
    Allow each atom to contain more specific quantum chemical properties
    than piratechem can currently handle.
    """
    def __init__(self, index, name, r):
        pc.atom.Atom.__init__(self, name, r)

        # storage for hyperfine values
        self.amatrix = np.array([[nan, nan, nan],
                                 [nan, nan, nan],
                                 [nan, nan, nan]])
        self.atensor = np.array([nan, nan, nan])
        self.aiso = nan
        self.index = index
        self.hyperfine = Hyperfine()
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

class Hyperfine:
    """
    Hold all of the fields that may be present in the output file
    from nuclear property calculations.
    """
    aiso = nan
    atensor = np.array([nan, nan, nan])
    amatrix = np.array([[nan, nan, nan],
                        [nan, nan, nan],
                        [nan, nan, nan]])

    def __str__(self):
        s = "Hyperfine([{0} {1} {2}]; {3})"
        return s.format(self.atensor[0],
                        self.atensor[1],
                        self.atensor[2],
                        self.aiso)
