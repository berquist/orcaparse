#!/usr/bin/env python2

import piratechem as pc
import mmap
import itertools
import numpy as np
import re
from piratechem.utils import *

class Atom(pc.atom.Atom):
    """
    Allow each atom to contain more specific quantum chemical properties
    than piratechem can currently handle.
    """
    def __init__(self, name, r):
        pc.atom.Atom.__init__(self, name, r)

        # storage for hyperfine values
        self.atensor = []
        self.aiso = 0.0

    def __str__(self):
        s = "Atom(%s, [%6.3f, %6.3f, %6.3f])"
        return s % (self.name, self.posx, self.posy, self.posz)

class GTensor():
    """
    Hold all of the fields that may be present in the output file
    from a g-tensor calculation.
    """
    def __init__(self):
        pass

class ATensor():
    """
    Hold all of the fields that may be present in the output file
    from nuclear property calculations.
    """
    def __init__(self):
        pass

class ORCAParser:
    """
    A parser for ORCA input and output files.

    :param file_name: the ORCA input file name
    :type file_name: str

    >>> molecule = OrcaParser('caffeine.inp')
    >>> print molecule.charge
    >>> molecule.scf['guess'] = 'pmodel'
    >>> molecule.inputblock += 'blyp def2-svp def2-svp/c'
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self.orcafile = []
        self.molecule = []
        self._has_coords = False
        self.keywords = []
        self.blocks = []

        self.gtensor = []
        self.giso = 0.0

    def __str__(self):
        # BUGS BUGS BUGS FIX FIX FIX
        for line in self.orcafile:
            print line

    def reset(self):
        pass

    def load(self):
        """
        Load the entire file into the object,
        without doing any processing on it.
        """
        handle = open(self.file_name, "r+b")
        self.orcafile = handle.readlines()

    def load_map(self):
        """
        Load the entire file as a memory-mapped object,
        without doing any processing on it.
        """
        handle = open(self.file_name, "r+b")
        self.orcamap = mmap.mmap(handle.fileno(), 0)
        self.orcamap.flush()

    def _parse_molecule(self):
        pass

    def _parse_inputblock(self):
        pass

    def _parse_fields(self):
        pass

    def _parse_charges(self):
        pass

    def _extract(self, keyword):
        pass

    def format_inputblock(self):
        pass

    def format_molecule(self):
        pass

    def format_multipole_fields(self):
        pass

    def to_string(self):
        pass

    def save(self, output_name):
        handle = open(output_name, "w")
        handle.write(self.to_string())
        handle.close()

    def _calc_interatomic_distance(self):
        distance_matrix = np.zeros((self.natoms, self.natoms))
        for j in xrange(0, self.natoms):
            for i in xrange(j+1, self.natoms):
                distance_matrix[i][j] = np.linalg.norm(self.molecule[i].r - self.molecule[j].r)
                distance_matrix[j][i] = distance_matrix[i][j]

        self.distance_matrix = distance_matrix

    def find_element(self, element_name):
        """
        Return a list of indices corresponding to all instances of an
        element in the molecule.
        """
        result_list = []
        for idx, atom in enumerate(self.molecule):
            if (atom.name == element_name):
                result_list.append(idx)
        return result_list

    def pair_distance(self, atomidx1, atomidx2):
        """
        Return the distance between the two atoms with the given indices.
        """
        return self.distance_matrix[atomidx1][atomidx2]

    def pair_distances(self, atoms_index):
        pass

    def interatomic_distance(self, atoms_index):
        """
        returns the pairwise distance array between two atoms in file

        :param atoms_index: contains index of positions of atoms of interest, 0 is first atom
        :type atoms_index: list
        """
        pass

class ORCAInputParser(ORCAParser):
    """
    A parser for ORCA input files.
    """
    def __init__(self, file_name):
        ORCAParser.__init__(self, file_name)

class ORCAOutputParser(ORCAParser):
    """
    A parser for ORCA output files.
    """
    def __init__(self, file_name):
        ORCAParser.__init__(self, file_name)
        self.load()
        self._get_coords()
        self._calc_interatomic_distance()

    def get_input_file(self):
        """
        """
        searchstr = "INPUT FILE"
        for i, line in enumerate(self.orcafile):
            if (line.find(searchstr) > -1):
                idxstart = i + 2

        # [00] NAME = SVP_decontract.inp
        # [01] |  1> ! uks pbe0 def2-svp def2-svp/jk ri rijk pmodel somf(1x) noautostart tightscf grid5 decontract
        # |  2>
        # |  3> %pal
        # |  4> nprocs 16
        # |  5> end
        # |  6>
        # |  7> * xyzfile 0 2 SVP_decontract.xyz *
        # |  8>
        # |  9> %eprnmr
        # | 10>  tol 1e-10
        # | 11>  gtensor 1
        # | 12>  ori -3
        # | 13>  nuclei = all N  { aiso, adip, aorb, fgrad, rho }
        # | 14>  nuclei = all Cu { aiso, adip, aorb, fgrad, rho }
        # | 15>  nuclei = all H  { aiso, adip, aorb, fgrad, rho }
        # | 16>  nuclei = all O  { aiso, adip, aorb, fgrad, rho }
        # | 17>  nuclei = all C  { aiso, adip, aorb, fgrad, rho }
        # | 18>  printlevel 5
        # | 19>  end
        # | 20>
        # | 21> %method
        # | 22>  z_tol 1e-10
        # | 23>  end
        # | 24>
        # | 25>
        # | 26>                          ****END OF INPUT****

        # search for the end of input so we don't loop like crazy
        searchstr = "****END OF INPUT****"
        for i, line in enumerate(self.orcafile[idxstart:]):
            if (line.find(searchstr) > -1):
                idxend = idxstart + i

        # deck = enumerate(line.strip().split() for line in self.orcafile[idxstart:idxend])
        # templist = list(deck)
        # decklist = [line[1] for line in templist[:]]
        inpdeck = [line.strip().split() for line in self.orcafile[idxstart:idxend]]

        # for line in inpdeck:
        #     # skip over lines that don't have anything
        #     if (len(line) == 2):
        #         continue
        #     # try and match the general input first
        #     if (line[2] == '!'):
        #         for word in line[3:]:
        #             self.keywords.append(word)
        #     # match an input block
        #     if (line[2][0] == '%'):
        #         blockname = line[2][1:].strip()
        #         print blockname
        #         self.blockname = dict()
        #         continue
        #         while (line[2] != 'end'):
        #             print line
        #             self.blockname[line[2]] = line[3]
        #             continue
        #         print self.blockname
        #         self.blocks['{}'.format(blockname)] = blockname
        # print self.blocks

    def get_energy(self, string_to_search = None):
        pass

    def _set_output_keywords(self):
        pass

    def check_method_type(self):
        pass

    def _get_coords(self):
        """
        Retrieve the molecular coordinates, in angstroms.
        The block looks something like this:

        ---------------------------------
        CARTESIAN COORDINATES (ANGSTROEM)
        ---------------------------------
        C     66.175587   43.218200   48.024107
        H     65.451645   42.402606   48.062464
        H     66.821636   43.047524   47.162541
        ...
        H     68.914798   48.777753   52.569627
        H     73.159381   46.517264   53.207480

        ----------------------------
        CARTESIAN COORDINATES (A.U.)
        ...
        """
        idx = 0
        searchstr = "CARTESIAN COORDINATES (ANGSTROEM)"
        for i, line in enumerate(self.orcafile):
            if (line.find(searchstr) > -1):
                # we add 2 to start at the first atom in the coordinate block
                idx = i + 2
                self._has_coords = True
                break

        for i, line in enumerate(self.orcafile[idx:]):
            if (line == "\n"):
                break
            atom = line.split()
            self.molecule.append(Atom(atom[0], np.array([atom[1], atom[2], atom[3]], dtype=np.float64)))
            self.molecule[i].index = i
        self.natoms = len(self.molecule)

    def get_gtensor(self):
        """
        Extract the electronic g-tensor from the output file.
        """
        idx = 0
        searchstr = "ELECTRONIC G-MATRIX"
        for i, line in enumerate(self.orcafile):
            if (line.find(searchstr) > -1):
                # we add 4 to start at the first row in the g-matrix
                idx = i + 4
                break
        if (idx == 0): return self.gtensor, self.giso

        # Here is a sample of what we would like to parse:
        # -------------------
        # ELECTRONIC G-MATRIX
        # -------------------

        #  The g-matrix:
        # [00]              2.1766588   -0.0419455    0.0785780
        # [01]             -0.0456024    2.1503399    0.0062106
        # [02]              0.0803984    0.0063831    2.0695626

        #  gel          2.0023193    2.0023193    2.0023193
        # [05] gRMC         2.0012821    2.0012821    2.0012821
        # [06] gDSO(1el)    0.0005885    0.0007469    0.0008681
        # [07] gDSO(2el)   -0.0002227   -0.0002837   -0.0003322
        # [08] gDSO(tot)    0.0003658    0.0004632    0.0005359
        # [09] gPSO(1el)    0.0334924    0.2332266    0.3909197
        # [10] gPSO(2el)   -0.0134555   -0.0947635   -0.1580697
        # [11] gPSO(tot)    0.0200369    0.1384632    0.2328500
        #            ----------   ----------   ----------
        # [13] g(tot)       2.0216853    2.1402089    2.2346690 iso=  2.1321877
        # [14] Delta-g      0.0193660    0.1378897    0.2323497 iso=  0.1298685
        #  Orientation:
        # [16]  X          -0.4921117    0.2594594   -0.8309675
        # [17]  Y          -0.2090765    0.8913857    0.4021425
        # [18]  Z           0.8450521    0.3716348   -0.3844145

        self.gel = 2.0023193
        # first, we parse the orientation-dependent g-matrix
        self.gmatrix = np.array([self.orcafile[idx+0].split(),
                                 self.orcafile[idx+1].split(),
                                 self.orcafile[idx+2].split()], dtype=np.float64)

        # then, we parse the decomposition of delta_g into its individual contributions
        self.grmc = np.asanyarray(self.orcafile[idx+5].split()[1:], dtype=np.float64)
        self.gdso1el = np.asanyarray(self.orcafile[idx+6].split()[1:], dtype=np.float64)
        self.gdso2el = np.asanyarray(self.orcafile[idx+7].split()[1:], dtype=np.float64)
        if (self.orcafile[idx+8].split()[1:] == ['----------', '----------']):
            print "derp, we fail when %eprnmr printlevel != 5"
        self.gdsotot = np.asanyarray(self.orcafile[idx+8].split()[1:], dtype=np.float64)
        self.gpso1el = np.asanyarray(self.orcafile[idx+9].split()[1:], dtype=np.float64)
        self.gpso2el = np.asanyarray(self.orcafile[idx+10].split()[1:], dtype=np.float64)
        self.gpsotot = np.asanyarray(self.orcafile[idx+11].split()[1:], dtype=np.float64)

        # then, we take the final tensors plus their isotropic values
        gtottmp = self.orcafile[idx+13].split()[1:]
        delgtmp = self.orcafile[idx+14].split()[1:]
        self.giso, self.dgiso = float(gtottmp[4]), float(delgtmp[4])
        self.gtensor = np.array([gtottmp[0], gtottmp[1], gtottmp[2]], dtype=np.float64)
        self.delgtensor = np.array([delgtmp[0], delgtmp[1], delgtmp[2]], dtype=np.float64)

        # finally, we take the orientation
        self.gori = np.array([self.orcafile[idx+16].split()[1:],
                              self.orcafile[idx+17].split()[1:],
                              self.orcafile[idx+18].split()[1:]], dtype=np.float64)

        return self.gtensor, self.giso

    def _get_nuclear_atom(self, atom):
        """
        Get all of the nuclear properties for a single atom.
        """

        # if the atom doesn't have coordinates, this isn't going to work
        if not self._has_coords: return ([], 0.0)

        # first, find the atom
        idx_nucleus = 0
        searchstr = "Nucleus\s+{}{}".format(atom.index, atom.name)
        for i, line in enumerate(self.orcafile):
            if (re.search(searchstr, line) is not None):
                # the first line we match contains information
                idx_nucleus = i
                break
        if (idx_nucleus == 0): return atom.atensor, atom.aiso

        # gather the hyperfine information
        idx = 0
        searchstr = "Raw HFC matrix (all values in MHz):"
        for i, line in enumerate(self.orcafile[idx_nucleus:]):
            if (line.find(searchstr) > -1):
                # start at the first row of the A-matrix
                idx = idx_nucleus + i + 1
                break
        if (idx == 0): return atom.atensor, atom.aiso

        # Raw HFC matrix (all values in MHz):
        # [00]             1.4411      -0.0380       0.0783
        # [01]            -0.0431       1.0183      -0.0205
        # [02]             0.0840      -0.0263       1.1855

        # [04] A(FC)       1.2026       1.2026       1.2026
        # [05] A(SD)      -0.1901      -0.0700       0.2602
        # [06] A(SO)       0.0005       0.0302       0.0065 A(PC) =     0.0124
        #         ----------   ----------   ----------
        # [08] A(Tot)      1.0130       1.1627       1.4693 A(iso)=     1.2150
        # Orientation:
        # [10]  X       0.0729906   -0.2793657    0.9574065
        # [11]  Y       0.9936103   -0.0624924   -0.0939856
        # [12]  Z       0.0860870    0.9581490    0.2730193

        # first, we parse the orientation-dependent A-matrix
        atom.amatrix = np.array([self.orcafile[idx+0].split(),
                                 self.orcafile[idx+1].split(),
                                 self.orcafile[idx+2].split()], dtype=np.float64)

        # then, we parse the decomposition of A into its individual contributions
        atom.afc = np.asanyarray(self.orcafile[idx+4].split()[1:], dtype=np.float64)
        atom.asd = np.asanyarray(self.orcafile[idx+5].split()[1:], dtype=np.float64)
        asotmp = self.orcafile[idx+6].split()[1:]
        atom.aso, atom.apc = np.array([asotmp[0], asotmp[1], asotmp[2]], dtype=np.float64), float(asotmp[-1])
        atottmp = self.orcafile[idx+8].split()[1:]
        atom.atensor, atom.aiso = np.array([atottmp[0], atottmp[1], atottmp[2]], dtype=np.float64), float(atottmp[-1])

        # finally, we take the orientation
        atom.aori = np.array([self.orcafile[idx+10].split()[1:],
                              self.orcafile[idx+11].split()[1:],
                              self.orcafile[idx+12].split()[1:]], dtype=np.float64)

        return atom.atensor, atom.aiso

    def _get_nuclear(self):
        """
        Get all of the fields that may be present in the output file
        from nuclear property calculations.
        """
        for atom in self.molecule:
            self.get_nuclear_atom(atom)

    def get_hyperfine(self, atom):
        """

        """
        return self._get_nuclear_atom(atom)

    def nitrogen_hyperfine_1st(self):
        """
        Return the hyperfine tensor of the closest nitrogen to the copper.
        """
        # Searching algorithm:
        # (1) Find all the nitrogens that are present and assume that there's only one copper
        # (2) Return the nitrogen-copper distances in the order of the indices from (1)
        # (3) The list position of the match in (2) corresponds to the list position
        #     of the id from (1)...phew
        # (4) Retrieve the stuff from this new id we've found

        # this lets us fail silently on bad output files
        if not self._has_coords: return ([], 0)
        ids_n = self.find_element("N")
        try:
            id_cu = self.find_element("Cu")[0]
        except IndexError:
            return
        distances = [self.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(one_smallest(distances))]
        hyperfine = self.get_hyperfine(self.molecule[idx])
        return hyperfine[0], idx

    def nitrogen_hyperfine_2nd(self):
        """
        Return the hyperfine tensor of the second-closest nitrogen to the copper.
        """

        if not self._has_coords: return ([], 0)
        ids_n = self.find_element("N")
        try:
            id_cu = self.find_element("Cu")[0]
        except IndexError:
            return
        distances = [self.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(two_smallest(distances)[1])]
        hyperfine = self.get_hyperfine(self.molecule[idx])
        return hyperfine[0], idx

    def nitrogen_hyperfine(self):
        """
        """
        pass

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="")
    parser.add_argument(dest="orcaname", metavar="<orca filename>", nargs="+", type=str, default=None, help="")
    args = parser.parse_args()

    namelist = args.orcaname

    s = "{:>34s} {:>28s} {:>3s} {:>25s} {:>3s} {:<s}"
    print s.format("g-tensor", "nitrogen_close", "id", "nitrogen_far", "id", "name")

    for name in namelist:

        orcafile = ORCAOutputParser(name)
        gtensor, giso = orcafile.get_gtensor()
        atensor_close, id_close = orcafile.nitrogen_hyperfine_1st()
        atensor_far, id_far = orcafile.nitrogen_hyperfine_2nd()

        s = "{:>28s} {:>28s} {:>3d} {:>25s} {:>3d} {:<s}"
        print s.format(gtensor, atensor_close, id_close, atensor_far, id_far, name)

