#!/usr/bin/env python2

import numpy as np
from numpy import nan
import mmap
import re
import os

import piratechem as pc
from piratechem.utils import one_smallest, two_smallest
from piratechem.utils import only_numerics

from Atom import Atom
from Molecule import Molecule

class ORCAParser:
    """
    A parser for ORCA input and output files.

    :param file_name: the ORCA input file name
    :type file_name: str

    >>> molecule = ORCAParser('caffeine.inp')
    >>> print molecule.charge
    >>> molecule.scf['guess'] = 'pmodel'
    >>> molecule.inputblock += 'blyp def2-svp def2-svp/c'
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self.orcafile = []
        self.molecule = Molecule(os.path.basename(os.path.splitext(self.file_name)[0]))
        self._has_coords = False
        self.keywords = []
        self.blocks = []

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

    def get_string_index(self, string_to_search, start_index=0):
        """
        Returns the index for the line containing the given string.
        (case-sensitive)
        """
        for idx, line in enumerate(self.orcafile[start_index:]):
            if (line.find(string_to_search) > -1):
                return idx
        return -1

    def get_regex_index(self, regex_to_search, start_index=0):
        """
        Returns the index for the line matching the given regular
        expression string. (case-sensitive)
        """
        for idx, line in enumerate(self.orcafile[start_index:]):
            if (re.search(regex_to_search, line) is not None):
                return idx
        return -1

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
        self._extract_input_file()
        self._extract_coords()
        if self._has_coords == True:
            self._calc_interatomic_distance()
            self._extract_gtensor()
            self._extract_molecule_euler()

    def _extract_input_file(self):
        """
        Extract the original input file from the output file.
        """
        searchstr = "INPUT FILE"
        idxstart = self.get_string_index(searchstr)
        idxstart += 2
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
        idxend = self.get_string_index(searchstr, idxstart)
        idxend += idxstart

        raw_inpdeck = (line.strip().split() for line in self.orcafile[idxstart:idxend])

        # the first line contains the short file name
        input_file_name = raw_inpdeck.next()[-1]
        print input_file_name

        for line in raw_inpdeck:
            # skip lines that don't contain anything
            if (len(line) == 2):
                continue
            print line
            # match general input
            if (line[2] == '!'):
                for word in line[3:]:
                    self.keywords.append(word)
            # # match an input block
            # if (line[2][0] == '%'):
            #     blockname = line[2][1:].strip()
            #     self.blockname = dict()
            #     raw_inpdeck.next()
            #     while (line[2] != 'end'):
            #         print line
            #         self.blockname[line[2]] = line[3]
            #         continue
            #     self.blocks['{}'.format(blockname)] = blockname
        print self.keywords
        print self.blocks

    def get_energy(self, string_to_search = None):
        pass

    def _set_output_keywords(self):
        pass

    def check_method_type(self):
        pass

    def _extract_coords(self):
        """
        Extract the molecular coordinates, in angstroms.
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
            self.molecule.append(Atom(i, atom[0], np.array([atom[1], atom[2], atom[3]], dtype=np.float64)))
        self.natoms = len(self.molecule)

    def return_gtensor(self):
        """
        Return the electronic g-tensor and isotropic g-value of the molecule
        from the output file.
        """
        return self.molecule.gtensor.gtensor, self.molecule.gtensor.giso

    def _extract_gtensor(self):
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
        if (idx == 0): return

        # Here is a sample of what we would like to parse (printlevel = 5):
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

        # printlevel = default:
        # -------------------
        # ELECTRONIC G-MATRIX
        # -------------------

        # The g-matrix:
        # [00] 2.1559400    0.0131394    0.0398422
        # [01] 0.0073318    2.2438835    0.0867736
        # [02] 0.0416022    0.0894444    2.0903767

        # gel          2.0023193    2.0023193    2.0023193
        # [05] gRMC         2.0013017    2.0013017    2.0013017
        # [06] gDSO(tot)    0.0003267    0.0004160    0.0005318
        # [07] gPSO(tot)    0.0390665    0.1584851    0.2874688
        # ----------   ----------   ----------
        # [09] g(tot)       2.0406959    2.1602051    2.2893042 iso=  2.1634017
        # [10] Delta-g      0.0383766    0.1578858    0.2869849 iso=  0.1610824
        # Orientation:
        # [12] X          -0.2796999    0.9393265    0.1985791
        # [13] Y          -0.3699085   -0.2963004    0.8805531
        # [14] Z           0.8859659    0.1728345    0.4303401

        # first, we parse the orientation-dependent g-matrix
        self.molecule.gtensor.gmatrix = np.array([self.orcafile[idx+0].split(),
                                                  self.orcafile[idx+1].split(),
                                                  self.orcafile[idx+2].split()], dtype=np.float64)

        # then, we parse the decomposition of delta_g into its individual contributions
        # we need to handle both the standard output level and %eprnmr printlevel = 5
        self.molecule.gtensor.grmc = np.asanyarray(self.orcafile[idx+5].split()[1:], dtype=np.float64)

        if (self.orcafile[idx+8].split()[1:] == ['----------', '----------']):
            self.molecule.gtensor.gdsotot = np.asanyarray(self.orcafile[idx+6].split()[1:], dtype=np.float64)
            self.molecule.gtensor.gpsotot = np.asanyarray(self.orcafile[idx+7].split()[1:], dtype=np.float64)
            gtottmp = self.orcafile[idx+9].split()[1:]
            delgtmp = self.orcafile[idx+10].split()[1:]
            self.molecule.gtensor.giso, self.molecule.gtensor.dgiso = float(gtottmp[4]), float(delgtmp[4])
            self.molecule.gtensor.gtensor = np.array([gtottmp[0], gtottmp[1], gtottmp[2]], dtype=np.float64)
            self.molecule.gtensor.delgtensor = np.array([delgtmp[0], delgtmp[1], delgtmp[2]], dtype=np.float64)
            self.molecule.gtensor.gori = np.array([self.orcafile[idx+12].split()[1:],
                                                   self.orcafile[idx+13].split()[1:],
                                                   self.orcafile[idx+14].split()[1:]], dtype=np.float64)

        else:
            self.molecule.gtensor.gdso1el = np.asanyarray(self.orcafile[idx+6].split()[1:], dtype=np.float64)
            self.molecule.gtensor.gdso2el = np.asanyarray(self.orcafile[idx+7].split()[1:], dtype=np.float64)
            self.molecule.gtensor.gdsotot = np.asanyarray(self.orcafile[idx+8].split()[1:], dtype=np.float64)
            self.molecule.gtensor.gpso1el = np.asanyarray(self.orcafile[idx+9].split()[1:], dtype=np.float64)
            self.molecule.gtensor.gpso2el = np.asanyarray(self.orcafile[idx+10].split()[1:], dtype=np.float64)
            self.molecule.gtensor.gpsotot = np.asanyarray(self.orcafile[idx+11].split()[1:], dtype=np.float64)
            # then, we take the final tensors plus their isotropic values
            gtottmp = self.orcafile[idx+13].split()[1:]
            delgtmp = self.orcafile[idx+14].split()[1:]
            self.molecule.gtensor.giso, self.molecule.gtensor.dgiso = float(gtottmp[4]), float(delgtmp[4])
            self.molecule.gtensor.gtensor = np.array([gtottmp[0], gtottmp[1], gtottmp[2]], dtype=np.float64)
            self.molecule.gtensor.delgtensor = np.array([delgtmp[0], delgtmp[1], delgtmp[2]], dtype=np.float64)
            # finally, we take the orientation
            self.molecule.gtensor.gori = np.array([self.orcafile[idx+16].split()[1:],
                                  self.orcafile[idx+17].split()[1:],
                                  self.orcafile[idx+18].split()[1:]], dtype=np.float64)

        return

    def _get_nuclear_atom(self, atom):
        """
        Get all of the nuclear properties for a single atom.
        """

        # first, find the atom
        searchstr = "Nucleus\s+{}{}".format(atom.index, atom.name)
        idx_nucleus = self.get_regex_index(searchstr)
        if (idx_nucleus == -1): return atom.hyperfine.atensor, atom.hyperfine.aiso

        # gather the hyperfine information
        searchstr = "Raw HFC matrix (all values in MHz):"
        idx = self.get_string_index(searchstr, idx_nucleus)
        # start searching from where we found the atom
        # the arrays start at idx_nucleus + idx + 1
        if (idx == -1): return atom.hyperfine.atensor, atom.hyperfine.aiso
        idx += (idx_nucleus + 1)

        ###############
        ### ORCA 3.0.0
        if (len(self.orcafile[idx].split()) == 1):
            # Raw HFC matrix (all values in MHz):
            # [00] ------------------------------
            # [01] 2.6153               0.1906               0.1416
            # [02] 0.1871               2.2165               0.2023
            # [03] 0.1370               0.1969               2.3443

            # [05] A(FC)           2.3761               2.3761               2.3761
            # [06] A(SD)          -0.3152              -0.0925               0.4078
            # [07] A(ORB+DIA)     -0.0013               0.0436               0.0057    A(PC) =   0.0160
            # [08] A(ORB)         -0.0014               0.0435               0.0058    A(PC) =   0.0160
            # [09] A(DIA)          0.0001               0.0000              -0.0000    A(PC) =   0.0000
            # ----------           ----------           ----------
            # [11] A(Tot)          2.0595               2.3271               2.7895    A(iso)=    2.3921
            # Orientation:
            # [13] X           0.1564822            0.5846922            0.7960203
            # [14] Y          -0.8421894           -0.3420471            0.4167983
            # [15] Z           0.5159752           -0.7356214            0.4388972

            atom.hyperfine.amatrix = np.array([self.orcafile[idx+1].split(),
                                               self.orcafile[idx+2].split(),
                                               self.orcafile[idx+3].split()], dtype=np.float64)
            atom.hyperfine.afc = np.asanyarray(self.orcafile[idx+5].split()[1:], dtype=np.float64)
            atom.hyperfine.asd = np.asanyarray(self.orcafile[idx+6].split()[1:], dtype=np.float64)
            atottmp = self.orcafile[idx+11].split()[1:]
            atom.hyperfine.atensor, atom.hyperfine.aiso = np.array([atottmp[0], atottmp[1], atottmp[2]], dtype=np.float64), float(atottmp[-1])
            atom.hyperfine.aori = np.array([self.orcafile[idx+13].split()[1:],
                                            self.orcafile[idx+14].split()[1:],
                                            self.orcafile[idx+15].split()[1:]], dtype=np.float64)

        ##############
        ### ORCA 2.9.1
        else:
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
            atom.hyperfine.amatrix = np.array([self.orcafile[idx+0].split(),
                                               self.orcafile[idx+1].split(),
                                               self.orcafile[idx+2].split()], dtype=np.float64)

            # then, we parse the decomposition of A into its individual contributions
            atom.hyperfine.afc = np.asanyarray(self.orcafile[idx+4].split()[1:], dtype=np.float64)
            atom.hyperfine.asd = np.asanyarray(self.orcafile[idx+5].split()[1:], dtype=np.float64)
            asotmp = self.orcafile[idx+6].split()[1:]
            atom.hyperfine.aso, atom.hyperfine.apc = np.array([asotmp[0], asotmp[1], asotmp[2]], dtype=np.float64), float(asotmp[-1])
            atottmp = self.orcafile[idx+8].split()[1:]
            atom.hyperfine.atensor, atom.hyperfine.aiso = np.array([atottmp[0], atottmp[1], atottmp[2]], dtype=np.float64), float(atottmp[-1])

            # finally, we take the orientation
            atom.hyperfine.aori = np.array([self.orcafile[idx+10].split()[1:],
                                            self.orcafile[idx+11].split()[1:],
                                            self.orcafile[idx+12].split()[1:]], dtype=np.float64)

        return atom.hyperfine.atensor, atom.hyperfine.aiso

    def _extract_molecule_euler(self):
        """
        Extract all of the Euler rotation angles from the output file.
        """
        endline = "-------------------------------------------------------------------"

        searchstr = "Euler rotation of hyperfine tensor to g-tensor"
        idx = self.get_string_index(searchstr)
        if (idx == -1): return
        idx += 7

        for line in self.orcafile[idx:]:
            tmpline = line.split()
            if (tmpline == []) or (tmpline[0] == endline): break
            idx_nucleus = int(only_numerics(tmpline[0]))
            tmpline[:] = [float(i) for i in tmpline[1:]]
            self.molecule[idx_nucleus].euler.hyperfine.alpha = tmpline[0]
            self.molecule[idx_nucleus].euler.hyperfine.beta = tmpline[1]
            self.molecule[idx_nucleus].euler.hyperfine.gamma = tmpline[2]
            self.molecule[idx_nucleus].euler.hyperfine.ax = tmpline[3]
            self.molecule[idx_nucleus].euler.hyperfine.ay = tmpline[4]
            self.molecule[idx_nucleus].euler.hyperfine.az = tmpline[5]

        searchstr = "Euler rotation of Electric Field Gradient (EFG) tensor to g-tensor"
        idx = self.get_string_index(searchstr)
        if (idx == -1): return
        idx += 7

        for line in self.orcafile[idx:]:
            tmpline = line.split()
            if (tmpline == []) or (tmpline[0] == endline): break
            idx_nucleus = int(only_numerics(tmpline[0]))
            tmpline[:] = [float(i) for i in tmpline[1:]]
            self.molecule[idx_nucleus].euler.efg.alpha = tmpline[0]
            self.molecule[idx_nucleus].euler.efg.beta = tmpline[1]
            self.molecule[idx_nucleus].euler.efg.gamma = tmpline[2]
            self.molecule[idx_nucleus].euler.efg.efgx = tmpline[3]
            self.molecule[idx_nucleus].euler.efg.efgy = tmpline[4]
            self.molecule[idx_nucleus].euler.efg.efgz = tmpline[5]

        return

    def _get_nuclear(self):
        """
        Get all of the fields that may be present in the output file
        from nuclear property calculations.
        """
        for atom in self.molecule:
            self._get_nuclear_atom(atom)

    def get_hyperfine(self, atom):
        """

        """
        return self._get_nuclear_atom(atom)

    def extract_multiple_jobs(self):
        """

        """
        pass

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
        if not self._has_coords: return (np.array([nan, nan, nan]), nan)
        ids_n = self.find_element("N")
        try:
            id_cu = self.find_element("Cu")[0]
        except IndexError:
            return
        distances = [self.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(one_smallest(distances))]
        hyperfine = self.get_hyperfine(self.molecule[idx])
        return hyperfine[0], idx, self.molecule[idx]

    def nitrogen_hyperfine_2nd(self):
        """
        Return the hyperfine tensor of the second-closest nitrogen to the copper.
        """

        if not self._has_coords: return (np.array([nan, nan, nan]), nan)
        ids_n = self.find_element("N")
        try:
            id_cu = self.find_element("Cu")[0]
        except IndexError:
            return
        distances = [self.pair_distance(id_cu, n) for n in ids_n]
        idx = ids_n[distances.index(two_smallest(distances)[1])]
        hyperfine = self.get_hyperfine(self.molecule[idx])
        return hyperfine[0], idx, self.molecule[idx]

# if __name__ == "__main__":
#     from orcaparse import *
#     import argparse

#     parser = argparse.ArgumentParser(description="A toolbox for creating and parsing ORCA input and output files from the command line or Python scripts.")
#     parser.add_argument(dest="namelist", metavar="<orca filename>", nargs="+", type=str, default=None, help="ORCA input or output files.")
#     args = parser.parse_args()
#     namelist = args.namelist
