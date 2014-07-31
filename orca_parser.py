#!/usr/bin/env python2

import numpy as np
from numpy import nan
import mmap
import re
import os

import piratechem as pc
from piratechem.utils import one_smallest, two_smallest
from piratechem.utils import only_numerics

from atom import Atom
from molecule import Molecule

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
        if self.orcafile == []:
            return "empty"
        else:
            return "loaded"

    def reset(self):
        pass

    def load(self):
        """
        Load the entire file into the object,
        without doing any processing on it.
        """
        handle = open(self.file_name, "r+b")
        self.orcafile = handle.readlines()
        handle.close()

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
        for j in range(0, self.natoms):
            for i in range(j+1, self.natoms):
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
        (case-sensitive, first match only)
        """
        if start_index < 0:
            start_index = 0
        for idx, line in enumerate(self.orcafile[start_index:]):
            if (line.find(string_to_search) > -1):
                return idx
        return -1

    def get_string_indices(self, string_to_search, start_index):
        """
        Returns a list of indices for the lines containing the given
        string. (case-sensitive, all matches)
        """
        if start_index < 0:
            start_index = 0
        idxs = []
        for idx, line in enumerate(self.orcafile[start_index:]):
            if (line.find(string_to_search) > -1):
                idxs.append(idx)
        return idxs

    def get_regex_index(self, regex_to_search, start_index=0):
        """
        Returns the index for the line matching the given regular
        expression string. (case-sensitive, first match only)
        """
        if start_index < 0:
            start_index = 0
        for idx, line in enumerate(self.orcafile[start_index:]):
            if (re.search(regex_to_search, line) is not None):
                return idx
        return -1

    def get_regex_indices(self, regex_to_search, start_index=0):
        """
        Returns a list of indices for the lines matching the given
        regular expression string. (case-sensitive, all matches)
        """
        if start_index < 0:
            start_index = 0
        idxs = []
        for idx, line in enumerate(self.orcafile[start_index:]):
            if (re.search(regex_to_search, line) is not None):
                idxs.append(idx)
        return idxs

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
            self.indices_hyperfine = []
            self.indices_nmr = []
            self._extract_molecule_nuclear()
            self._extract_molecule_euler()
            self._extract_spin_contamination()

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

        # remove lines that only have "| num>"
        trimmed_inpdeck = (line[2:] for line in raw_inpdeck if len(line) > 2)

        # the first line contains the short file name
        input_file_name = trimmed_inpdeck.next()[0]

        # for line in trimmed_inpdeck:
        #     # match general input
        #     if (line[0] == '!'):
        #         for word in line[3:]:
        #             self.keywords.append(word)
        #     # match an input block
        #     # if (line[0][0] == '%'):
        #     #     blockname = line[0][1:].strip()
        #     #     self.blockname = dict()
        #     #     raw_inpdeck.next()
        #     #     while (line[0] != 'end'):
        #     #         print line
        #     #         self.blockname[line[0]] = line[1]
        #     #         continue
        #     #     self.blocks['{}'.format(blockname)] = blockname
        # print self.keywords
        # print self.blocks

    def get_energy(self, string_to_search = "Total Energy"):
        # Total Energy       :         -814.97149364 Eh          -22176.50177 eV
        idx = self.get_string_index(string_to_search)
        energy = float(self.orcafile[idx].split()[3])
        return energy

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
        idxs = []
        searchstr = "ELECTRONIC G-MATRIX"
        for i, line in enumerate(self.orcafile):
            if (line.find(searchstr) > -1):
                # we add 4 to start at the first row in the g-matrix
                idx = i + 4
                idxs.append(idx)
        if (idxs == []): return

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

        for idx in idxs:
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

    def return_atom_hyperfine(self, atom):
        """
        Return the hyperfine tensor and isotropic hyperfine value for the given atom.
        """
        return atom.hyperfine.atensor, atom.hyperfine.aiso

    def return_atom_nqi(self, atom):
        """
        Return the nuclear quadrupolar interaction (NQI) tensor,
        the NQ coupling constant, and the asymmetry parameter eta
        for the given atom.
        """
        return atom.efg.p, atom.efg.nqcc, atom.efg.eta

    def _extract_molecule_nuclear(self):
        """
        Get all of the fields that may be present in the output file
        from nuclear property calculations.
        """
        for atom in self.molecule:
            self._extract_atom_hyperfine(atom)
            self._extract_atom_shifts(atom)

    def _extract_atom_hyperfine(self, atom):
        """
        Get all of the hyperfine properties for a single atom.
        """

        # first, find the atom
        searchstr = "Nucleus\s+{}{}".format(atom.index, atom.name)
        idxs_nucleus = self.get_regex_indices(searchstr)
        if (idxs_nucleus == []): return
        self.indices_hyperfine.append(atom.index)

        for idx_nucleus in idxs_nucleus:
        # gather the hyperfine information
            searchstr = "Raw HFC matrix (all values in MHz):"
            idx_hyp = self.get_string_index(searchstr, idx_nucleus)
            # start searching from where we found the atom
            # the arrays start at idx_nucleus + idx + 1
            if (idx_hyp == -1): return
            idx_hyp += (idx_nucleus + 1)

            ###############
            ### ORCA 3.0.0+
            if (len(self.orcafile[idx_hyp].split()) == 1):
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

                atom.hyperfine.amatrix = np.array([self.orcafile[idx_hyp+1].split(),
                                                   self.orcafile[idx_hyp+2].split(),
                                                   self.orcafile[idx_hyp+3].split()], dtype=np.float64)
                atom.hyperfine.afc = np.asanyarray(self.orcafile[idx_hyp+5].split()[1:], dtype=np.float64)
                atom.hyperfine.asd = np.asanyarray(self.orcafile[idx_hyp+6].split()[1:], dtype=np.float64)
                atom.hyperfine.aso = np.asanyarray(self.orcafile[idx_hyp+8].split()[1:4], dtype=np.float64)
                atottmp = self.orcafile[idx_hyp+11].split()[1:]
                atom.hyperfine.atensor, atom.hyperfine.aiso = np.array([atottmp[0], atottmp[1], atottmp[2]], dtype=np.float64), float(atottmp[-1])
                atom.hyperfine.aori = np.array([self.orcafile[idx_hyp+13].split()[1:],
                                                self.orcafile[idx_hyp+14].split()[1:],
                                                self.orcafile[idx_hyp+15].split()[1:]], dtype=np.float64)

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

                atom.hyperfine.amatrix = np.array([self.orcafile[idx_hyp+0].split(),
                                                   self.orcafile[idx_hyp+1].split(),
                                                   self.orcafile[idx_hyp+2].split()], dtype=np.float64)

                atom.hyperfine.afc = np.asanyarray(self.orcafile[idx_hyp+4].split()[1:], dtype=np.float64)
                atom.hyperfine.asd = np.asanyarray(self.orcafile[idx_hyp+5].split()[1:], dtype=np.float64)
                asotmp = self.orcafile[idx_hyp+6].split()[1:]
                atom.hyperfine.aso, atom.hyperfine.apc = np.array([asotmp[0], asotmp[1], asotmp[2]], dtype=np.float64), float(asotmp[-1])
                atottmp = self.orcafile[idx_hyp+8].split()[1:]
                atom.hyperfine.atensor, atom.hyperfine.aiso = np.array([atottmp[0], atottmp[1], atottmp[2]], dtype=np.float64), float(atottmp[-1])

                atom.hyperfine.aori = np.array([self.orcafile[idx_hyp+10].split()[1:],
                                                self.orcafile[idx_hyp+11].split()[1:],
                                                self.orcafile[idx_hyp+12].split()[1:]], dtype=np.float64)

            atom.hyperfine._calc_eff_spin_params()

            searchstr = "Raw EFG matrix (all values in a.u.**-3):"
            # the EFG output comes after the hyperfine output
            idx_efg = self.get_string_index(searchstr, idx_hyp)
            # the arrays start at idx + 1
            if (idx_efg == -1): return
            idx_efg += idx_hyp

            # Raw EFG matrix (all values in a.u.**-3):
            # [01]           -0.0099       0.2487       0.1679
            # [02]            0.2487       0.0368      -0.1736
            # [03]            0.1679      -0.1736      -0.0269

            # [05] V(El)      -0.4732      -0.1364       0.6097
            # [06] V(Nuc)      0.6032       0.4010      -1.0042
            #         ----------   ----------   ----------
            # [08] V(Tot)      0.1300       0.2646      -0.3945
            # Orientation:
            # [10] X       0.4690732    0.6390650   -0.6095623
            # [11] Y      -0.2975070    0.7642073    0.5722559
            # [12] Z       0.8315407   -0.0870808    0.5485954

            # Quadrupole tensor eigenvalues (in MHz;Q= 0.0193 I=  1.0)
            # [16] e**2qQ            =    -1.792 MHz
            # [17] e**2qQ/(4I*(2I-1))=    -0.448 MHz
            # [18] eta               =     0.341
            #  NOTE: the diagonal representation of the SH term I*Q*I = e**2qQ/(4I(2I-1))*[-(1-eta),-(1+eta),2]

            atom.efg.vmatrix = np.array([self.orcafile[idx_efg+1].split(),
                                         self.orcafile[idx_efg+2].split(),
                                         self.orcafile[idx_efg+3].split()], dtype=np.float64)

            atom.efg.vel = np.asanyarray(self.orcafile[idx_efg+5].split()[1:], dtype=np.float64)
            atom.efg.vnuc = np.asanyarray(self.orcafile[idx_efg+6].split()[1:], dtype=np.float64)
            atom.efg.vtot = np.asanyarray(self.orcafile[idx_efg+8].split()[1:], dtype=np.float64)

            atom.efg.vori = np.array([self.orcafile[idx_efg+10].split()[1:],
                                      self.orcafile[idx_efg+11].split()[1:],
                                      self.orcafile[idx_efg+12].split()[1:]], dtype=np.float64)

            atom.efg.nqcc = float(self.orcafile[idx_efg+16].split()[-2])
            atom.efg.k = float(self.orcafile[idx_efg+17].split()[-2])
            atom.efg.eta = float(self.orcafile[idx_efg+18].split()[-1])

            atom.efg._calc_nqi_tensor()

    def _extract_atom_shifts(self, atom):
        """
        Extract all the NMR chemical shift information for a single atom.
        """

        # first, find the section
        searchstr = "CHEMICAL SHIFTS"
        idx_section = self.get_string_index(searchstr)
        if (idx_section == -1): return
        self.indices_nmr.append(atom.index)

        # now, find the atom
        searchstr = "Nucleus\s+{}{}".format(atom.index, atom.name)
        idx = self.get_regex_index(searchstr, idx_section)
        if (idx == -1): return
        idx += idx_section

        # --------------
        # [00] Nucleus   2O :
        # --------------
        # Raw-matrix : 
        # [03]        0.0001559    0.0001446    0.0000972
        # [04]        0.0000626    0.0003228   -0.0000362
        # [05]        0.0000767   -0.0000402    0.0003296

        # Diagonalized sT*s matrix:
        # [08] sDSO    0.0005205    0.0005961    0.0005458  iso=   0.0005541
        # [09] sPSO   -0.0004465   -0.0002325   -0.0001751  iso=  -0.0002847
        #        ---------------  ---------------  ---------------
        # [11] Total   0.0000739    0.0003636    0.0003707  iso=   0.0002694

        # Orientation:
        # [14]  X     0.8882944    0.1926421    0.4169197
        # [15]  Y    -0.3206033   -0.3899065    0.8632418
        # [16]  Z    -0.3288564    0.9004787    0.2845901

        # --------------
        # Nucleus   3S :
        # --------------
        # Raw-matrix : 
        #         0.0002778   -0.0000741   -0.0000508
        #        -0.0000840    0.0005165    0.0000449
        #        -0.0000388    0.0000294    0.0006807

        # Diagonalized sT*s matrix:
        # sDSO    0.0011586    0.0011695    0.0011147  iso=   0.0011476
        # sPSO   -0.0009070   -0.0006440   -0.0004167  iso=  -0.0006559
        #        ---------------  ---------------  ---------------
        # Total   0.0002516    0.0005255    0.0006979  iso=   0.0004917

        # Orientation:
        #   X     0.9570388    0.2501512   -0.1466324
        #   Y     0.2808378   -0.9255185    0.2540581
        #   Z     0.0721581    0.2843234    0.9560091

        atom.nmr.shiftmat = np.array([self.orcafile[idx+3].split(),
                                      self.orcafile[idx+4].split(),
                                      self.orcafile[idx+5].split()], dtype=np.float64)
        sdso_tmp = self.orcafile[idx+8].split()[1:]
        spso_tmp = self.orcafile[idx+9].split()[1:]
        total_tmp = self.orcafile[idx+11].split()[1:]
        atom.nmr.sdso, atom.nmr.sdso_iso = np.asanyarray(sdso_tmp[0:3], dtype=np.float64), float(sdso_tmp[-1])
        atom.nmr.spso, atom.nmr.spso_iso = np.asanyarray(spso_tmp[0:3], dtype=np.float64), float(spso_tmp[-1])
        atom.nmr.shiftpri, atom.nmr.shiftiso = np.asanyarray(total_tmp[0:3], dtype=np.float64), float(total_tmp[-1])
        atom.nmr.shiftori = np.array([self.orcafile[idx+14].split()[1:],
                                      self.orcafile[idx+15].split()[1:],
                                      self.orcafile[idx+16].split()[1:]], dtype=np.float64)
        atom.nmr._scale()
        atom.nmr._diag()

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

    def _extract_spin_contamination(self):
        """
        For spin-unrestricted jobs, find the ideal and actual values of <S**2>.
        """
        searchstr = "Expectation value of <S**2>"
        idx = self.get_string_index(searchstr)
        if (idx == -1): return

        self.molecule.ssq_actual = float(self.orcafile[idx].split()[-1])
        self.molecule.ssq_ideal = float(self.orcafile[idx+1].split()[-1])
        self.molecule.ssq_deviation = float(self.orcafile[idx+2].split()[-1])

        return

    def extract_multiple_jobs(self):
        """

        """
        pass

    def extract_timings_modules(self):
        """
        """
        searchstr = "Timings for individual modules:"
        idx = self.get_string_index(searchstr)
        if (idx == -1): return
        idx += 2

        for line in self.orcafile[idx:]:
            print line

        return
