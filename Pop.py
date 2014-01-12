#!/usr/bin/env python2

# from piratechem.utils import get_string_index
# from piratechem.utils import get_regex_index

import numpy as np

class Pop:
    """
    """
    type_string = ""
    def __init__(self, orcafile, type_string):
        self.idx_section = get_string_index(header)
        return

    def _find_section_start(header):
        """
        """
        idx_section = orcafile.get_string_index(header)

    def _get_chg_spin_atomic(self):
        section_string = " ATOMIC CHARGES AND SPIN POPULATIONS"
        header = type_string + section_string
        idx_increment = orcafile.get_string_index(header, idx_section)
        idx_start = idx_section + idx_increment + 2
        
        return

    def _get_chg_spin_orbital(self):
        section_string = " ORBITAL CHARGES AND SPIN POPULATIONS"
        header = type_string + section_string
        return

    def _get_chg_spin_orbital_reduced(self):
        section_string = " REDUCED ORBITAL CHARGES AND SPIN POPULATIONS"
        header = type_string + section_string
        return

    def _get_bond_orders(self):
        section_string = " BOND ORDERS (THRESH 0.05)"
        header = type_string + section_string
        return

    def _get_pop_mo_atom(self):
        section_string = " ATOM POPULATIONS PER MO"
        header = type_string + section_string
        return

    def _get_pop_mo_orbital(self):
        section_string = " ORBITAL POPULATIONS PER MO"
        header = type_string + section_string
        return

    def _get_pop_mo_orbital_reduced(self):
        section_string = " REDUCED ORBITAL POPULATIONS PER MO"
        header = type_string + section_string
        return

class Mulliken(Pop):
    """
    """
    type_string = "MULLIKEN"
    def __init__(self, orcafile, type_string):
        Pop.__init__(self, orcafile, type_string)
        return

class Loewdin(Pop):
    """
    """
    type_string = "LOEWDIN"
    def __init__(self, orcafile, type_string):
        Pop.__init__(self, orcafile, type_string)
        return

class Mayer(Pop):
    """
    """
    type_string = "MAYER"
    def __init__(self, orcafile, type_string):
        Pop.__init__(self, orcafile, type_string)
        return

class Hirshfeld(Pop):
    """
    """
    type_string = "HIRSHFELD"
    def __init__(self, orcafile, type_string):
        Pop.__init__(self, orcafile, type_string)
        return

class Density:
    """
    Parse and store density matrix-related information all in one place.
    """
    def __init__(self, orcafile):
        self._get_basis_info(orcafile)
        self._parse_density(orcafile)
        return

    # initialize with fake values until we have parsed the file
    primitive_gaussian_shells = -1
    primitive_gaussian_functions = -1
    contracted_shells = -1
    contracted_basis_functions = -1

    density = np.nan
    density_alpha = np.nan
    density_beta = np.nan
    density_spin = np.nan

    def _get_basis_info(self, orcafile):
        """
        Parse and gather some basic information required for later parsing of
        basis-dependent properties such as matrices.
        """
        searchstr = "BASIS SET STATISTICS AND STARTUP INFO"
        idx = orcafile.get_string_index(searchstr)
        idx += 2
        self.primitive_gaussian_shells = int(orcafile.orcafile[idx+0].split()[-1])
        self.primitive_gaussian_functions = int(orcafile.orcafile[idx+1].split()[-1])
        self.contracted_shells = int(orcafile.orcafile[idx+2].split()[-1])
        self.contracted_basis_functions = int(orcafile.orcafile[idx+3].split()[-1])
        return

    def _parse_density(self, orcafile):
        """
        From an ORCA output file with "%output print[p_density] 1 end", gather
        and store the alpha and beta spin densities.
        """
        searchstr = "DENSITY\n"
        idxstart = orcafile.get_string_index(searchstr) + 3
        # the dimension of a density matrix is the 
        # number of contracted basis functions
        cbf = self.contracted_basis_functions
        self.density_alpha = -np.ones((cbf, cbf))
        self.density_beta = -np.ones((cbf, cbf))
        # ORCA prints 6 columns at a time, starting from 0
        # alpha density        
        file_counter = idxstart
        row_counter = 0
        col_counter = 0
        width = 6
        while col_counter < cbf:
            while row_counter < cbf:
                self.density_alpha[row_counter][col_counter:col_counter+width] = [float(x) for x in orcafile.orcafile[file_counter].split()[1:]]
                s = "alpha r: {0:3d} c: {1:3d} f: {2:3d}"
                print(s.format(row_counter, col_counter, file_counter-idxstart))
                row_counter += 1
                file_counter += 1
            file_counter += 1
            row_counter = 0
            col_counter += 6
        # beta density
        file_counter += 1
        row_counter = 0
        col_counter = 0
        while col_counter < cbf:
            while row_counter < cbf:
                self.density_beta[row_counter][col_counter:col_counter+width] = [float(x) for x in orcafile.orcafile[file_counter].split()[1:]]
                s = "beta r: {0:3d} c: {1:3d} f: {2:3d}"
                print(s.format(row_counter, col_counter, file_counter-idxstart))
                row_counter += 1
                file_counter += 1
            file_counter += 1
            row_counter = 0
            col_counter += 6
        return

    def _calc_derived_densities(self):
        """
        Assuming that the alpha and beta spin densities have been gathered
        and stored, calculate the total density and the spin density.
        """
        self.density = self.density_alpha + self.density_beta
        self.density_spin = self.density_alpha - self.density_beta
        return
