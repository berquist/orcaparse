#!/usr/bin/env python2

import numpy as np

class Pop(object):
    """
    Hold all of the information associated with different types of population
    analyses.
    """
    def __init__(self, job, type_string = ""):
        self.job = job
        self.type_string = type_string
        self.charges = []

        self.idx_section = self.job.get_string_index(type_string +
                                                     " POPULATION ANALYSIS")
        return

    def _get_chg_atomic(self):
        section_string = " ATOMIC CHARGES"
        header = self.type_string + section_string
        idx_increment = self.job.get_string_index(header, self.idx_section)
        idx_start = self.idx_section + idx_increment + 2
        for line in self.job.orcafile[idx_start : idx_start + self.job.natoms]:
            self.charges.append(float(line.split()[-1]))
        return

    def _get_chg_spin_atomic(self):
        section_string = " ATOMIC CHARGES AND SPIN POPULATIONS"
        header = self.type_string + section_string
        idx_increment = self.job.get_string_index(header, self.idx_section)
        idx_start = self.idx_section + idx_increment + 2
        return

    def _get_chg_spin_orbital(self):
        section_string = " ORBITAL CHARGES AND SPIN POPULATIONS"
        header = self.type_string + section_string
        return

    def _get_chg_orbital_reduced(self):
        section_string = " REDUCED ORBITAL CHARGES"
        header = self.type_string + section_string
        return

    def _get_chg_spin_orbital_reduced(self):
        section_string = " REDUCED ORBITAL CHARGES AND SPIN POPULATIONS"
        header = self.type_string + section_string
        return

    def _get_bond_orders(self):
        section_string = " BOND ORDERS (THRESH 0.05)"
        header = self.type_string + section_string
        return

    def _get_pop_mo_atom(self):
        section_string = " ATOM POPULATIONS PER MO"
        header = self.type_string + section_string
        return

    def _get_pop_mo_orbital(self):
        section_string = " ORBITAL POPULATIONS PER MO"
        header = self.type_string + section_string
        return

    def _get_pop_mo_orbital_reduced(self):
        section_string = " REDUCED ORBITAL POPULATIONS PER MO"
        header = self.type_string + section_string
        return

class Mulliken(Pop):
    """
    Hold results from a Mulliken population analysis.
    """
    def __init__(self, job, type_string = "MULLIKEN"):
        Pop.__init__(self, job, type_string)
        if self.idx_section > -1:
            self._get_chg_atomic()
        return

class Loewdin(Pop):
    """
    Hold results from a Loewdin population analysis.
    """
    def __init__(self, job, type_string = "LOEWDIN"):
        Pop.__init__(self, job, type_string)
        if self.idx_section > -1:
            self._get_chg_atomic()
        return

class Mayer(Pop):
    """
    Hold results from a Mayer population analysis.
    """
    def __init__(self, job, type_string = "MAYER"):
        Pop.__init__(self, job, type_string)
        if self.idx_section > -1:

            self.NA = []
            self.QA = []
            self.VA = []
            self.BVA = []
            self.FA = []

        return

class Hirshfeld(Pop):
    """
    Hold results from a Hirshfeld population analysis.
    """
    def __init__(self, job, type_string = "HIRSHFELD"):
        Pop.__init__(self, job, type_string)
        if self.idx_section > -1:
            pass
        return

class ChElPG(Pop):
    """
    Hold results from a ChElPG (atomic) charge analysis.
    """
    def __init__(self, job, type_string = "CHELPG"):
        Pop.__init__(self, job, type_string)
        self.idx_section = self.job.get_string_index("ORCA CHELPG CHARGES GENERATION")
        if self.idx_section > -1:
            self._get_chg_atomic()
        return

    def _get_chg_atomic(self):
        section_string = " Charges"
        header = self.type_string + section_string
        idx_increment = self.job.get_string_index(header, self.idx_section)
        idx_start = self.idx_section + idx_increment + 2
        for line in self.job.orcafile[idx_start : idx_start + self.job.natoms]:
            self.charges.append(float(line.split()[-1]))
        return

class Density:
    """
    Parse and hold density matrix-related information all in one place.
    """
    def __init__(self, job):
        self._parse_basis_info(job)
        self._parse_density(job)
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

    def _parse_basis_info(self, job):
        """
        Parse and gather some basic information required for later parsing of
        basis-dependent properties such as matrices.
        """
        searchstr = "BASIS SET STATISTICS AND STARTUP INFO"
        idx = job.get_string_index(searchstr)
        idx += 2
        self.primitive_gaussian_shells = int(job.orcafile[idx+0].split()[-1])
        self.primitive_gaussian_functions = int(job.orcafile[idx+1].split()[-1])
        self.contracted_shells = int(job.orcafile[idx+2].split()[-1])
        self.contracted_basis_functions = int(job.orcafile[idx+3].split()[-1])
        return

    def _parse_density(self, job):
        """
        From an ORCA output file with "%output print[p_density] 1 end", gather
        and store the alpha and beta spin densities.
        """
        searchstr = "DENSITY\n"
        idxstart = job.get_string_index(searchstr) + 3
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
                self.density_alpha[row_counter][col_counter:col_counter+width] = [float(x) for x in job.orcafile[file_counter].split()[1:]]
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
                self.density_beta[row_counter][col_counter:col_counter+width] = [float(x) for x in job.orcafile[file_counter].split()[1:]]
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
