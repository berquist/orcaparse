#!/usr/bin/env python2

import numpy as np
import numpy.linalg as npl
import argparse
import orca_parser
import pyqchem

# parser = argparse.ArgumentParser(description="")

# parser.add_argument(dest="orca_filename", metavar="<orca filename", type=str, default=None, help="Path to ORCA output file.")
# parser.add_argument(dest="qchem_filename", metavar="<qchem filename>", type=str, default=None, help="Path to Q-Chem output file.")

# args = parser.parse_args()

# orca_filename = args.orca_filename
# qchem_filename = args.qchem_filename

orca_filename = "/home/eric/Desktop/epr_tests/cpscfman_test_jobs/nmr/chiral.bsl.orca.out"
qchem_filename = "/home/eric/Desktop/epr_tests/cpscfman_test_jobs/nmr/chiral.bsl.qchem.out"

orca_job = orca_parser.ORCAOutputParser(orca_filename)
qchem_job = pyqchem.QChemOutputParser(qchem_filename)

qchem_nmr = pyqchem.nmr.NMR(qchem_job)

# compare results on an atom-by-atom basis

natom = orca_job.molecule.num_atoms()

for atom_num in range(natom):
    orca_shiftmat = orca_job.molecule.atoms[atom_num].nmr.shiftmat
    orca_eigvals_calc = orca_job.molecule.atoms[atom_num].nmr.eigvals
    orca_eigvals_print = orca_job.molecule.atoms[atom_num].nmr.shiftpri
    orca_iso_calc = orca_job.molecule.atoms[atom_num].nmr.iso
    orca_iso_print = orca_job.molecule.atoms[atom_num].nmr.shiftiso
    qchem_shiftmat = qchem_nmr.totshifts[atom_num]
    qchem_eigvals = qchem_nmr.totshift_eigvals[atom_num]
    qchem_iso_calc = qchem_nmr.totshift_iso[atom_num]
    qchem_iso_print = qchem_nmr.traces_iso_tot[atom_num]

    print "======================================================================"
    print "  ORCA atom:", orca_job.molecule.atoms[atom_num]
    print "Q-Chem atom:", qchem_job.atoms[atom_num], qchem_job.positions[atom_num]
    print "  ORCA shiftmat:"
    print orca_shiftmat
    print "Q-Chem shiftmat:"
    print qchem_shiftmat
    print "Absolute difference:"
    print np.absolute(orca_shiftmat - qchem_shiftmat)
    print "  ORCA eigenvalues (calc):"
    print orca_eigvals_calc
    print "  ORCA eigenvalues (print):"
    print orca_eigvals_print
    print "Q-Chem eigenvalues:"
    print qchem_eigvals
    print "Absolute difference (calc, print) and norm:"
    print np.absolute(orca_eigvals_calc - qchem_eigvals), npl.norm(orca_eigvals_calc - qchem_eigvals)
    print np.absolute(orca_eigvals_print - qchem_eigvals), npl.norm(orca_eigvals_print - qchem_eigvals)
    print "  ORCA isotropic (calc, printed): {0:12.8f} {1:12.8f}".format(orca_iso_calc, orca_iso_print)
    print "Q-Chem isotropic (calc, printed): {0:12.8f} {1:12.8f}".format(qchem_iso_calc, qchem_iso_print)
    print "======================================================================\n"

