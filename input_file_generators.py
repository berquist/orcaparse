#!/usr/bin/env python2

class EPRFile:
    """
    """
    def __init__(self, charge=0, xyzfile=""):
        pass

    def default(self, charge, xyzfile):
            """
            Generate the ORCA input file, here for a geometry optimization,
            followed by an EPR calculation (g-tensor and N/Cu hyperfine)
            """
            return """! uks def2-tzvpp def2-tzvpp/j ri rijcosx somf(1x) tightscf tightopt grid5

* xyzfile {0} 2 {1}.xyz *

%pal
 nprocs 16
 end

%scf
 maxiter 512
 guess pmodel
 end

%geom
 maxiter 512
 trust 0.3
 inhess almloef
 end

%method
 functional pbe0
 end

%eprnmr
 tol 1e-10
 gtensor 1
 ori -3
 nuclei = all N  {{ aiso, adip, aorb, fgrad, rho }}
 nuclei = all Cu {{ aiso, adip, aorb, fgrad, rho }}
 printlevel 5
 end

""".format(charge, xyzfile)

class OptFile:
    """
    """
    def __init__(self):
        pass

class EPROptFile:
    """
    """
    def __init__(self):
        pass

class PBSFile:
    """
    """
    def __init__(self):
        pass

    def default(self, xyzfile):
        """
        """
        return """#!/bin/bash

#PBS -N {0}
#PBS -q ishared_large
#PBS -l nodes=1:ppn=16
#PBS -l walltime=144:00:00
#PBS -j oe
#PBS -l qos=low

module purge
module load intel/2013.0
module load openmpi/1.6.3-intel13
module load orca/3.0.0

cp $PBS_O_WORKDIR/{0}.inp $LOCAL
cp $PBS_O_WORKDIR/{0}.xyz $LOCAL
cd $LOCAL

run_on_exit() {{
    set -v
    cp $LOCAL/* $PBS_O_WORKDIR
}}

trap run_on_exit EXIT

`which orca` {0}.inp >& $PBS_O_WORKDIR/{0}.out
""".format(xyzfile)

if __name__ == "__main__":
    import argparse
    import subprocess as sp
    import os

    parser = argparse.ArgumentParser(description="Script to output default ORCA input files and PBS/Torque job files to run the corresponding calculations.")
    parser.add_argument(dest="xnamelist", metavar="<xyzfile>", type=str, nargs="+", help="XYZ files; create one job per file")
    parser.add_argument("--charge", dest="charge", metavar="<charge>", type=int, default=0, help="charge on the system")
    args = parser.parse_args()

    xnamelist = args.xnamelist
    charge = args.charge

    for xname in xnamelist:
        stub = os.path.splitext(xname)[0]
        orcahandle = stub + ".inp"
        jobhandle  = stub + ".pbs"
        orcafile = open(orcahandle, "w")
        jobfile  = open(jobhandle, "w")
        
        print >> orcafile, EPRFile.default(charge, stub)
        print >> jobfile,  PBSFile.default(stub)
        
        orcafile.close()
        jobfile.close()


