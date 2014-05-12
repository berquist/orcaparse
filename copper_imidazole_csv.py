#!/usr/bin/env python2

import orca_parser
from copper_imidazole_analysis import CopperImidazoleAnalysis
import argparse
import csv

cia = CopperImidazoleAnalysis()

parser = argparse.ArgumentParser(description="Given pathnames of ORCA output files, make a dump of certain spectroscopic parameters to a CSV file.")

parser.add_argument("--csvname", dest="csvname", metavar="<CSV output root name>", type=str, default="output.csv", help="optional name for the CSV output file")
parser.add_argument(dest="namelist", metavar="<ORCA filename>", nargs="+", type=str, default=None, help="ORCA output files")

args = parser.parse_args()
namelist = args.namelist

with open(args.csvname, 'wb') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')

    csvwriter.writerow(["g-tensor",
                        "id_copper",
                        "A_copper",
                        "euler_copper",
                        "nqcc_copper",
                        "id_nitrogen_far",
                        "A_nitrogen_far",
                        "euler_nitrogen_far",
                        "nqcc_nitrogen_far",
                        "filename"])

    for name in namelist:
        orcafile = orca_parser.ORCAOutputParser(name)
        gtensor, giso = orcafile.return_gtensor()
        id_cu = cia.copper_id(orcafile)
        id_far = cia.nitrogen_far_id(orcafile)
        atensor_cu = cia.hyperfine(orcafile, id_cu)
        atensor_far = cia.hyperfine(orcafile, id_far)
        euler_cu = cia.euler(orcafile, id_cu)
        euler_far = cia.euler(orcafile, id_far)
        nqi_cu, nqcc_cu, eta_cu = cia.nqi(orcafile, id_cu)
        nqi_far, nqcc_far, eta_far = cia.nqi(orcafile, id_far)

        csvwriter.writerow([gtensor,
                            id_cu,
                            atensor_cu,
                            euler_cu,
                            nqcc_cu,
                            id_far,
                            atensor_far,
                            euler_far,
                            nqcc_far,
                            name])
