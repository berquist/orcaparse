#!/usr/bin/env python2

import orca_parser
from copper_imidazole_analysis import CopperImidazoleAnalysis
import argparse
import csv

cia = CopperImidazoleAnalysis()

parser = argparse.ArgumentParser(description="Given pathnames of ORCA output files, make a dump of all nitrogen parameters to a CSV file.")

parser.add_argument("--csvname", dest="csvname", metavar="<CSV output root name>", type=str, default="nitrogen.csv", help="optional name for the CSV output file")
parser.add_argument(dest="namelist", metavar="<ORCA filename>", nargs="+", type=str, default=None, help="ORCA output files")

args = parser.parse_args()
namelist = args.namelist

with open(args.csvname, 'wb') as csvfile:
    csvwriter = csv.writer(csvfile, delimiter=',')

    for name in namelist:
        csvwriter.writerow([name])
        csvwriter.writerow(["g-tensor",
                            "id_copper",
                            "A_copper (MHz)",
                            "euler_copper (deg.)",
                            "NQCC_copper (MHz)",
                            "eta_copper"])

        orcafile = orca_parser.ORCAOutputParser(name)
        gtensor, giso = orcafile.return_gtensor()
        id_copper = cia.copper_id(orcafile)
        atensor_copper = cia.hyperfine(orcafile, id_copper)
        euler_copper = cia.euler(orcafile, id_copper)
        nqi_copper, nqcc_copper, eta_copper = cia.nqi(orcafile, id_copper)

        csvwriter.writerow([gtensor,
                            id_copper,
                            atensor_copper,
                            euler_copper,
                            nqcc_copper,
                            eta_copper])

        csvwriter.writerow(["",
                            "id_nitrogen",
                            "A_nitrogen (MHz)",
                            "euler_nitrogen (deg.)",
                            "NQCC_nitrogen (MHz)",
                            "eta_nitrogen",
                            "Cu_N_distance (Angstroms)"])

        nitrogen_list = orcafile.find_element("N")

        for id_nitrogen in nitrogen_list:
            atensor_nitrogen = cia.hyperfine(orcafile, id_nitrogen)
            euler_nitrogen = cia.euler(orcafile, id_nitrogen)
            nqi_nitrogen, nqcc_nitrogen, eta_nitrogen = cia.nqi(orcafile, id_nitrogen)
            cu_n_dist = orcafile.pair_distance(id_copper, id_nitrogen)

            csvwriter.writerow(["",
                                id_nitrogen,
                                atensor_nitrogen,
                                euler_nitrogen,
                                nqcc_nitrogen,
                                eta_nitrogen,
                                cu_n_dist])
