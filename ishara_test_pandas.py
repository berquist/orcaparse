#!/usr/bin/env python

i1_e = {'g_para': 2.22, 'g_perp': 2.06, 'A_para': 535, 'A_perp': 39, 'A_iso': 1.70, 'T_dip': 0.14, 'kappa': 1.70, 'eta': 0.67, 'count_nitrogens': 5, 'count_imidazoles': 1, 'origin': 'expt'}
i2_e = {'g_para': 2.30, 'g_perp': 2.05, 'A_para': 443, 'A_perp': 42, 'A_iso': 1.75, 'T_dip': 0.12, 'kappa': 1.65, 'eta': 0.71, 'count_nitrogens': 4, 'count_imidazoles': 2, 'origin': 'expt'}
i4_e = {'g_para': 2.26, 'g_perp': 2.06, 'A_para': 510, 'A_perp': 39, 'A_iso': 1.72, 'T_dip': 0.14, 'kappa': 2.80, 'eta': 0.67, 'count_nitrogens': 8, 'count_imidazoles': 4, 'origin': 'expt'}

def determine_origin(stub):
    '''Based on a filename stub, return a string that identifies the origin
    of the structure.'''
    if 'full' in stub:
        return 'qm_opt'
    else:
        return 'qm_crystal'

if __name__ == '__main__':

    import argparse
    import os
    import pandas as pd

    import orcaparse.orca_parser as orca_parser
    from orcaparse.copper_imidazole_analysis import CopperImidazoleAnalysis

    cia = CopperImidazoleAnalysis()

    parser = argparse.ArgumentParser()
    parser.add_argument('namelist', nargs='+')
    args = parser.parse_args()
    namelist = args.namelist

    all_results = dict()

    for filename in namelist:
        # parse the filename for:
        # 1. number of imidazoles (always the first character)
        # 2. where's the structure from? 
        #  'expt', 'qm_crystal' (as-is or constrained opt), or 'qm_opt' (full opt)
        stub = os.path.splitext(os.path.basename(filename))[0]
        count_imidazoles = int(stub[0])

        orcafile = orca_parser.ORCAOutputParser(filename)

        job_results = dict()
        gtensor, giso = orcafile.return_gtensor()
        g_para = gtensor[2]
        g_perp = (gtensor[0] + gtensor[1])/2

        id_copper = cia.copper_id(orcafile)
        ids_nitrogen = []

        atom_copper = orcafile.molecule.atoms[id_copper]
        atoms_nitrogen = []

        # if we are within the distance cutoff, we must be directly bonded, and do not
        # contribute to the ESEEM results, so leave those out
        cutoff = 3.0
        for id_n in cia.nitrogen_ids(orcafile):
            if orcafile.pair_distance(id_copper, id_n) > cutoff:
                ids_nitrogen.append(id_n)
                atoms_nitrogen.append(orcafile.molecule.atoms[id_n])

        A_para = atom_copper.hyperfine.atensor[2]
        A_perp = (atom_copper.hyperfine.atensor[0] + atom_copper.hyperfine.atensor[1])/2

        # need to average A_iso, T_dip, kappa, and eta
        vals_A_iso = []
        vals_T_dip = []
        vals_kappa = []
        vals_eta = []
        for id_n in ids_nitrogen:
            hyp = orcafile.molecule.atoms[id_n].hyperfine
            efg = orcafile.molecule.atoms[id_n].efg

            vals_A_iso.append(hyp.aiso)
            vals_T_dip.append(hyp.tdip)
            vals_kappa.append(efg.nqcc)
            vals_eta.append(efg.eta)

        A_iso = sum(vals_A_iso)/len(vals_A_iso)
        T_dip = sum(vals_T_dip)/len(vals_T_dip)
        kappa = sum(vals_kappa)/len(vals_kappa)
        eta = sum(vals_eta)/len(vals_eta)

        spin_contam_pct = (orcafile.molecule.ssq_deviation / orcafile.molecule.ssq_ideal) * 100
        alpha_sq = cia.determine_copper_covalency(orcafile)

        print("   name:", filename)
        print(" g_para:", g_para)
        print(" g_perp:", g_perp)
        print(" A_para:", A_para)
        print(" A_perp:", A_perp)
        print("========")
        print("  atoms:", ids_nitrogen)
        print("  A_iso:", vals_A_iso, A_iso)
        print("  T_dip:", vals_T_dip, T_dip)
        print("  kappa:", vals_kappa, kappa)
        print("    eta:", vals_eta, eta)
        print("========")
        print(" ssq_pct:", spin_contam_pct)
        print("alpha_sq:", alpha_sq)
        print("\n")

        job_results['filename'] = os.path.abspath(filename)
        job_results['g_para'] = g_para
        job_results['g_perp'] = g_perp
        job_results['A_para'] = abs(A_para)
        job_results['A_perp'] = abs(A_perp)
        job_results['nitrogen_atoms'] = ids_nitrogen
        job_results['A_iso_vals'] = vals_A_iso
        job_results['A_iso'] = A_iso
        job_results['T_dip_vals'] = vals_T_dip
        job_results['T_dip'] = T_dip
        job_results['kappa_vals'] = vals_kappa
        job_results['kappa'] = abs(kappa)
        job_results['eta_vals'] = vals_eta
        job_results['eta'] = eta
        job_results['count_imidazoles'] = count_imidazoles
        job_results['count_nitrogens'] = len(cia.nitrogen_ids(orcafile))
        job_results['origin'] = determine_origin(stub)
        # job_results['basis_copper'] = 
        # job_results['basis_full'] = 
        job_results['spin_contam_pct'] = spin_contam_pct
        job_results['alpha_sq'] = alpha_sq

        all_results[filename] = job_results

    all_results['expt1'] = i1_e
    all_results['expt2'] = i2_e
    all_results['expt4'] = i4_e

    all_results_df = pd.DataFrame(all_results).transpose()
    all_results_df.to_json('results.json')
