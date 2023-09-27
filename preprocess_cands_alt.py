# A script to load in and preprocess candidates using ACCEL_sift machinery
# Remove duplicates for whole list, track them, and print out things about the hits of the points left
# Collate a candlist

# Processing and saving utilities
import pickle
import numpy as np
import pandas as pd
from optparse import OptionParser
import os, sys
import re

#from presto.presto import fourierprops, get_rzw_cand                     # for this, initialize and pass a fourierprops object for cand to read binary file
import carmen_sifting as sifting
import copy

from multiprocessing import Pool
import multiprocessing

from argparse import ArgumentParser

from datetime import datetime

import traceback

def unpack_to_pandas(candlist):

    cand_frame = pd.DataFrame(columns=["cand_num", "sigma", "sum_pow", "c_pow", 
                                       "num_harm", "period", "frequency", 
                                       "rbin", "zbin", "path", "DMstr", "DM", 
                                       "tobs", "note"])

    allcands = candlist.get_all_goodcands()

    cand_frame["cand_num"] = np.array([c.candnum for c in allcands])
    cand_frame["sigma"] = np.array([c.sigma for c in allcands])
    cand_frame["sum_pow"] = np.array([c.ipow_det for c in allcands])
    cand_frame["c_pow"] = np.array([c.cpow for c in allcands])
    cand_frame["num_harm"] = np.array([c.numharm for c in allcands])
    cand_frame["period"] = np.array([c.p for c in allcands])
    cand_frame["frequency"] = np.array([c.f for c in allcands])
    cand_frame["rbin"] = np.array([c.r for c in allcands])
    cand_frame["zbin"] = np.array([c.z for c in allcands])
    cand_frame["path"] = np.array([os.path.join(c.path, c.filename) for c in allcands])
    cand_frame["DMstr"] = np.array([c.DMstr for c in allcands])
    cand_frame["DM"] = np.array([c.DM for c in allcands])
    cand_frame["tobs"] = np.array([c.T for c in allcands])
    cand_frame["note"] = np.array([c.note for c in allcands])


    return cand_frame

if __name__ == "__main__":
    usage = "usage: %prog [options]"
    parser = ArgumentParser()

    parser.add_argument("-d", "--dir", type=str, dest="directory", default='.', required=True,
                        help="Directory with candidate files in it")
    parser.add_argument("-m", "--mjd", type=str, dest="mjd", default=None, required=True,
                        help="Text file with MJDs to process")
    parser.add_argument("-o", "--outdir", type=str, dest="outdir", default="/users/cchoza/chunk_search", required=True,
                        help="Directory for output files")

    args = parser.parse_args()

    cand_dir = args.directory
    mjd_file = args.mjd

    prelim_reject = True

    # In how many DMs must a candidate be detected to be considered "good"
    min_num_DMs = 2
    # Lowest DM to consider as a "real" pulsar
    low_DM_cutoff = 2.0
    # Ignore candidates with a sigma (from incoherent power summation) less than this
    sifting.sigma_threshold = 4.0
    # Ignore candidates with a coherent power less than this
    sifting.c_pow_threshold = 100.0

    # The following are all defined in the sifting module.
    # But if we want to override them, uncomment and do it here.
    # You shouldn't need to adjust them for most searches, though.

    # How close a candidate has to be to another candidate to                
    # consider it the same candidate (in Fourier bins)
    sifting.r_err = 1.1
    # Shortest period candidates to consider (s)
    sifting.short_period = 0.0005
    # Longest period candidates to consider (s)
    sifting.long_period = 15.0
    # Ignore any candidates where at least one harmonic does exceed this power
    sifting.harm_pow_cutoff = 8.0

    # If the birds file works well, the following shouldn't
    # be needed at all...  If they are, add tuples with the bad
    # values and their errors.
    #                (ms, err)
    sifting.known_birds_p = []
    #                (Hz, err)
    sifting.known_birds_f = []

    par_dir = '/users/sransom/timing/Ter5'


    par_files = [f for f in os.listdir(par_dir) if '.par' in f]

    lock = multiprocessing.Lock()
    with lock:
        all_cand_files = [os.path.join(cand_dir, f) for f in os.listdir(cand_dir) if '.cand' not in f and '.inf' not in f and 'ACCEL' in f]

    for file in par_files:
        with open(os.path.join(par_dir, file), 'r') as f:
            lines = f.readlines()
            pulsar_name = [line for line in lines if 'PSR' in line][0].split()[1]
            spin_freq = [line for line in lines if 'F0' in line][0].split()[1]
            pulse_period = 1/float(spin_freq)
            if 'ad' not in pulsar_name and 'aw' not in pulsar_name:
                sifting.known_birds_f.append((float(spin_freq), 0.1))
                sifting.known_birds_p.append((pulse_period, 0.0001))
                continue

    curr_time = datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
    run_path = os.path.join("runs", curr_time)
    os.makedirs(run_path)

    def process_day(day):
        day = day.strip()

        out_path = os.path.join(run_path, f"{day}_88M_FIT_DUPLICATES.txt")
        with open(out_path, "w") as f:
            try:
                f.write("Started\n")
                
                cand_files = [os.path.join(cand_dir, f) for f in all_cand_files if day.split()[0] in f]

                # Check to see if this is from a short search
                if len(re.findall("_[0-9][0-9][0-9]M_", cand_files[0])):
                    dmstrs = [x.split("DM")[-1].split("_")[0] for x in cand_files]
                dms = list(map(float, dmstrs))
                dms.sort()
                dmstrs = ["%.2f"%x for x in dms]
                print(day.split()[0])

                f.write("Loading files and removing known pulsars\n")
                
                filtered_cands = []

                for cand_file in cand_files:
                    # Read in candidates
                    cands = sifting.candlist_from_candfile(cand_file, f)
                    filtered_cands.extend(np.copy(cands.get_all_cands()))

                for c in filtered_cands:
                    if os.path.join(cand_dir, c.filename) not in cand_files:
                        f.write(f"WRONG FILE: {c.filename}\n")

                # Create huge candlist for all files with MJD
                urcands = sifting.Candlist(filtered_cands, trackbad=False, trackdupes=True)

                number_of_cands = len(urcands.get_all_goodcands())
                # print(f"NUMBER OF CANDIDATES BEFORE FILTERING: {number_of_cands}", flush=True)
                f.write(f"Number of candidates before filtering: {number_of_cands}\n")

                if prelim_reject:
                        urcands.default_rejection(f)

                number_of_cands = len(urcands.get_all_goodcands())
                # print(f"NUMBER OF CANDIDATES AFTER DEFAULT REJECTION: {number_of_cands}", flush=True)
                f.write(f"Number of candidates after default rejection: {number_of_cands}\n")

                ### ACCEL_sift on the candidates

                # f.write(f"{np.unique(np.array([c.hits[0][0] for c in urcands.get_all_cands()]))}\n")

                # print(f"[REMOVING DUPLICATES] Process {os.getpid()} | Day {day}")
                f.write(f"Removing duplicates\n")

                # Remove duplicates for full list
                if len(urcands):
                    urcands = sifting.remove_duplicate_candidates(urcands, f)

                # f.write(f"{np.unique(np.array([c.hits[0][0] for c in urcands.get_all_cands()]))}\n")
                f.write(f"Removing DM problems\n")
                # f.write(f"{len(cand_files)}\n")

                # Remove DM problems    
                if len(urcands):
                    urcands = sifting.remove_DM_problems(urcands, min_num_DMs, dmstrs, low_DM_cutoff, verbosity=1)

                number_of_cands = len(urcands.get_all_goodcands())
                f.write(f"Number of candidates after removing DM problems: {number_of_cands}\n")

                if len(urcands):
                    urcands = sifting.remove_low_DM_variance(urcands, f, verbosity=1)

                number_of_cands = len(urcands.get_all_goodcands())
                f.write(f"Number of candidates after removing DM problems: {number_of_cands}\n")

                f.write(f"Removing harmonics\n")
                # Remove harmonics
                if len(urcands):
                    urcands = sifting.remove_harmonics(urcands)

                number_of_cands = len(urcands.get_all_goodcands())
                f.write(f"Number of candidates after removing harmonics: {number_of_cands}\n")

                f.write("Pickling\n")
                # Pickle urcands Candlist object
                candlist_outname = os.path.join(args.outdir, f"candlist_{day.split()[0]}_88M_FIT_DUPLICATES.p")
                pickle.dump(urcands, open(candlist_outname, "wb"))

                f.write("Unpacking to pandas\n")
                # Unpack the candlist to pandas
                ur_cand_frame = unpack_to_pandas(urcands)
                # Discard those above 1000 Hz
                ur_cand_frame = ur_cand_frame.loc[ur_cand_frame["frequency"]<1000]
                number_of_cands = len(ur_cand_frame)
                f.write(f"Number of candidates after cutting to 1000 Hz: {number_of_cands}\n")

                f.write("Binning in frequency\n")
                ur_cand_frame["bin_number"] = [-1]*len(ur_cand_frame)
                frequencies = np.sort(ur_cand_frame.frequency.tolist())
                hist, bin_edges = np.histogram(frequencies, bins=10000, range=[0.0, 1000.0])
                labels = range(len(hist))
                ur_cand_frame["bin_number"] = pd.cut(ur_cand_frame["frequency"], bins=bin_edges, labels=labels)

                f.write("Pickling panda\n")

                csv_outname = os.path.join(args.outdir, f"candlist_{day.split()[0]}_88M_FIT_DUPLICATES_panda.p")
                pickle.dump(ur_cand_frame, open(csv_outname, "wb"))

                f.write("Finished\n")

            except Exception as e:
                f.write("\n")
                f.write(str(e))
                f.write(traceback.format_exc())



# Try changing to starmap with [(tuples, of), (days, day_files)]
    with open(mjd_file, 'r') as days_to_process:
        days = days_to_process.readlines()
        with Pool(10) as p:
            p.map(process_day, days)