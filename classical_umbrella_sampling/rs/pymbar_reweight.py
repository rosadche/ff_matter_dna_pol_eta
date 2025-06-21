#===============================================================================
#                                   IMPORTS
#===============================================================================
import pymbar  # multistate Bennett acceptance ratio
from pymbar import timeseries  # timeseries analysis

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import os

#===============================================================================
#                                   IMPORTS
#===============================================================================
def amber_umbrella_bias(cv_val, cv_details):
    """
    """
    r1, r2, r3, r4, k1, k2 = cv_details
    
    if cv_val <= r1:
        U = k1 * (r1 - r2) * cv_val
    elif r1 < cv_val <= r2:
        U = 0.5 * k1 * np.power((cv_val - r2), 2)
    elif r2 < cv_val <= r3:
        U = 0
    elif r3 < cv_val <= r4:
        U = 0.5 * k2 * np.power((cv_val - r3), 2)
    elif r4 <= cv_val:
        U = k2 * (r4 - r3) * cv_val
    else: 
        "Error, not possible"
        exit()
    
    #print(r1, r2, r3, r4, k1, k2, U)
    return U


def read_amber_umbrella_1d(filename):
    """
    returns cv_vals as an array where each column is that CV value over time
    returns an array where each column represents cv_number and rows are 
    r1, r2, r3, r4, k1, k2
    """
    # = NFE%PMD ==================================================================
    #   << anchor(1) : position = -98.000000, 2.000000, 2.000000, 102.000000
    #                  strength = 250.000000, 250.000000 >>
    # -----------------------------------------------------------------------------
    # MD time (ps), CV(1:1)
    # =============================================================================

    with open(filename, "r") as f:
        cv_details = []
        
        line = f.readline().rstrip("\n")
        while line[0] == "#":
            if "anchor" in line:
                line = ''.join(c for c in line if (c.isalnum() or c.isspace()) or (c in [".", "-"]))
                split = [x for x in line.split(" ") if x]
                cv_details += [float(x) for x in split[2:6] ]
            if "strength" in line:
                line = ''.join(c for c in line if (c.isalnum() or c.isspace()) or (c in [".", "-"]))
                split = [x for x in line.split(" ") if x]
                cv_details += [float(x) for x in split[1:]  ]
                
            line = f.readline().rstrip("\n")
                
            
    data = np.loadtxt(filename)
    cv_vals = data[:, 1]
    
    cv_details = np.asarray(cv_details)
    return cv_vals, cv_details
    

#===============================================================================
#                             SETUP FILES & MBAR
#===============================================================================
#------ File Details --------------------------------------------------------------
base_dir = "/expanse/lustre/scratch/rosadche/temp_project/dnap_string/wt_2mg_m1264/comdist1/"

K = 109                                  # Represents number of umbrella windows

umbrella_file = base_dir + "image_{}/umbrella_{}.txt"
files = [umbrella_file.format(x, x) for x in range(1, K+1) ]


#------ Constants --------------------------------------------------------------
kB = 1.3806503 * 6.0221415 / 4184.0     # Boltzmann constant in kcal/mol/K
temperature = 303.15                     # All windows run at a single temperature
T_k = np.ones(K, float) * temperature   # initial temperatures are all equal
beta_k = 1.0 / (kB * T_k)               # inverse temperature of simulations in (kcal/mol)^-1
kBT = kB * temperature

#------ Parameters --------------------------------------------------------------
N_max = 2000         # maximum number of snapshots/simulation
nbins = 200        # number of bins for 1D free energy profile

#------ Allocate Storage --------------------------------------------------------------
N_k = np.zeros([K], dtype=int)      # N_k[k] is the number of snapshots from umbrella simulation k
restraint_k = np.zeros([K,6])         # Retraint_k[k] is the Umbrella spring constant and center vals for simualtion k: r1, k1
cv_mat_kn = np.zeros([K, N_max])    # cv_mat_kn[k,n] is the CV value for snapshot n from umbrella simulation k
cv_mat_kn[:] = np.nan
u_kn = np.zeros([K, N_max])         # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k (only needed if T is not constant)
u_kln = np.zeros([K, K, N_max])     # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l
g_k = np.zeros([K])                 # statistical inefficiency of simualtion k 

#------ Load Data --------------------------------------------------------------
print("Loading Data...")
print("Image |   g   |   frames   ")
for k in range(K):
    
    file = files[k]
    cv_vals, cv_details = read_amber_umbrella_1d(file)
    cv_mat_kn[k, :] = cv_vals
    restraint_k[k, :]  = cv_details
    
    # Compute correlation times for cv_val timeseries
    g = timeseries.statistical_inefficiency(cv_vals) # compute statistical inefficiency
    start = 0
    g = 1
    indices = timeseries.subsample_correlated_data(cv_vals[start:], g=g) # compute indices of uncorrelated timeseries
    print( f"{k}     | {g:.2f} | {len(indices)}" )
    
    # Subsample data.
    N_k[k] = len(indices)
    u_kn[k, 0 : N_k[k]]         = u_kn[k, start:][indices]
    cv_mat_kn[k, 0 : N_k[k]]    = cv_vals[start:][indices]

cv_min = np.nanmin(cv_mat_kn)
cv_max = np.nanmax(cv_mat_kn)

#------ Set Up For Final MBAR --------------------------------------------------------------
N_max = np.max(N_k)  # shorten the array size

# Set zero of u_kn -- this is arbitrary.
u_kn -= u_kn.min()

# compute bin centers
bin_center_i = np.zeros([nbins])
bin_edges = np.linspace(cv_min, cv_max, nbins + 1)
for i in range(nbins):
    bin_center_i[i] = 0.5 * (bin_edges[i] + bin_edges[i + 1])

N = np.sum(N_k)
cv_n = pymbar.utils.kn_to_n(cv_mat_kn, N_k=N_k)

# Evaluate reduced energies in all umbrellas
print("Evaluating reduced potential energies...")
for k in range(K):
    for n in range(N_k[k]):
        for l in range(K):
            # Compute energy of snapshot n from simulation k in umbrella potential l
            u_kln[k, l, n] = u_kn[k, n] + beta_k[k] * amber_umbrella_bias(cv_mat_kn[k,n], restraint_k[l])

# initialize free energy profile with the data collected
print("Creating FES Object...")
fes = pymbar.FES(u_kln, N_k, verbose=False)

# Compute free energy profile in unbiased potential (in units of kT).
histogram_parameters = {}
histogram_parameters["bin_edges"] = bin_edges

print("Generating FES...")
fes.generate_fes(u_kn, cv_n, fes_type="histogram", histogram_parameters=histogram_parameters)
results = fes.get_fes(bin_center_i, reference_point="from-lowest", uncertainty_method="analytical")
center_f_i  = kBT * results["f_i"]
center_df_i = kBT * results["df_i"]

# Write out free energy profile
text = "#free energy profile (kcal/mol), from histogramming" + "\n"
text += f"{'bin':>8s} {'f':>8s} {'df':>8s}" + "\n"
for i in range(nbins):
    text += f"{bin_center_i[i]:8.3f} {center_f_i[i]:8.3f} {center_df_i[i]:8.3f}" + "\n"

print(text)

with open("fes.dat", "w") as f:
    f.write(text)
