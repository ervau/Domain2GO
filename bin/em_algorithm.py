import pandas as pd
from math import log
from time import time
import os
from math import exp
from collections import defaultdict
import sys

alpha = 1.0
beta = 1.0
p_random = .001
t0 = 0.0
t1 = 0.0

cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])

output_folder = "{}/em_algorithm".format(output_file_path)
log_path = "{}/EM.log".format(output_folder)
logfile = open(log_path, 'w')

def calculate_initial_theta(mapping, mapping_w_scores):
  subset_mapping = mapping_w_scores[mapping_w_scores["s"] >= 0.1]
  result_df = pd.merge(mapping, subset_mapping, on=("GO", "Interpro"))[["GO", "Uniprot", "Interpro"]]
  result_df["possible"] = result_df.groupby(by=['GO','Interpro'])['GO'].transform('size')
  curated_data_file = "{}/interpro2go_2021.txt".format(input_file_path)
  curated_data = pd.read_csv(curated_data_file, delimiter = " ")
  curated_data.columns = ["GO", "Interpro"]
  reported = result_df.merge(curated_data, how='left', indicator=True)[result_df.merge(curated_data, how='left', indicator=True)["_merge"] == "both"]
  result_df["reported"] = 0
  result_df.loc[reported.index, "reported"] = 1
  result_df["initial_theta"] = (result_df["reported"] + 1) / (result_df["possible"] + 2)

  return result_df


def initialize_theta(theta, M, N, Z, pairs):
  for p,g,i in pairs:
    theta[(g,i)] = M[(g,i)] / (M[(g,i)] + N[(g,i)] + Z[(g,i)])


# C:dict of dicts
# example: C["p1"] # p1 is a protein
# {('i', 'j'): 1} # i and j is a domain-go pair co-occuring on this protein, 1 is the C value of this pair
def initialize_C(C, pairs):
  for p,g,i in pairs:
    C[p] = {}
  for p,g,i in pairs:
    C[p][(g,i)] = 1.0


def initialize_MN(M, N, pairs):
  for p,g,i in pairs:
    M[(g,i)] = 0.0
    N[(g,i)] = 0.0

def define_RS(R, S, M, initial_theta):
  for idx, key in enumerate(M):
      R[key] = initial_theta["reported"][idx]

  for idx, key in enumerate(M):
    S[key] = initial_theta["possible"][idx]

def initialize_ddi2edge(ddi2edge, pairs):
  for p,g,i in pairs:
    ddi2edge[(g,i)].append(p)


def initialize_Q(Q, M):
  for key in M:
      Q.append(key)

def initialize_E(E, Q):
  for i,j in Q:
    E[(i,j)] = 0.0



# M and N are populated based on C value of pairs
def populate_MN(C, M, N):
  for e in C.keys():
    for i,j in C[e]:
      M[(i,j)] += C[e][(i,j)]
      N[(i,j)] += 1 - C[e][(i,j)]


def update_theta(Q, theta, M, N, Z, alpha, beta):
  for i,j in Q:
    theta[(i,j)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)] + Z[(i,j)] + alpha + beta)
    #theta[(i,j)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)]  + beta)
    #theta[(j,i)] = (M[(i,j)] + alpha) / (M[(i,j)] + N[(i,j)]  + beta)
    if theta[(i,j)] < 0.0 or theta[(i,j)] > 1.0:
      print ("ERROR: theta[(%s,%s)] = %.30f in update" % (i, j, theta[(i,j)]))
      sys.exit()


def update_C(C, theta):
    for e in C.keys():
        # First the denominator, P(m <-> n)
        prod = 1.0
        for k,l in C[e]:
            prod *= (1.0 - theta[(k,l)])
        #p = round(1.0 - prod, 4)
        p = 1.0 - prod
        for i,j in C[e]:
            #q =  round(theta[(i,j)], 4) / p
            q =  round(theta[(i,j)] / p, 4)
            C[e][(i,j)] = q


def update_MN(C, M, N):
  # Zero them out
  for e in C.keys():
    for i,j in C[e]:
      M[(i,j)] = 0.0
      N[(i,j)] = 0.0
  # Now count them
  for e in C.keys():
    for i,j in C[e]:
      M[(i,j)] += C[e][(i,j)]
      N[(i,j)] += (1 - C[e][(i,j)])


def initialize_Z(Z, M):
  z_scores_file = "{}/em_z_scores.txt".format(input_file_path)
  z_scores_data = pd.read_csv(z_scores_file, delimiter = " ")
  for idx, key in enumerate(M):
    Z[key] = z_scores_data["Z"][idx]



def calc_likelihood(Q, theta, M, N, Z, alpha, beta):
  l = 0.0
  bottom = 10e-30
  for i,j in Q:
    a = max(pow(theta[(i,j)], M[(i,j)] + alpha), bottom)
    b = max(pow(1-theta[(i,j)], N[(i,j)] + Z[(i,j)] + beta), bottom)
    p = log(a)
    q = log(b)
    l += (p + q)
  return l


def calc_Eij_null(i, j, ddi2edge, C, theta, theta_null, p_random, E):
  for e in ddi2edge[(i,j)]:
    # Numerator
    nprod = 1.0
    for k,l in C[e]:
      nprod *= (1 - theta[(k,l)])
    u = max(1 - nprod, p_random)
    # Denominator
    dprod = 1.0
    for k,l in C[e]:
      #if (k,l) != (i,j) and (k,l) != (j,i):
      dprod *= (1 - theta_null[(k,l)])
    o = max(1 - dprod, p_random)
    E[(i,j)] += log(u / o)


def run_em_algorithm_theta_calc(mapping, mapping_w_scores, only_theta="yes"):
  
  if only_theta == "yes":
    output_file = "{}/em_only_theta_output.txt".format(output_folder)
    output = open(output_file, 'w')
  
  print("Starting EM algorithm: user can check the progress from log file: {}".format(log_path))

  C = {}
  M = {}
  N = {}
  Z = {}
  R = {}
  S = {}
  Q=[]
  theta = {}

  initial_theta = calculate_initial_theta(mapping, mapping_w_scores)
  pairs = list(zip(initial_theta["Uniprot"], initial_theta["GO"],  initial_theta["Interpro"]))
  
  initialize_C(C, pairs)
  initialize_MN(M, N, pairs)
  populate_MN(C, M, N)
  define_RS(R, S, M, initial_theta)
  initialize_Q(Q, M)
  initialize_Z(Z, M)
  initialize_theta(theta, M, N, Z, pairs)

  t0 = time()
  logfile.write("Running EM to calculate theta values (Initial EM)\n")
  logfile.flush()
  L = calc_likelihood(Q, theta, M, N, Z, alpha, beta)
  t1 = time()
  logfile.write("Initial likelihood: %.2f\n" % (L))
  logfile.flush()

  for iter in range(150):
    update_theta(Q, theta, M, N, Z, alpha, beta)
    update_C(C, theta)
    update_MN(C, M, N)
    L = calc_likelihood(Q, theta, M, N, Z, alpha, beta)
    logfile.write("  %d\t%f\n" % (iter, L))
    logfile.flush()
  
  t1 = time()
  logfile.write("Initial EM finished:\t%.2f\n" % (t1 - t0))
  logfile.flush()

  if only_theta == "yes":
    logfile.write("Writing output file\n")
    logfile.flush()
    output.write("GO\tInterpro\tinitial_theta\ttheta\n")
    output.flush()
    for i,j in Q:
      output.write("%s\t%s\t%f\t%f\n" % (i, j, float(R[(i,j)] + alpha) / (S[(i,j)] + alpha + beta), theta[(i,j)]))
      output.flush()
      
  else:
    return initial_theta, pairs, R, S, theta

def run_em_algorithm_e_calc(mapping, mapping_w_scores):
  print("WARNING: Running EM in full mode, user can check the progress from log file: {}. calculating E values can take several hours, you can choose to calculate only theta values by setting the '--em' parameter to 'only_theta', or to skip EM algorithm by setting the '--em' parameter to 'skip' \n".format(log_path))
  initial_theta, pairs, R, S, theta = run_em_algorithm_theta_calc(mapping, mapping_w_scores, only_theta= "no")

  output_file = "{}/em_full_output.txt".format(output_folder)
  output = open(output_file, 'w')

  C = {}
  M = {}
  N = {}
  Z = {}
  E={}
  Q=[]
  ddi2edge=defaultdict(list)

  initialize_C(C, pairs)
  initialize_MN(M, N, pairs)
  populate_MN(C, M, N)
  initialize_Q(Q, M)
  initialize_Z(Z, M)
  initialize_ddi2edge(ddi2edge, pairs)
  initialize_E(E, Q)

  t0 = time()
  logfile.write("Running EM to calculate E values \n")
  logfile.flush()

  for a,b in Q:
    logfile.write("\tRecalculating excluding %s <-> %s\n" % (a,b))
    logfile.flush()
    
    theta_null={}
    initialize_C(C, pairs)
    initialize_MN(M, N, pairs)
    populate_MN(C, M, N)
    initialize_theta(theta_null, M, N, Z, pairs)
    L = calc_likelihood(Q, theta_null, M, N, Z, alpha, beta)
    logfile.write("\t\tInitial likelihood: %.2f\n" % (L))

    for iter in range(10):
      update_theta(Q, theta_null, M, N, Z, alpha, beta)
      theta_null[(a,b)] = p_random
      update_C(C, theta_null)
      update_MN(C, M, N)
      L_curr = calc_likelihood(Q, theta_null, M, N, Z, alpha, beta)
      logfile.write("\t\t%f\n" % (L_curr))
    
    calc_Eij_null(a, b, ddi2edge, C, theta, theta_null, p_random, E)

  t1 = time()
  logfile.write("All parameters evaluated:\t%.2f\n" % (t1 - t0))
  logfile.flush()
  logfile.write("Writing output file\n")
  logfile.flush()

  output.write("GO\tInterpro\tinitial_theta\ttheta\tE\n")
  output.flush()
  for i,j in Q:
    output.write("%s\t%s\t%f\t%f\t%f\n" % (i, j, float(R[(i,j)] + alpha) / (S[(i,j)] + alpha + beta), theta[(i,j)], E[(i,j)]))

  logfile.write("Finished\n")
  logfile.flush()

