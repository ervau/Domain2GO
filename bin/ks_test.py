import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import kstest

cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/output".format(cwd.split("Domain2GO")[0])

def ks_test_calc(original, random):
    # score combinations for s and n values
    s_range=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    n_range=[1, 2, 3, 4, 5]

    column_names= ["s<" + str(s) for s in s_range]
    index_names=["n>=" + str(s) for s in n_range]
    ks_test_results=pd.DataFrame(columns = column_names, index=index_names)

    # creating histograms with 20 bins for each score (s & n) combination
    # then performing ks test with the value distribution of the histogram bins
    # comparing original and random co-occurence similarity distributions

    for x, y in [(x,y) for x in s_range for y in n_range]:
        original_subset = original.loc[(original["n"] >= y ) & (original["s"] < x )]
        random_subset = random.loc[(random["n"] >= y ) & (random["s"] < x )]
        original_bin_values = np.histogram(original_subset["s"], bins=np.arange(0, x+0.001, 0.005))[0]
        random_bin_values = np.histogram(random_subset["s"], bins=np.arange(0, x+0.001, 0.005))[0] / 10 # dividing to 10 because 10 different random mapping sets were created to avoid bias

        ks_test_result = kstest(original_bin_values, random_bin_values)[1]
        ks_test_results.at["n>=" + str(y) , "s<" + str(x)] = round(ks_test_result, 10)

    
    output_file = "{}/ks_test_results.txt".format(output_file_path)
    ks_test_results.to_csv(output_file)

    return output_file

