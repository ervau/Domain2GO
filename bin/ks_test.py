import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.stats import kstest

cwd = os.getcwd()
input_file_path = "{}Domain2GO/input_data".format(cwd.split("Domain2GO")[0])
output_file_path = "{}Domain2GO/outputs".format(cwd.split("Domain2GO")[0])

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


# PLOTLARI ÇİZDİRME KISIMLARINI YAPSAK MI? 
# def plot_cooccurrence_distributions(original_mapping, n_treshold, , folder_name, s_threshold=1, random_mapping=None):
#     plt.ioff()
#     fig = plt.figure()

#     if random_mapping:
#         random_co_occs = random_mapping[(random_mapping["n"] < n_treshold) & (random_mapping["s"] < s_threshold)]["s"]
#         bins = np.arange(0, 1, 0.025)
#         plt.hist(random_co_occs, bins, alpha=0.5, log=True, edgecolor= "white", facecolor= "#7A72BD", label="Randomized Mapping") 
#         folder_name = "cooccurence_dist_comparison_plots"
#     else: 
#         bins = np.arange(0, 0.2, 0.0075)
        
    
#     original_co_occs = original_mapping[(original_mapping["n"] < n_treshold) & (original_mapping["s"] < s_threshold)]["s"]
#     plt.hist(original_co_occs, bins, alpha=1, log=True, edgecolor= "white", facecolor= "#4AA3A0", label="Original Mapping") 

#     plt.suptitle("Co-occurence Similarity Distribution (n≥" + str(n_treshold)+ ", S<" + str(s_threshold), fontsize=12)
#     plt.legend(loc='upper right')
#     plt.xlabel('Co-occurence similarity bins', fontsize=12)
#     plt.ylabel('# of unique mappings', fontsize=12)
#     fig.set_figwidth(8)
#     fig.set_figheight(5)

#     output_file = "{}/{}/s{}n{}.png".format(output_file_path, folder_name, s_threshold, n_treshold)
#     plt.savefig(output_file)
#     plt.close(fig)

# def plot_original_distributions():
#     folder_name = "cooccurence_dist_plots"
#     output_folder = "{}/{}".format(output_file_path, folder_name)
#     os.mkdir(output_folder)

#     n_values = []


