from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd
from numpy import mean

def plot_hist(value_list,title0,x0,y0):
    plt.hist(value_list)
    plt.title(title0)
    plt.xlabel(x0)
    plt.ylabel(y0)
    save(title0,'png')

def count_shared(listy, listx):
    n_shared = 0
    for frag0 in listx:
        for fragy in listy:
            if frag0 in fragy or Levenshtein.ratio(frag0,fragy) > 0.9 or fragy in frag0:
            #if frag0 in fragy or fragy in frag0:
                n_shared += 1
                break
    return n_shared

mhc1_matrix = []
mhc2_matrix = []
MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict','r'))
pid_list = MCL_data['pid']['pid']
#the format: Row X, and Column Y is percentage of patient X is shared with patient Y
#mhc1 vs. mhc2
#MCL_data[x]['MHC1_frag'][i]
set2 = set()
for x in pid_list:
    set2 = set2 | set(MCL_data[x]['MHC2_frag'])
list2 = list(set2)
len2 = []
for x in list2:
    len2.append(len(x))
plot_hist(len2,'Distribution of peptide length among MHCII','AA Length','Frequency')
    