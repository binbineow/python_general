from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd

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
for mhc0 in ['MHC1','MHC2']:
    #row
    matrix0 = []
    for x in pid_list: 
        # column
        print x
        n0 = len(MCL_data[x][mhc0+'_frag'])
        row0 = []
        for y in pid_list:
            #y on x, put y first in the function
            n_shared = count_shared(MCL_data[y][mhc0+'_frag'],MCL_data[x][mhc0+'_frag'])
            row0.append(float(n_shared)/n0)      
        matrix0.append(row0)
    if mhc0 == 'MHC1':
        mhc1_matrix = matrix0
    else:
        mhc2_matrix = matrix0
print('MHC1')
print(mhc1_matrix)
print('MHC2')
print(mhc2_matrix)
pickle.dump(mhc1_matrix,open('heatmap_mhc1_L.mat','w+'))
pickle.dump(mhc2_matrix,open('heatmap_mhc2_L.mat','w+'))
    
    
    
    

             