from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd

def count_shared(listy, listx):
    n_shared = 0
    common0 = set()
    for fragx in listx:
        for fragy in listy:
            if fragx in fragy or fragy in fragx:
                common0.add(fragx)
                common0.add(fragy)
    x_only = set(listx) - common0
    y_only = set(listy) - common0
    jaccard0 = float(len(common0))/(len(common0)+len(x_only)+len(y_only))
    return jaccard0

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
            jaccard0 = count_shared(MCL_data[y][mhc0+'_frag'],MCL_data[x][mhc0+'_frag'])
            row0.append(jaccard0)      
        matrix0.append(row0)
    if mhc0 == 'MHC1':
        mhc1_matrix = matrix0
    else:
        mhc2_matrix = matrix0
print('MHC1')
print(mhc1_matrix)
print('MHC2')
print(mhc2_matrix)
pickle.dump(mhc1_matrix,open('heatmap_mhc1_J.mat','w+'))
pickle.dump(mhc2_matrix,open('heatmap_mhc2_J.mat','w+'))
    
    
    
    

             