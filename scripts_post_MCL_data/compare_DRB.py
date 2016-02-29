from utilities import *
#dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd

def get_HLA_name(str0):
    str0 = str0.split('-')[1]
    if ':' in str0:
        str1 = str0.split(':')[0]+':'+str0.split(':')[1]
    else:
        str1 = str0
    #print str1
    return str1

def count_shared(listx, listy):
    x1 = dict_DRB[get_HLA_name(listx[0])]
    x2 = dict_DRB[get_HLA_name(listx[1])]
    y1 = dict_DRB[get_HLA_name(listy[0])]
    y2 = dict_DRB[get_HLA_name(listy[1])]
    dist1 = Levenshtein.ratio(x1,y1) + Levenshtein.ratio(x2,y2)
    dist2 = Levenshtein.ratio(x1,y2) + Levenshtein.ratio(x2,y1)
    dist0 = max(dist1,dist2)/float(2)    
    return dist0

matrix0 = []
MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict','r'))
dict_DRB = pickle.load(open('DRB1_seq.dict','r'))
#print dict_DRB
pid_list = MCL_data['pid']['pid']
#the format: Row X, and Column Y is percentage of patient X is shared with patient Y
#mhc1 vs. mhc2
#MCL_data[x]['MHC1_frag'][i]
print(pid_list)
for x in pid_list: 
        # column
    print x
    #n0 = len(MCL_data[x][mhc0+'_frag'])
    row0 = []
    for y in pid_list:
            #y on x, put y first in the function
        shared0 = count_shared(MCL_data[x]['HLA_typing'][-2:],MCL_data[y]['HLA_typing'][-2:])
        print ([x,y,shared0])
        row0.append(shared0)      
    matrix0.append(row0)
    

#print(matrix0)
pickle.dump(matrix0,open('heatmap_drb.mat','w+'))

    
    
    
    

             