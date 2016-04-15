from utilities import *
#dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
path_MCL = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
path_IEDB = '/scratch/users/bchen45/HLA_prediction/IEDB/raw_data/'
path_train = '/scratch/users/bchen45/HLA_prediction/RNN_data/core/'
import math
import pandas as pd
from numpy import mean
from numpy import std
target_mhc = 'HLA-DRB1'


###########get peptides from IEDB ############
def get_IEDB_pep_set():
    #1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
    list0 = []
    print path_IEDB
    for line0 in open(path_IEDB+'mhc_ligand_full.csv','r'):
        line0=line0.rstrip()
        line0=line0.split('"')      
        if len(line0) > 127:
            #print line0[127]
            if target_mhc in line0[127]:
                #print line0[127]
                list0.append(line0[23])
    return set(list0)

def get_shared(listy, listx):
    n_shared = 0
    common0=set()
    for fragx in listx:
        for fragy in listy:
            if fragx in fragy or fragy in fragx:
                common0.add(fragx)
                common0.add(fragy)
    return common0

##########get peptides from MCL##############
MCL_data = pickle.load(open(path_MCL+'MCL_data11_18_2015v1.1.dict','r'))
pid_list = MCL_data['pid']['pid']
#the format: Row X, and Column Y is percentage of patient X is shared with patient Y
#mhc1 vs. mhc2
#MCL_data[x]['MHC1_frag'][i]
set2 = set()
for x in pid_list:
    set2 = set2 | set(MCL_data[x]['MHC2_frag'])
set_iedb = get_IEDB_pep_set()
print('IEDB peptide number='+str(len(set_iedb)))
print('MCL peptide number='+str(len(set2)))
overlap0 = get_shared(list(set2),list(set_iedb))
print('Overlap number='+str(len(overlap0)))
#set2 =set3  #stop here
set2 = set2 | set_iedb
list2 = list(set2)
len2 = []
for x in list2:
    len2.append(len(x))
list_len_2 = zip(list2,len2)
list2,len2 = zip(*sorted(list_len_2,key=lambda x: x[1]))
#print list2
list2_tested = [True]*len(list2)
cluster_list = []
len_d_ave = []
len_d_max = []
len_std = []
for n0 in range(0,len(list2)):
    if list2_tested[n0]:
        list0 = []
        len_d0 = []
        len_l0 = []
        list0.append(list2[n0])
        for m0 in range(n0+1,len(list2)):
            if list2_tested[m0]:
                if list2[n0] in list2[m0]:
                    list2_tested[m0] = False
                    list0.append(list2[m0])
                    len_d0.append(len(list2[m0])-len(list2[n0]))
        if len(len_d0) >= 1:
            len_d_ave.append(mean(len_d0))
            len_d_max.append(max(len_d0))
            len_d0.append(0)
            len_std.append(std(len_d0))
        cluster_list.append(list0)
core_list = []
for list0 in cluster_list:
    core_list.append(list0[0])
pickle.dump(cluster_list,open(path_train+'drb1_MCL_IEDB_cluster.list','w+'))
pickle.dump(core_list,open(path_train+'drb1_MCL_IEDB_core.list','w+'))
print(str(len(core_list)))
    
    
    

             