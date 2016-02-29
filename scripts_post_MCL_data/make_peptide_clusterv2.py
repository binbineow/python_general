from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd
from numpy import mean
from numpy import std

def plot_hist(value_list,title0,x0,y0,control0=False):
    plt.hist(value_list,25)
    plt.title(title0)
    plt.xlabel(x0)
    plt.ylabel(y0)
    if control0:
        plt.axis([0,12,0,1800])
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
print('total classII peptides='+str(len(list2)))
print('total classII pep clusterings='+str(len(cluster_list)))
n0 = 0
for x in cluster_list:
    if len(x)>1:
        n0 += 1
print('total classII pep Non-singular clusterings='+str(n0))
#print(cluster_list)
print('Average of string mean difference among all clusters'+str(mean(len_d_ave)))
plot_hist(len_d_ave,'Mean string length difference among all Class II peptide clusters','Mean difference in one cluster','Cluster number')
print('Average of string max difference among all clusters'+str(mean(len_d_max)))
plot_hist(len_d_max,'Max string length difference among all Class II peptide clusters','Max difference in one cluster','Cluster number')
print('Average of string length std among all clusters'+str(mean(len_std)))
plot_hist(len_std,'String length standard deviations among all Class II peptide clusters','Standard deviation in one cluster','Cluster number')
#print(len_d_ave)
#print(len_d_max)
        

    
    
    

             