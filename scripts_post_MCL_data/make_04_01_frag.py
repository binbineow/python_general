from utilities import *
#dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd

def get_shared(listy, listx):
    n_shared = 0
    common0=set()
    for fragx in listx:
        for fragy in listy:
            if fragx in fragy or fragy in fragx:
                common0.add(fragx)
                common0.add(fragy)
    return common0


MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict','r'))
#pid_list = MCL_data['pid']['pid']
#MCL005,008,0022,041
pid_04_01 = ['MCL005','MCL008','MCL022','MCL041']
set_04_01=set()
for n0 in range(0,len(pid_04_01)):
    for m0 in range(n0+1,len(pid_04_01)):
        print([pid_04_01[n0],pid_04_01[m0]])
        common0=get_shared(MCL_data[pid_04_01[n0]]['MHC2_frag'],MCL_data[pid_04_01[m0]]['MHC2_frag'])
        set_04_01=set_04_01 | common0
#print(len(set_04_01))

############make clusters
list2 = list(set_04_01)
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
        cluster_list.append(list0)
print('total classII peptides='+str(len(list2)))
print('total classII pep clusterings='+str(len(cluster_list)))
n0 = 0
for x in cluster_list:
    if len(x)>1:
        n0 += 1
print('total classII pep Non-singular clusterings='+str(n0))
print cluster_list
  