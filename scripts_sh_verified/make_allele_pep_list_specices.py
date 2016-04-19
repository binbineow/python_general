#make a dictionary = HLA-type to list of peptides
#code to convert a list into cluster 
from utilities import *
#dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
path0 = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
import math
import random
import re

#not allow substring or meta-string to be added into validation set
def check_val(pep_train, pep0):
    return0 = True
    for pep_t0 in pep_train0:
        if pep0 in pep_t0 or pep_t0 in pep0:
            return0 = False
            break
    return return0

#make a list of lists/culsters, within each cluster, strings are substring among each other within a cluster
def make_cluster(list0):
    list2 = list(list0)
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
    return cluster_list


#return a dictionary of d' which is d without key
def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

##substring between two patients
def get_shared(listy, listx):
    n_shared = 0
    common0=set()
    for fragx in listx:
        for fragy in listy:
            if fragx in fragy or fragy in fragx:
                common0.add(fragx)
                common0.add(fragy)
    return common0

#make a dictionary = HLA-type to list of patient IDs
def get_pid_for_each_hla(MCL_data,pid_list):
    dict_hla_pid = dumb()
    for pid0 in pid_list:
        #print MCL_data[pid0]['HLA_typing']
        hla1 = MCL_data[pid0]['HLA_typing'][-1]
        hla2 = MCL_data[pid0]['HLA_typing'][-2]
        dict_hla_pid[hla1].append(pid0)
        dict_hla_pid[hla2].append(pid0)
    return dict_hla_pid

#get IEDB peptides
def get_IEDB_pep_dict():
    #1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
    list_IEDB= []
    list_IEDB_MCL = []
    list_human = []
    list_IEDB_MCL_human = []
    dict_s = dict()
    human_c = 0
    path_IEDB = '/scratch/users/bchen45/HLA_prediction/IEDB/raw_data/'
    #print path_IEDB
    for line0 in open(path_IEDB+'mhc_ligand_full.csv','r'):
        line0=line0.rstrip()
        line0=line0.split('"')      
        if len(line0) > 127:
            #print line0[127]
            #process the HLA allele name
            if 'HLA-DRB1' in line0[127] and (not '/' in line0 [127]) and (not '(' in line0[23]) :
                line0[127] = line0[127][0:len('HLA-DRB1*08:02')+1]
                if line0[127] in list_MCL_hla:
                    list_IEDB_MCL.append(line0[23])
                dict_s[line0[23]] = line0[37]
                list_IEDB.append(line0[23])
                if not line0[37] == '9606':
                    human_c +=1
                    list_human.append(line0[23])
                    if line0[127] in list_MCL_hla:
                        list_IEDB_MCL_human.append(line0[23])
                       
    print('Total human='+str(human_c))           
    return [set(list_IEDB),set(list_human), set(list_IEDB_MCL), set(list_IEDB_MCL_human), dict_s]

#make peptide list from a set of pids
def get_pep_for_each_hla(dict_hla_pid,MCL_data):
    dict_hla_pep = dumb()
    for key, value in dict_hla_pid.iteritems():
        print key
        if len(value) > 1:
            set0 = set()
            for n0 in range(0,len(value)):
                for m0 in range(n0+1,len(value)):
                    common0 = get_shared(MCL_data[value[n0]]['MHC2_frag'],MCL_data[value[m0]]['MHC2_frag'])
                    set0 = set0 | common0
            #consider to add a cut-off
            dict_hla_pep[key] = list(set0)
    return dict_hla_pep

def print_d_list(dict0):
    for key,value in dict0.iteritems():
        print(key+': '+str(len(value)))

def shuffle_list(list0):
    list0 = random.sample(list0,len(list0))
    return list0
        

#this process merge patient MHC2_frag into the first patient scanned and delete the redundant pid
def del_sameHLA(MCL_data0):
    pid = MCL_data0['pid']['pid']
    to_del0 = set()
    for i in range(0,len(pid)):
        pid0 = pid[i]
        hla1 = MCL_data0[pid0]['HLA_typing'][-1]
        hla2 = MCL_data0[pid0]['HLA_typing'][-2]
        for j in range(i+1,len(pid)):
            pid1 = pid[j]
            hla1_1 = MCL_data0[pid1]['HLA_typing'][-1]
            hla2_1 = MCL_data0[pid1]['HLA_typing'][-2]
            set0 = set()
            set0.add(hla1)
            set0.add(hla2)
            set1 = set()
            set1.add(hla1_1)
            set1.add(hla2_1)
            if set0 == set1:
                MCL_data0[pid0]['MHC2_frag'].extend(MCL_data[pid1]['MHC2_frag'])
                to_del0.add(pid1)
    print(to_del0)
    for x in to_del0:
        del MCL_data0[x]
    MCL_data0['pid']['pid'] = list(set(MCL_data0['pid']['pid'])- to_del0)
    return MCL_data0

def get_shared_v2(str0,listy, listx):
    n_shared = 0
    common0=set()
    for fragx in listx:
        for fragy in listy:
            if fragx in fragy or fragy in fragx:
                common0.add(fragx)
                common0.add(fragy)
    print(str0+' N='+str(len(common0)))
    print(str0+' Percentage of 2nd List='+str(float(len(common0))/len(listx)))
    return common0

def tell_me_length(str0,list0):
    print(str0+' N='+str(len(list0)))

####make MCL peptide dictionary
# MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
# MCL_data = del_sameHLA(MCL_data)
# pid_list = MCL_data['pid']['pid']
# dict_hla_pid = get_pid_for_each_hla(MCL_data,pid_list)
# print_d_list(dict_hla_pid) 
# dict_hla_pep = get_pep_for_each_hla(dict_hla_pid,MCL_data)
# print_d_list(dict_hla_pep) 
# pickle.dump(dict_hla_pep,open(path0+'MCL_pep_by_DRB1.dict','w+'))

#####set up MCL peptides
dict_MCL = pickle.load(open(path0+'MCL_pep_by_DRB1.dict','r'))
list_MCL_hla = []
MCL_pep = set()
for key,value in dict_MCL.iteritems():
    list_MCL_hla.append(key)
    MCL_pep = MCL_pep | set(value)
#MCL_pep = list(MCL_pep)
###MCL_pep now is a set containing all peptides from our patients
###list_MCL_hla contains all HLA types in our patients
###THis program will:
##create a dictionary from sequence to organism specices ID
##what percentage of human peptides in IEDB appear in our MCL peptides
##what percenate of human peptides with the correct alleles appear in our MCL peptides


###make IEDB peptide dictionary
[set_IEDB, set_human, set_IEDB_MCL, set_IEDB_MCL_human, dict_IEDB_s] = get_IEDB_pep_dict()
tell_me_length('All MCL DRB1 peptides',MCL_pep)
tell_me_length('All IEDB DRB1 peptides',set_IEDB)
tell_me_length('All IEDB DRB1 with human peptides',set_human)
tell_me_length('All IEDB with alleles in our MCL patients (DRB1)',set_IEDB_MCL)
tell_me_length('All IEDB with human peptides plus MCL-specific alleles',set_IEDB_MCL_human)

#general common = not restricted to human
#human common = restricted to common
#allele common = restricted to specific alleles
general_common = get_shared_v2('general_common',MCL_pep,set_IEDB)
human_common = get_shared_v2('human_common',MCL_pep,set_human)
allele_common = get_shared_v2('Allele_common',MCL_pep,set_IEDB_MCL)
human_allele_common = get_shared_v2('human+allele',MCL_pep,set_IEDB_MCL_human)

