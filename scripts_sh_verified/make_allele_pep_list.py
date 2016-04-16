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
    dict0 = dict()
    #print path_IEDB
    for line0 in open(path_IEDB+'mhc_ligand_full.csv','r'):
        line0=line0.rstrip()
        line0=line0.split('"')      
        if len(line0) > 127:
            #print line0[127]
            if line0[127] in dict0:
                #print line0[127]
                dict0[line0[127]].append(line0[23])
            else:
                dict0[line0[127]] = [line0[23]]
    return dict0

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
        
#0 - negative training 1 - positive training 2 - negative validation 3 - positive validation        
def make_training(path_save,hla_name0,list_len,cluster_list,version0,t_ratio,v_ratio):
    one_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/human_proteinome_oneline.str'
    onegenestr = pickle.load(open(one_gene_path,'r'))
    len_one = len(onegenestr)
    file_out = open(path_save+hla_name0+version0+'_tr_'+str(t_ratio)+'_val.csv','w+')
    train_list = []
    val_list = []
    val_goal = list_len*v_ratio
    #split training and validation data based on clusters
    cluster_num_train = 0
    cluster_num_val = 0
    cluster_list = shuffle_list(cluster_list)
    for cluster0 in cluster_list:
        if len(val_list) < val_goal:
            cluster_num_val += 1
            val_list.extend(cluster0)
        else:
            cluster_num_train += 1
            train_list.extend(cluster0)
    #report training and validation split
    print(hla_name0+' training cluster_n: '+str(cluster_num_train)+' validation cluster_n: '+str(cluster_num_val))
    print(hla_name0+' training pep_n: '+str(len(train_list))+' validation pep_n: '+str(len(val_list)))
    #generate negative for each positive peptides in both training and validation
    for pos0 in val_list:
        #making validation
        file_out.write(pos0+'\t'+'3\n')
        rand0 = random.randint(0,len_one)
        neg0 = onegenestr[rand0:rand0+len(pos0)]
        file_out.write(neg0+'\t'+'2\n')
        neg0 = ''.join(random.sample(pos0,len(pos0)))
        file_out.write(neg0+'\t'+'2\n')
    for pos0 in train_list:
        file_out.write(pos0+'\t'+'1\n')
        for i in range(0,t_ratio):
            rand0 = random.randint(0,len_one)
            neg0 = onegenestr[rand0:rand0+len(pos0)]
            file_out.write(neg0+'\t'+'0\n')
            neg0 = ''.join(random.sample(pos0,len(pos0)))
            file_out.write(neg0+'\t'+'0\n')
    file_out.close()
        
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


####make MCL peptide dictionary
MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
MCL_data = del_sameHLA(MCL_data)
pid_list = MCL_data['pid']['pid']
dict_hla_pid = get_pid_for_each_hla(MCL_data,pid_list)
print_d_list(dict_hla_pid) 
dict_hla_pep = get_pep_for_each_hla(dict_hla_pid,MCL_data)
print_d_list(dict_hla_pep) 
pickle.dump(dict_hla_pep,open(path0+'MCL_pep_by_allele.dict','w+'))

###make IEDB peptide dictionary
dict_IEDB = get_IEDB_pep_dict()
pickle.dump(dict_IEDB,open(path0+'IEDB_pep_by_allele.dict','w+'))


