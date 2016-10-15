#make a dictionary = HLA-type to list of peptides
#code to convert a list into cluster 
from utilities import *
from predict_netmhciipan import *
import math
import random
import re
import subprocess
from collections import defaultdict
import os
from matplotlib.font_manager import pickle_dump

####this script will be run on Sherlock#########
path_save = '/share/PI/rbaltman/bchen45/data/ig_specific/'
path0 = '/share/PI/rbaltman/bchen45/data/ig_specific/'
one_gene_path = '/share/PI/rbaltman/bchen45/software/IEDB/test0/human_proteinome_oneline.str'

def make_neg(pep_list):
    list_out = []
    for x in pep_list:
        neg0 = ''.join(random.sample(x,len(x)))
        rand0 = random.randint(0,len_one)
        neg1= onegenestr[rand0:rand0+len(x)]
        list_out.append(neg0)
        list_out.append(neg1)
    return list_out

#requires every elements in the list >= 9
def keep_long(list0):
    list_out = []
    for x in list0:
        if len(x)>8:
            list_out.append(x)
    return list_out

def update_dict(dict1,dict2):
    dict0 = dict()
    ##merge two allele information by getting the better rank score from NetMHCIIpan
    for key0, value0 in dict1.iteritems():
        dict1[key0] = min(dict1[key0][2],dict2[key0][2])
    return dict0
        
def get_frag(str0,n0):
    list_out = []
    for i in range(0,len(str0)-n0+1):
        list_out.append(str0[i:i+n0]) 
    return list_out  


def make_predict_dict(gene_seq,hla1,hla2):
    dict_run_pos = dict()
    pep_list = get_frag(gene_seq,n_frag)
    dict_run_pos = predict_netmhciipan(hla1,pep_list)
    dict_second = predict_netmhciipan(hla2,pep_list)
     ##merge two allele information by getting the better rank score from NetMHCIIpan
    dict_run_pos = update_dict(dict_run_pos,dict_second)
    return dict_run_pos

def get_list_from_file(file_name0):
    file0 = open(file_name0,'r')
    list_out = []
    for x in file0:
        x = x.rstrip()
        list_out.append(x)
    return list_out



#######loading data#######################
file_pid = 'target_patient.txt'
file_name_pid = path0+ file_pid
patient_target = get_list_from_file(file_name_pid)
dict_name = 'MCL_data11_18_2015_10_9_2016v1.1.dict'
MCL_data = pickle.load(open(path0+dict_name,'r'))

#########Parameters######################
n_frag = 15
#fragment >8 to run netMHCIIpan
#dict_pos = dict()
#ict_neg = dict()
#write training data into a txt file
if len(patient_target)<1:
    patient_target = MCL_data['pid']['pid']
for pid0 in patient_target:
    if 1==1:
        hla1 = MCL_data[pid0]['HLA_typing'][-1]
        hla2 = MCL_data[pid0]['HLA_typing'][-2]
        print hla1
        ###making gene list
        gene_seq = ''
        gene_seq = gene_seq + MCL_data['Constant']['IGHM']
        gene_seq = gene_seq + MCL_data[pid0]['Variable_h_seq']
        gene_seq = gene_seq + MCL_data[pid0]['Variable_l_seq']
        seq_test = gene_seq[0:15]
        print('length of fragments to be predicted= '+str(len(gene_seq)))
        dict_out = make_predict_dict(gene_seq,hla1,hla2)
        print(pid0+', sample score:'+str(dict_out[seq_test]))
        MCL_data[pid0]['NetMCHIIpan_dict'] = dict_out
        #save
        #pickle.dump(dict_pos,open(path_save+'netmhc_predict_'+pid0+'.pos.dict','w+'))
        #pickle.dump(dict_neg,open(path_save+'netmhc_predict_'+pid0+'.neg.dict','w+')) 

pickle_dump(MCL_data,open(path0+dict_name+'with_netmhcii','w+'))
  
