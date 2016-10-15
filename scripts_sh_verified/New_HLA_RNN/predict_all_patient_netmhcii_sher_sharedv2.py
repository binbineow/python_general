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

path_save = '/share/PI/rbaltman/bchen45/software/IEDB/MCL_netmhc_predict_results/'
path0 = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
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
    for key0, value0 in dict1.iteritems():
        dict1[key0] = dict1[key0] + dict2[key0]
    return dict1
        

def make_predict_dict(pid0,hla1,hla2):
    dict_run_pos = dict()
    dict_run_neg = dict()
    pep_list = keep_long(MCL_data[pid0]['MHC2_frag'])
    pep_list_neg = make_neg(pep_list)
    dict_run_pos = predict_netmhciipan(hla1,pep_list)
    dict_second = predict_netmhciipan(hla2,pep_list)
    dict_run_pos = update_dict(dict_run_pos,dict_second)
    dict_run_neg = predict_netmhciipan(hla1,pep_list_neg)
    dict_second = predict_netmhciipan(hla2,pep_list_neg)
    dict_run_neg = update_dict(dict_run_neg,dict_second)
    return [dict_run_pos,dict_run_neg]

#generating the random peptide sequence
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
#read in data
#patient_val = ['MCL041','MCL128','MCL019']
patient_target = []
done_list = ['MCL128','MCL019']
MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
#dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
#initiate the training set
#set_train = set()
#dict_pos = defaultdict(defaultdict)
#dict_neg = defaultdict(defaultdict)
dict_pos = dict()
dict_neg = dict()
#write training data into a txt file
if len(patient_target)<1:
    patient_target = MCL_data['pid']['pid']
for pid0 in patient_target:
    if not pid0 in done_list:
        hla1 = MCL_data[pid0]['HLA_typing'][-1]
        hla2 = MCL_data[pid0]['HLA_typing'][-2]
        print hla1
        [dict_pos,dict_neg] = make_predict_dict(pid0,hla1,hla2)
        print dict_pos
        #save
        pickle.dump(dict_pos,open(path_save+'netmhc_predict_'+pid0+'.pos.dict','w+'))
        pickle.dump(dict_neg,open(path_save+'netmhc_predict_'+pid0+'.neg.dict','w+')) 

  
