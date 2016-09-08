#make a dictionary = HLA-type to list of peptides
#code to convert a list into cluster 
from utilities import *
from predict_netmhciipan import *
import math
import random
import re
import subprocess
from collections import defaultdict


#folder for the main peptide data
path0 = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/'
#folder for encoding dictionary
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#file for random peptide sequence
one_gene_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/human_proteinome_oneline.str'
#training and validation data save path
path_save = '/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/'
#RNASeq file if needed
#dictRNA_file = path0+'MCLRNASeq_ave.dict'
hla_dict_file = 'DRB1_pseudo_seq.dict'
version0 = '_generalv1_x_'
#v2 contains training examples with both allele 1,2 and allele 2,1
out_file_name = 'hla_ii_train_val'
note_label = 'val_note.txt'
t_ratio = 1



def make_neg(pep_list):
    list_out = []
    for x in pep_list:
        neg0 = ''.join(random.sample(x,len(x)))
        rand0 = random.randint(0,len_one)
        neg1= onegenestr[rand0:rand0+len(x)]
        list_out.append(neg0)
        list_out.append(neg1)

#requires every elements in the list >= 9
def keep_long(list0):
    list_out = []
    for x in list0:
        if len(x)>8:
            list_out.append(x)
    return list_out

def update_dict(dict1,dict2):
    for key0, value0 in dict1.iteritems():
        dict1[key0] = dict[key0] + dict[key0]
    return dict1
        

def make_predict_dict(pid0,hla1,hla2):
    dict_run_pos = dict()
    dict_run_neg = dict()
    pep_list = keep_long(MCL_data[pid0]['MHC2_frag'])
    pep_list_neg = make_neg(pep_list)
    dict_run_pos = predict_netmhciipan(hla1,pep_list)
    dict_second = predict_netmhciipan(hla2,peplist)
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

MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
#initiate the training set
set_train = set()
dict_pos = defaultdict(defaultdict)
dict_neg = defaultdict(defaultdict)
#write training data into a txt file
pid_list = MCL_data['pid']['pid']
for pid0 in pid_list:
    hla1 = MCL_data[pid0]['HLA_typing'][-1]
    hla2 = MCL_data[pid0]['HLA_typing'][-2]
    [dict_pos[pid0],dict_neg[pid0]] = make_predict_dict(pid0,hla1,hla2)

#save
pickle.dump(dict_pos,open(path_save+'netmhc_predict_mcl.pos.dict'))
pickle.dump(dict_neg,open(path_save+'netmhc_predict_mcl.neg.dict'))   
