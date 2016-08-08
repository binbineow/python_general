#make a dictionary = HLA-type to list of peptides
#code to convert a list into cluster 
from utilities import *
#dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
path0 = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
path_encoding = '/scratch/users/bchen45/code/python_general/python_general/encoding_dict/'
one_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/human_proteinome_oneline.str'
path_save = '/scratch/users/bchen45/HLA_prediction/RNN_data/training_files/psedu_seq_train/'
hla_dict_file = 'DRB1_pseudo_seq.dict'
version0 = 'psedo_seqv1'
t_ratio = 1

import math
import random
import re
import subprocess



        
def make_training2(path_save,version0,pid0,set_train):        
    #set up file
    cmd = path_save+'hla_ii_training_'+version0+'.txt'
    cmd0 = subprocess.Popen(cmd,shell=True)      
    cmd0.wait() 
    file_out = open(path_save+'hla_ii_training'+version0+'.txt','a')
    #generate the hla type sequence
    hla_seq = dict_hla[MCL_data[pid0]['HLA_typing'][-1]] + dict_hla[MCL_data[pid0]['HLA_typing'][-2]]
    #write in training file line by line
    for pos0 in MCL_data[pid0]['MHC2_frag']:
        set_train.add(pos0)
        file_out.write(hla_seq+pos0+'\t'+'1\n')
        for i in range(0,t_ratio):
            rand0 = random.randint(0,len_one)
            neg0 = onegenestr[rand0:rand0+len(pos0)]
            file_out.write(hla_seq+neg0+'\t'+'0\n')
            neg0 = ''.join(random.sample(pos0,len(pos0)))
            file_out.write(hla_seq+neg0+'\t'+'0\n')
    file_out.close()
    
def val_judge(type0,pos0):
    b_return = False
    if type0 == 'all':
        b_return = True
    if type0 == 'identical' and not pos0 in set_train:
        b_return = True
    if type0 == 'substring':
        b_sub = False
        for x in set_train:
            if pos0 in x:
                b_sub = True
                break
        if not b_sub:
            b_return = True
    return b_return
    

def make_val2(path_save,version0,pid0):        
    #set up file
    for type0 in ['all','identical','substring']:
        cmd = path_save+'hla_ii_val'+type0+'_'+version0+'.txt'
        cmd0 = subprocess.Popen(cmd,shell=True)      
        cmd0.wait() 
        file_out = open(path_save+'hla_ii_training'+version0+'.txt','a')
        #generate the hla type sequence
        hla_seq = dict_hla[MCL_data[pid0]['HLA_typing'][-1]] + dict_hla[MCL_data[pid0]['HLA_typing'][-2]]
        #write in training file line by line
        for pos0 in MCL_data[pid0]['MHC2_frag']:
            if val_judge(type0,pos0):               
                file_out.write(hla_seq+pos0+'\t'+'3\n')
                for i in range(0,t_ratio):
                    rand0 = random.randint(0,len_one)
                    neg0 = onegenestr[rand0:rand0+len(pos0)]
                    file_out.write(hla_seq+neg0+'\t'+'2\n')
                    neg0 = ''.join(random.sample(pos0,len(pos0)))
                    file_out.write(hla_seq+neg0+'\t'+'2\n')
        file_out.close()
    
    




onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
patient_val = ['MCL041','MCL128','MCL019']
MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
set_train = set()
#MCL_data = del_sameHLA(MCL_data)
pid_list = MCL_data['pid']['pid']
for pid0 in pid_list:
    if not pid0 in patient_val:
        make_training2(path_save,version0,pid0,set_train)
for pid0 in patient_val:
    make_val2(path_save, version0, pid0)
    
        
#dict_hla_pid = get_pid_for_each_hla(MCL_data,pid_list)
#print_d_list(dict_hla_pid) 
#dict_hla_pep = get_pep_for_each_hla(dict_hla_pid,MCL_data)
#print_d_list(dict_hla_pep) 
