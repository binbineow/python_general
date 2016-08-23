#make a dictionary = HLA-type to list of peptides
#code to convert a list into cluster 
from utilities import *
import math
import random
import re
import subprocess


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
version0 = '_generalv2_x_'
#v2 contains training examples with both allele 1,2 and allele 2,1
out_file_name = 'hla_ii_train_val'
note_label = 'val_note.txt'
t_ratio = 1

def make_training2(path_save,version0,pid0,set_train):        
    #set up file
    touch_file(path_save+out_file_name+version0+'.txt')
    file_out = open(path_save+out_file_name+version0+'.txt','a')
    #generate the hla type sequence
    hla_seq = dict_hla[MCL_data[pid0]['HLA_typing'][-1]] + dict_hla[MCL_data[pid0]['HLA_typing'][-2]]
    hla_seq2 = dict_hla[MCL_data[pid0]['HLA_typing'][-2]] + dict_hla[MCL_data[pid0]['HLA_typing'][-1]]
    #write in training file line by line
    for pos0 in MCL_data[pid0]['MHC2_frag']:
        set_train.add(pos0)
        file_out.write(hla_seq+pos0+'\t'+'1\n')
        file_out.write(hla_seq2+pos0+'\t'+'1\n')
        for i in range(0,t_ratio):
            rand0 = random.randint(0,len_one)
            neg0 = onegenestr[rand0:rand0+len(pos0)]
            neg0 = neg0.upper()
            if random.random() > 0.5:
                hla_seq0 = hla_seq
            else:
                hla_seq0 = hla_seq2
            file_out.write(hla_seq0+neg0+'\t'+'0\n')
            neg0 = ''.join(random.sample(pos0,len(pos0)))
            if random.random() > 0.5:
                hla_seq0 = hla_seq
            else:
                hla_seq0 = hla_seq2
            file_out.write(hla_seq0+neg0+'\t'+'0\n')
    file_out.close()

#generating validation examples using all peptides
#generating a numpy list of non-identical peptides and non-substring peptides 


#how to use it in the validation step
'''
mask_non_i = list_val == 2
# non_i includes non_sub
len_non_i = sum(list_val >= 2)
mask_non_sub = list_val == 3
len_non_sub = sum(list_val == 3)
p_predicted = model.predict_classes(X_train_p)
#recall = tp/(total positive by gold standard)
recall_non_i = sum(p_predicted[mask_non_i])
recall_non_sub = sum(p_predicted[mask_non_sub])
'''
    
# 1 = share identical string to the set_train
# 2 = non-identical but contains substring in the set_train
# 3 = non-substring
def val_judge(pos0):
    label0 = 3
    
    for x in set_train:
        if pos0 in x:
            label0 = label0 -1
            break
    if pos0 in set_train:
        label0 = label0 - 1
    
    return label0
    

def make_val2(path_save,version0,pid0):        
    #set up file for non_i and non_sub list
    touch_file(path_save+out_file_name+version0+'val_note.txt')
    note_out = open(path_save+out_file_name+version0+note_label,'a')
    # continue writing into the train and val output file
    file_out = open(path_save+out_file_name+version0+'.txt','a')
    #generate the hla type sequence
    hla_seq = dict_hla[MCL_data[pid0]['HLA_typing'][-1]] + dict_hla[MCL_data[pid0]['HLA_typing'][-2]]
    #write in training file line by line
    for pos0 in MCL_data[pid0]['MHC2_frag']:
        #output sequence and classes
        file_out.write(hla_seq+pos0+'\t'+'3\n')
        for i in range(0,t_ratio):
            rand0 = random.randint(0,len_one)
            neg0 = onegenestr[rand0:rand0+len(pos0)]
            neg0 = neg0.upper()
            file_out.write(hla_seq+neg0+'\t'+'2\n')
            neg0 = ''.join(random.sample(pos0,len(pos0)))
            file_out.write(hla_seq+neg0+'\t'+'2\n')
        #output label in order
        note_out.write(str(val_judge(pos0))+'\n')
    file_out.close()
    note_out.close()

#generating the random peptide sequence
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
#read in data
patient_val = ['MCL041','MCL128','MCL019']
MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
#initiate the training set
set_train = set()
#write training data into a txt file
pid_list = MCL_data['pid']['pid']
for pid0 in pid_list:
    if not pid0 in patient_val:
        make_training2(path_save,version0,pid0,set_train)
#write validation into the previous text file and output a list 
for pid0 in patient_val:
    make_val2(path_save, version0, pid0)
    
        
#dict_hla_pid = get_pid_for_each_hla(MCL_data,pid_list)
#print_d_list(dict_hla_pid) 
#dict_hla_pep = get_pep_for_each_hla(dict_hla_pid,MCL_data)
#print_d_list(dict_hla_pep) 
