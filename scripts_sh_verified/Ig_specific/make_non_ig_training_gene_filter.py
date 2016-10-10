##########
#this script takes in 
#MCL_all_nonIg.list as training
#MCL_all_IG_constant.list as validation
##MCL_all_V_heavy.list and MCL_all_V_light.list as held-out test set 
#(don't need to processe here yet)
#generate 1-> pos training 0 -> neg training
#generate 3-> pos valdiation 2 -> neg validation
#filtering against peptides shorter 9 (preparing for combining NetMHCpanII in the future)
######### 
from utilities import *
import math
import random
import re
import subprocess


#folder for the main peptide data
path0 = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/'
#folder for encoding dictionary
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#file for random peptide sequence
one_gene_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/human_proteinome_oneline.str'
#training and validation data save path
path_save = '/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/'
#RNASeq file if needed
#dictRNA_file = path0+'MCLRNASeq_ave.dict'
#hla_dict_file = 'DRB1_pseudo_seq.dict'
version0 = '_nonIG_v3'
#v2 contains training examples with both allele 1,2 and allele 2,1
out_file_name = 'hla_ii_train_val'
#note_label = 'val_note.txt'
t_ratio = 1
#input file names
file_nonIG= 'MCL_all_nonIg.list'
file_nonIG_gene =  'MCL_all_nonIG_gene.list'
file_constant = 'MCL_all_IG_constant.list'
filter9 = True

gene_filter = ['IGLL1', 'A2J1N5', 'HEL180', 'IGLC3', 'B1N7B6','DKFZp686C15213']

#generating the random peptide sequence
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
#read in data
list_nonig=pickle.load(open(path0+file_nonIG))
list_nonig_gene = pickle.load(open(path0+file_nonIG_gene))
list_constant=pickle.load(open(path0+file_constant))

def gene_filter_classifier(str0):
    out0 = False
    for x in gene_filter:
        if str0 == x:
            out0 = True
            break
    return out0

def make_training(path_save,version0,list_nonig):        
    #set up file
    touch_file(path_save+out_file_name+version0+'.txt')
    file_out = open(path_save+out_file_name+version0+'.txt','a')
    n_short = 0
    n_train = 0
    #generate the hla type sequence
    #hla_seq = dict_hla[MCL_data[pid0]['HLA_typing'][-1]] + dict_hla[MCL_data[pid0]['HLA_typing'][-2]]
    #hla_seq2 = dict_hla[MCL_data[pid0]['HLA_typing'][-2]] + dict_hla[MCL_data[pid0]['HLA_typing'][-1]]
    #write in training file line by line
    for n0 in range(0,len(list_nonig)):
        pos0 = list_nonig[n0]
        gene0 = list_nonig_gene[n0]
        if (not filter9 or len(pos0)>8) and (not len(gene_filter)>0 or not gene_filter_classifier(gene0)):
            file_out.write(pos0+'\t'+'1\n')
            n_train +=1
            for i in range(0,t_ratio):
                rand0 = random.randint(0,len_one)
                neg0 = onegenestr[rand0:rand0+len(pos0)]
                neg0 = neg0.upper()
                file_out.write(neg0+'\t'+'0\n')
                neg0 = ''.join(random.sample(pos0,len(pos0)))
                file_out.write(neg0+'\t'+'0\n')
        else:
            n_short +=1
    file_out.close()
    print('Total_positive_examples='+str(n_train))
    print('Total_examples_excluded_due_to_short_length='+str(n_short))
    
def make_validation(path_save,version0,list_nonig):        
    #set up file
    touch_file(path_save+out_file_name+version0+'.txt')
    file_out = open(path_save+out_file_name+version0+'.txt','a')
    n_short = 0
    n_train = 0
    #generate the hla type sequence
    #hla_seq = dict_hla[MCL_data[pid0]['HLA_typing'][-1]] + dict_hla[MCL_data[pid0]['HLA_typing'][-2]]
    #hla_seq2 = dict_hla[MCL_data[pid0]['HLA_typing'][-2]] + dict_hla[MCL_data[pid0]['HLA_typing'][-1]]
    #write in training file line by line
    for pos0 in list_nonig:
        if not filter9 or len(pos0)>8:
            file_out.write(pos0+'\t'+'3\n')
            n_train +=1
            for i in range(0,t_ratio):
                rand0 = random.randint(0,len_one)
                neg0 = onegenestr[rand0:rand0+len(pos0)]
                neg0 = neg0.upper()
                file_out.write(neg0+'\t'+'2\n')
                neg0 = ''.join(random.sample(pos0,len(pos0)))
                file_out.write(neg0+'\t'+'2\n')
        else:
            n_short +=1
    file_out.close()
    print('Total_positive_examples_for_validaiton='+str(n_train))
    print('Total_examples_excluded_due_to_short_length_for_validation='+str(n_short))


    
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



#patient_val = ['MCL041','MCL128','MCL019']
#MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
#dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
#initiate the training set
#set_train = set()
#write training data into a txt file
#pid_list = MCL_data['pid']['pid']
make_training(path_save,version0,list_nonig)
#write validation into the previous text file and output a list 
make_validation(path_save, version0, list_constant)
    
        
#dict_hla_pid = get_pid_for_each_hla(MCL_data,pid_list)
#print_d_list(dict_hla_pid) 
#dict_hla_pep = get_pep_for_each_hla(dict_hla_pid,MCL_data)
#print_d_list(dict_hla_pep) 
