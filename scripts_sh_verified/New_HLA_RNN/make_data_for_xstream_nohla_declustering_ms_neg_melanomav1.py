#10/30/2016
#make RNN training and validation file
#declustering
#no hla information
#on xstream

from utilities import *
import math
import random
import re
import subprocess

#make a list of lists/culsters, within each cluster, strings are substring among each other within a cluster
def make_cluster(list0):
    list2 = list(list0)
    len2 = []
    for x in list2:
        len2.append(len(x))
    list_len_2 = zip(list2,len2)
    list2,len2 = zip(*sorted(list_len_2,key=lambda x: x[1]))
    #print list2
    #print list2
    list2_tested = [True]*len(list2)
    cluster_list = []
#     len_d_ave = []
#     len_d_max = []
#     len_std = []
    n_element = 0
    for n0 in range(0,len(list2)):
        if list2_tested[n0]:   
            list0 = []             
            list0.append(list2[n0])
            for m0 in range(n0+1,len(list2)):
                if list2_tested[m0]:
                    if list2[n0] in list2[m0]:
                        list2_tested[m0] = False
                        list0.append(list2[m0])
            cluster_list.append(list0)
            n_element += len(list0)
    print('total cluster='+str(len(cluster_list)))
    print('total elements='+str(n_element))
    return cluster_list

def get_cluster_list(list_pep_total):
    #get clusters
    list_culsters = make_cluster(list_pep_total)
    #make the cluster random
    random.shuffle(list_culsters)
    return list_culsters

def split_cluster(list_clusters,n0,split0):
    train_target = n0*(1-split0)
    print('training_target='+str(train_target))
    list_train = []
    list_val = []
    add_to_train = True
    for cluster0 in list_clusters:
        if add_to_train:
            list_train.extend(cluster0)
        else:
            list_val.extend(cluster0)
        if len(list_train) >= train_target:
            add_to_train = False
    return list_train,list_val

def write_data_with_neg(list0,path_save,file_name0,neg_int,shuffle0):
    #touch_file is in utilities
    touch_file(path_save+out_file_name+version0+'.txt')
    file_out = open(path_save+out_file_name+version0+'.txt','a')
    for pos0 in list0:
        #output sequence and classes
        if not 'X' in pos0:
            file_out.write(pos0+'\t'+str(neg_int+1)+'\n')
            if shuffle0:
                for i in range(0,t_ratio):
                    rand0 = random.randint(0,len_one_ms)
                    neg0 = onegenestr_ms[rand0:rand0+len(pos0)]
                    neg0 = neg0.upper()
                    file_out.write(neg0+'\t'+str(neg_int)+'\n')
                    neg0 = ''.join(random.sample(pos0,len(pos0)))
                    file_out.write(neg0+'\t'+str(neg_int)+'\n')
            elif not mixed0:
                for i in range(0,t_ratio*2):
                    rand0 = random.randint(0,len_one_ms)
                    neg0 = onegenestr_ms[rand0:rand0+len(pos0)]
                    neg0 = neg0.upper()
                    file_out.write(neg0+'\t'+str(neg_int)+'\n')
            else:
                for i in range(0,t_ratio):
                    rand0 = random.randint(0,len_one_ms)
                    neg0 = onegenestr_ms[rand0:rand0+len(pos0)]
                    neg0 = neg0.upper()
                    file_out.write(neg0+'\t'+str(neg_int)+'\n')
                    rand0 = random.randint(0,len_one)
                    neg0 = onegenestr[rand0:rand0+len(pos0)]
                    neg0 = neg0.upper()
                    file_out.write(neg0+'\t'+str(neg_int)+'\n')

                
            #output label in order
    file_out.close()
    
#filter out peptide length <len0    
def limit_length(list0,len0):
    list_out = []
    n0 = 0
    for x in list0:
        if len(x) < len0:
            n0 += 1
        else:
            list_out.append(x)
    print('Peptides shorter than '+str(len0)+' ='+str(n0))
    return list_out

#filter out peptide length <len0    
def limit_long(list0,len0):
    list_out = []
    n0 = 0
    for x in list0:
        if len(x) > len0:
            n0 += 1
        else:
            list_out.append(x)
    print('Peptides shorter than '+str(len0)+' ='+str(n0))
    return list_out

#folder for the main peptide data
path0 = '/home/stanford/rbaltman/users/bchen45/data/protein_general/'
#folder for encoding dictionary
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#file for random peptide sequence
one_gene_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/human_proteinome_oneline.str'
#from ms data of jeko and 128
one_ms_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/jeko_128_ms.str'
#training and validation data save path
path_save = '/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/'
#RNASeq file if needed
#dictRNA_file = path0+'MCLRNASeq_ave.dict'
#hla_dict_file = 'DRB1_pseudo_seq.dict'
version0 = '_melanoma_deculster_ms_plusrandom_short'
#mix random peptide types or not
mixed0 = True
#v2 contains training examples with both allele 1,2 and allele 2,1
out_file_name = 'hla_i_train_val'
#note_label = 'val_note.txt'
t_ratio = 1
#validation split
split0 = 0.1 #use 20% of data as the validation
#use shuffled positive peptides as negative or not
shuffle0 = False
#MHC group
mhc0 = 'MHC1_frag'
#generating the random peptide sequence
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
onegenestr_ms = pickle.load(open(one_ms_path,'r'))
len_one_ms = len(onegenestr_ms)
#read in data
list_pep_total = pickle.load(open(path0+'nature_Bassani_2016_mhc2.list'))
#list_pep_total = list(set(list_pep_total))
print('redundant set of peptide of '+mhc0+' ='+str(len(list_pep_total)))
#make it unique
list_pep_total = list(set(list_pep_total))
#print(list_pep_total[0:15])
list_pep_total = limit_length(list_pep_total,8)
list_pep_total = limit_long(list_pep_total,13)
print('After filtering, peptide n='+str(len(list_pep_total)))
list_clusters = get_cluster_list(list_pep_total)
[list_train,list_val] = split_cluster(list_clusters,len(list_pep_total),split0)
#0 indicates the code for negative data, 0 for training, 2 for validation
del_file(path_save+out_file_name+version0+'.txt')
print('Save the file to '+path_save+out_file_name+version0+'.txt')
write_data_with_neg(list_train,path_save,out_file_name+version0,0,shuffle0)
write_data_with_neg(list_val,path_save,out_file_name+version0,2,shuffle0)
print('Training positive = '+str(len(list_train)))
print('Validation positive = '+str(len(list_val)))


