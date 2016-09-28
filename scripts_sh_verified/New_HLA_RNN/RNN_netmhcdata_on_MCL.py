from __future__ import print_function
from keras.models import Sequential
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent, Merge
from keras.callbacks import ModelCheckpoint
from utilities import *
from keras.models import model_from_json
from scipy.stats import percentileofscore
from keras.regularizers import l2, activity_l2
from sklearn.metrics import roc_auc_score
from collections import defaultdict
from random import shuffle
import random

#parameters
#folder where pid -> strings dicts are saved
path_pep = '/home/stanford/rbaltman/users/bchen45/results/MCL_netmhc_predict_results/'
#example file netmhc_predict_MCL034.neg.dict
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
mhc_dict_file = 'DRB1_pseudo_seq.dict'
#path where the model is saved
path_model = '/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/'
model_name0 = 'netMHCIIpan_train1.tab_chems.txtn64_final_hnn0_l20.1_d0.2_reluv2_model.json'
weight_name0 = 'netMHCIIpan_train1.tab_chems.txtn64_final_hnn0_l20.1_d0.2_reluv2_weight.h5'
#patients excluded
#length_max
max0 = 74
#aa encoding
dict_name='Blosum50_sparse.dict'
dict_aa = pickle.load(open(path_encoding+dict_name,'r'))
###determine the encoding size
chars = dict_aa['A']
##batch_size
b_size = 128

#ouptu file
file_name_out = path_pep+'rnn_nemhc_data_on_mcl_v2_sub.txt'

#mhc_pseudosequence_dict
mhc_dic = pickle.load(open(path_encoding+mhc_dict_file,'r'))

##functions
def encoding_line(str0, max_len):
    #print(type(dict_aa['A']))
    #print(type(list(dict_aa['A'])))
    #print(type(max_len))
    if len(str0) == 1:
        coded0 = np.zeros(2)
        if str0 == '0':
            coded0[0] = 1
        else:
            coded0[1] = 1
    else:
        coded0 = np.zeros((max_len,len(list(dict_aa['A']))))
        for i,char0 in enumerate(str0):
            coded0[i,:] = dict_aa[char0] 
    #print(str0)
    #print(coded0)
    return coded0

def encoding(matrix0, input0, len0):
    for i, sentence in enumerate(input0):
        matrix0[i] = encoding_line(sentence, len0)
    return matrix0

def encoding_data(list0,MAXLEN):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, len(chars)))
    X_encoded = encoding(X_0_m,list0,MAXLEN)
    return X_encoded

def import_model(path_model,model_name0,weight_name0):
    model_name0 = path_model+ model_name0
    weight_name0 = path_model + weight_name0
    model0 = model_from_json(open(model_name0).read())
    model0.load_weights(weight_name0)
    return model0

##list1,2 are two possible binding probabilities give two potential MHC types
def cal_percentile_from2lists(list1,list2,distr1,distr2):
    distr1 = list(distr1.flat)
    distr2 = list(distr2.flat)
    list_out = []
    for n0 in range(0,len(list1)):
        #val1 = percentileofscore(np.array(distr1), list1[n0])
        #val1 = percentileofscore(distr1, list1[n0])
        #print(val1)
        #print(list1[n0])
        max0 = max(percentileofscore(distr1, list1[n0]),percentileofscore(distr2, list2[n0]))
        list_out.append(max0)
    return list_out

def get_max_list_from2lists(list1,list2):
    list_out = []
    for n0 in range(0,len(list1)):
        list_out.append(max(list1[n0],list2[n0]))
    return list_out
        

def cal_auc_from2lists(post_list,neg_list):
    list_values = np.concatenate((post_list,neg_list))
    list_true = np.concatenate((np.ones(len(post_list)),np.zeros(len(neg_list))))
    auc_val = roc_auc_score(list_true, list_values)
    return auc_val

def add_mhc_to_peplist(mhc0,list0):
    list_out = []
    mhc_seq0 = mhc_dic[mhc0]
    for x in list0:
        list_out.append(mhc_seq0+x)
    return list_out

def make_random_dict(model0,mhc_set,len_list,file_name0=path_save+'random_pep_by_mhc.dict'):
    dict_random_mhc = dict()
    onegenestr = pickle.load(open(one_gene_path,'r'))
    len_one = len(onegenestr)
    list_random = []
    for len0 in len_list:
        rand0 = random.randint(0,len_one)
        neg0 = onegenestr[rand0:rand0+len0]
        list_random.append(neg0)
    for mhc0 in mhc_set:
        list_with_seq = encoding_data(add_mhc_to_peplist(mhc0, list_random),max0)
        list_val_p = model0.predict_proba(list_with_seq,batch_size=b_size)
        dict_random_mhc[mhc0] = list_val_p
    pickle.dump(dict_random_mhc,open(file_name0,'w+'))
    print('dcit_random_mhc is saved at '+file_name0)
    return dict_random_mhc

def get_key_list_from_dict(dict0):
    list_out = []
    for key0,_ in dict0.iteritems():
        list_out.append(key0)
    return list_out

def predict_with_rnn(model0,list_pos,list_neg,mhc1,mhc2):
    list_pos_1 = add_mhc_to_peplist(mhc1,list_pos)
    list_pos_2 = add_mhc_to_peplist(mhc2,list_pos)
    list_neg_1 = add_mhc_to_peplist(mhc1,list_neg)
    list_neg_2 = add_mhc_to_peplist(mhc2,list_neg)
    list_pos_1 = encoding_data(list_pos_1, max0)
    list_pos_2 = encoding_data(list_pos_2, max0)
    list_neg_1 = encoding_data(list_neg_1, max0)
    list_neg_2 = encoding_data(list_neg_2, max0)
    val_pos_1 = list(model0.predict_proba(list_pos_1,batch_size=b_size).flat)
    val_pos_2 = list(model0.predict_proba(list_pos_2,batch_size=b_size).flat)
    val_neg_1 = list(model0.predict_proba(list_neg_1,batch_size=b_size).flat)
    val_neg_2 = list(model0.predict_proba(list_neg_2,batch_size=b_size).flat)
    val_pos = get_max_list_from2lists(val_pos_1, val_pos_2)
    val_neg = get_max_list_from2lists(val_neg_1, val_neg_2)
    return val_pos,val_neg,val_pos_1,val_pos_2,val_neg_1,val_neg_2

def get_len(list0):
    list_out = []
    for x in list0:
        list_out.append(len(x))
    return list_out

def clean_list(list0):
    list_out = []
    for x in list0:
        if not 'o' in x:
            list_out.append(x)
    return list_out

def process_data(model_rnn,file_name0):
    MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
    #dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
    #initiate the training set
    #set_train = set()
    #dict_pos = defaultdict(defaultdict)
    #dict_neg = defaultdict(defaultdict)
    #write training data into a txt file
    patient_target = []
    done_list = ['MCL019']
    patient_target = []
    output_string = ''
    if len(patient_target)<1:
        patient_target = MCL_data['pid']['pid']
    #get mhc_set   
    #create 10000 random peptide sequence lenth 
    '''
    mhc_set = set()
    len_list = []
    for pid0 in MCL_data['pid']['pid']:
        mhc_set.add(MCL_data[pid0]['HLA_typing'][-1])
        mhc_set.add(MCL_data[pid0]['HLA_typing'][-2])
        len_list.extend(get_len(MCL_data[pid0]['MHC2_frag']))
    #shuffle len_list
    shuffle(len_list)
    #get to first 10,000
    n0 = 10000
    len_list = len_list[0:n0]
    
    #create random peptide dictioanry based on mhc type
    dict_random = make_random_dict(model_rnn,mhc_set,len_list)
    '''
    dict_random = pickle.load(open('/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/random_pep_by_mhc.dict','r'))
    print('start prediction')
    for pid0 in patient_target:
        if not pid0 in done_list:
            #print(pid0)
            mhc1 = MCL_data[pid0]['HLA_typing'][-1]
            mhc2 = MCL_data[pid0]['HLA_typing'][-2]
            dict_pos = pickle.load(open(path_pep+'netmhc_predict_'+pid0+'.pos.dict','r'))
            dict_neg = pickle.load(open(path_pep+'netmhc_predict_'+pid0+'.neg.dict','r'))
            if len(dict_pos)>1:
                list_pos = get_key_list_from_dict(dict_pos)
                list_pos = clean_list(list_pos)
                list_neg = get_key_list_from_dict(dict_neg)
                list_neg = clean_list(list_neg)
                [val_pos,val_neg,val_pos_1,val_pos_2,val_neg_1,val_neg_2] = predict_with_rnn(model_rnn,list_pos,list_neg,mhc1,mhc2)
                auc_raw = cal_auc_from2lists(val_pos,val_neg)
                val_pos_per = cal_percentile_from2lists(val_pos_1,val_pos_2,dict_random[mhc1],dict_random[mhc2])
                val_neg_per = cal_percentile_from2lists(val_neg_1,val_neg_2,dict_random[mhc1],dict_random[mhc2])
                auc_percent = cal_auc_from2lists(val_pos_per,val_neg_per)
                str0 = pid0+','+mhc1+','+mhc2+','+str(auc_raw)+','+str(auc_percent)+'\n'
                print(str0)
                output_string = output_string + str0
    return output_string

def main(file_name_out):
    model_rnn = import_model(path_model, model_name0, weight_name0)
    output_str = process_data(model_rnn,file_name_out)
    file_out = open(file_name_out,'w+')
    file_out.write(output_str)
    file_out.close()

main(file_name_out)