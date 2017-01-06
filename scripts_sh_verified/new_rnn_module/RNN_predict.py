from __future__ import print_function
from keras.models import Sequential
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent, Merge
from keras.callbacks import ModelCheckpoint
import cPickle as pickle
import numpy as np
#from utilities import *
from keras.models import model_from_json
from scipy.stats import percentileofscore
from keras.regularizers import l2, activity_l2
from sklearn.metrics import roc_auc_score


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
path_model = '/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/mcl_model/'
model_name0 = 'hla_ii_train_val_generalv1_x_deculster.txthla2_noshuffle_s_only_n64_hnn32_layer1_d0.3_l2_0.1_forP_model.json'
weight_name0 = 'hla_ii_train_val_generalv1_x_deculster_noshuffle.txthla2_decluster_noshuffle_Sonly_n64_hnn32_layer1_d0.3_l2_0.1_weight.h5'
#patients excluded
#length_max
max0 = 40

#aa encoding
dict_name='aa_21_sparse_encoding.dict'
dict_aa = pickle.load(open(path_encoding+dict_name,'r'))
###determine the encoding size
chars = dict_aa['A']
##batch_size
b_size = 512


##functions
def encoding_line(str0, max_len,dict0):
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
        coded0 = np.zeros((max_len,len(list(dict0['A']))))
        for i,char0 in enumerate(str0):
            coded0[i,:] = dict0[char0] 
    #print(str0)
    #print(coded0)
    return coded0

def encoding(matrix0, input0, len0,dict0):
    for i, sentence in enumerate(input0):
        matrix0[i] = encoding_line(sentence, len0,dict0)
    return matrix0

def encoding_data(list0,MAXLEN,dict0):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, len(list(dict0['A']))))
    X_encoded = encoding(X_0_m,list0,MAXLEN,dict0)
    return X_encoded


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


def predict_with_rnn(model0,list_pos_1,dict0,max_l):
    list_pos_1 = encoding_data(list_pos_1, max_l,dict0)
    val_pos_1 = model0.predict_proba(list_pos_1,batch_size=b_size,verbose=0)[:,1]
    class_pos_1 = model0.predict_classes(list_pos_1,batch_size=b_size)
    return val_pos_1,class_pos_1



def import_model(path_model, model_name0,weight_name0):
    model_name0 = path_model+ model_name0
    weight_name0 = path_model + weight_name0
    model0 = model_from_json(open(model_name0).read())
    model0.load_weights(weight_name0)
    return model0

def cal_sensi(list0,cut_off):
    n0 = 0
    for x in list0:
        if x>= cut_off:
            n0 +=1
    return float(n0)/len(list0)



def RNN_predict(list0, cut_off=0.4, path0 = path_model,model0=model_name0,weight0=weight_name0,dict0=dict_aa,max_l=max0):
    model0 = import_model(path0,model0,weight0)
    [list_scores,list_class] = predict_with_rnn(model0,list0,dict0,max_l)
    sensitivity_cut = cal_sensi(list_scores,cut_off)
    #dict_random = pickle.load(open('/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/random_pep_by_mhc.dict','r'))
    #list_random = dict_random['HLA-DRB1*01:01']
    #neg_list_scores = run_rnn_model(model0,list_random)
    sensitivity0 = np.sum(list_class)/float(len(list0))
    return list_scores,list_class, sensitivity0, sensitivity_cut
    
    
    