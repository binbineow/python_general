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
pathig = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/'
file_name0 = 'ig_constant_region.txt'
mhc_file0 = 'MCL_all_IG_constant.list'
#folder for encoding dictionary
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#file for random peptide sequence
one_gene_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/human_proteinome_oneline.str'
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
#training and validation data save path
path_save = '/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/'
#RNASeq file if needed
#dictRNA_file = path0+'MCLRNASeq_ave.dict'
mhc_dict_file = 'DRB1_pseudo_seq.dict'
#path where the model is saved
path_model = '/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/'
model_name0 = 'hla_ii_train_val_nonIG_v3.txtnon_ig_v3_Sonly_n64_layer1_d0.3_l2_0.1_hnn32_model.json'
weight_name0 = 'hla_ii_train_val_nonIG_v3.txtnon_ig_v3_Sonly_n64_layer1_d0.3_l2_0.01_weight.h5'
#sb model
#model_name0 = 'hla_ii_train_val_nonIG_v1.txtnon_ig_v1_BSonly_n64_layer1_d0.3_l2_0.1_hnn32_model.json'
#weight_name0 = 'hla_ii_train_val_nonIG_v1.txtnon_ig_v1_BSonly_n64_layer1_d0.3_l2_0.1_hnn32_weight.h5'
#sonly
#model_name0 = 'hla_ii_train_val_nonIG_v2.txt_model.json'
#weight_name0 = 'hla_ii_train_val_nonIG_v2.txtnon_ig_v2_Sonly_n64_layer1_d0.3_l2_0.01_weight.h5'

#length of framments
n_frag = 15

#patients excluded
#length_max
max0 = 74
#max0 = 56

#verbose
vb0=1
#aa encoding
dict_name_for_output = 'SparseOnly'
dict_name='aa_21_sparse_encoding.dict'
#dict_name_for_output = 'Sparse+BLOSUM50'
#dict_name='Blosum50_sparse.dict'
dict_aa = pickle.load(open(path_encoding+dict_name,'r'))
###determine the encoding size
chars = dict_aa['A']
##batch_size
b_size = 128

#ouptu file
file_name_out = path_model+'perf'+dict_name_for_output+'_variable_region_test.csv'

#mhc_pseudosequence_dict
#mhc_dic = pickle.load(open(path_encoding+mhc_dict_file,'r'))

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
    #print(mhc_seq0)
    for x in list0:
        #print(mhc_seq0+x)
        list_out.append(mhc_seq0+x)
    return list_out


def get_key_list_from_dict(dict0):
    list_out = []
    for key0,_ in dict0.iteritems():
        list_out.append(key0)
    return list_out

def predict_with_rnn(model0,list_pos):

    list_pos_1 = encoding_data(list_pos, max0)

    #list_val_p = model.predict_proba(X_val_p,verbose=vb0)[:,1]
    val_pos_1 = list(model0.predict_proba(list_pos_1,batch_size=b_size,verbose=vb0)[:,1])


    return val_pos_1

def get_len(list0):
    list_out = []
    for x in list0:
        list_out.append(len(x))
    return list_out

def clean_list(list0):
    list_out = []
    for x in list0:
        if not 'o' in x and len(x)<=max0-19:
            list_out.append(x)
    return list_out

#helping Function
def get_list_from_file(file_name0):
    file0 = open(file_name0,'r')
    list_out = []
    for x in file0:
        x = x.rstrip()
        list_out.append(x)
    return list_out


def get_v_pos_list(list0):
    list_out = []
    for file0 in list0:
        list_out.extend(pickle.load(open(file0)))
    return list_out

def get_neg_from_pos(pos_data):
    list_neg = []
    for pos0 in pos_data:
        rand0 = random.randint(0,len_one)
        neg0 = onegenestr[rand0:rand0+len(pos0)]
        neg0 = neg0.upper()
        list_neg.append(neg0)  
        neg0 = ''.join(random.sample(pos0,len(pos0)))
        list_neg.append(neg0) 
    return list_neg        

def get_line_from_list(file_name0,str0):
    read0 = False
    for line0 in open(file_name0,'rU'):
        line0 = line0.rstrip()
        if read0:
            out_line = line0
            break
        if str0 in line0:
            read0 = True
    return out_line

def get_frag(str0,n0):
    list_out = []
    for i in range(0,len(str0)-n0+1):
        list_out.append(str0[i:i+n0]) 
    return list_out  

def mark_list(list0,index0,len0):
    for i in range(index0,index0+len0):
        list0[i] += 1     

def get_reco_map(mhc2_pos_data,seq_data):
    list_out = [0]*len(seq_data)
    for x in mhc2_pos_data:
        index0 = seq_data.find(x)
        if index0 > 0:
            mark_list(list_out,index0,len(x))
    return list_out
            
def print_list0(list0,del0):
    str0 = ''
    for x in list0:
        str0=str(x)+del0
    return str0

def process_data(file_name0,mhc_file0,model_rnn):
    seq_data = get_line_from_list(file_name0,'IGHM')
    ighm_pred = [0]*len(seq_data)
    seq_frag = get_frag(seq_data,n_frag)
    #pos_data = list(set(pos_data))
    mhc2_pos_data = pickle.load(open(mhc_file0))
    ighm_reco = get_reco_map(mhc2_pos_data,seq_data)
    print('Positive Ig Variable peptides recovered from MHCII='+str(len(mhc_file0)))
    print('start prediction')
    val_pos = predict_with_rnn(model_rnn,seq_frag)
    for i in range(0,len(val_pos)):
        ighm_pred = val_pos[i]
    print('Model_used='+model_name0)
    print('Weight_used='+weight_name0)
    print_list0(list(ighm_pred),',')
    print_list0(list(ighm_reco),',')
    return list(ighm_pred), list(ighm_reco)


def main(file_name0,mhc_file0,file_out):
    model_rnn = import_model(path_model, model_name0, weight_name0)
    [ighm_pred,ighm_reco] = process_data(file_name0,mhc_file0, model_rnn)
    pickle.dump(ighm_pred,open(path_ig+file_out+'predicted.list'))
    pickle.dump(ighm_reco,open(path_ig+file_out+'recovered.list'))

file_names_v = [pathig+'MCL_all_V_heavy.list',pathig+'MCL_all_V_light.list']
file_out = 'constantv1.txt'
main(pathig+file_name0,pathig+mhc_file0,file_out)
      