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
from sklearn.metrics import roc_curve, auc
from collections import defaultdict
from random import shuffle
import random
import pandas as pd


#parameters
#folder where pid -> strings dicts are saved
path_pep = '/home/stanford/rbaltman/users/bchen45/results/MCL_netmhc_predict_results/'
#example file netmhc_predict_MCL034.neg.dict
#folder for the main peptide data
path0 = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/'
pathig = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/'
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

#cut_off for positive peptide calling
cut_off = 0.34

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

def Find_Optimal_Cutoff(target, predicted):
    """ Find the optimal probability cutoff point for a classification model related to event rate
    Parameters
    ----------
    target : Matrix with dependent or target data, where rows are observations

    predicted : Matrix with predicted data, where rows are observations

    Returns
    -------     
    list type, with optimal cutoff value

    """
    fpr, tpr, threshold = roc_curve(target, predicted)
    i = np.arange(len(tpr)) 
    roc = pd.DataFrame({'tf' : pd.Series(tpr-(1-fpr), index=i), 'threshold' : pd.Series(threshold, index=i)})
    roc_t = roc.ix[(roc.tf-0).abs().argsort()[:1]]

    return list(roc_t['threshold']) 


def cal_pos_acc(list0,cut_off):
    n_pos = 0
    for x in list0:
        if x>= cut_off:
            n_pos += 1
    print('Recall='+str(n_pos/float(len(list0))))
    
def cal_neg_acc(name0,list0,cut_off):
    n_neg = 0
    for x in list0:
        if x< cut_off:
            n_neg += 1
    print('Specificity for ' +name0 +' = '+str(n_neg/float(len(list0))))

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
        

def cal_auc_from2lists(post_list,neg_list,neg_list_shuffle):
    list_values = np.concatenate((post_list,neg_list))
    list_true = np.concatenate((np.ones(len(post_list)),np.zeros(len(neg_list))))
    auc_val = roc_auc_score(list_true, list_values)
    optimal0 = Find_Optimal_Cutoff(list_true, list_values)
    print('Optimal threshold='+str(optimal0))
    cal_pos_acc(post_list, cut_off)
    cal_neg_acc('Random peptides', neg_list, cut_off)
    cal_neg_acc('Shuffle peptides',neg_list_shuffle,cut_off)
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

def predict_with_rnn(model0,list_pos,list_neg,list_neg2):

    list_pos_1 = encoding_data(list_pos, max0)
    list_neg_1 = encoding_data(list_neg, max0)
    list_neg_2 = encoding_data(list_neg2, max0)
    #list_val_p = model.predict_proba(X_val_p,verbose=vb0)[:,1]
    val_pos_1 = list(model0.predict_proba(list_pos_1,batch_size=b_size,verbose=vb0)[:,1])
    val_neg_1 = list(model0.predict_proba(list_neg_1,batch_size=b_size,verbose=vb0)[:,1])
    val_neg_2 = list(model0.predict_proba(list_neg_2,batch_size=b_size,verbose=vb0)[:,1])
    return val_pos_1,val_neg_1,val_neg_2

def get_val_with_rnn(model0,list0):
    list0 = encoding_data(list0, max0)
    val_list0 = list(model0.predict_proba(list0,batch_size=b_size,verbose=vb0)[:,1])
    return val_list0

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
    list_neg_pro = []
    list_neg_shuffle = []
    for pos0 in pos_data:
        rand0 = random.randint(0,len_one)
        neg0 = onegenestr[rand0:rand0+len(pos0)]
        neg0 = neg0.upper()
        list_neg_pro.append(neg0)  
        neg0 = ''.join(random.sample(pos0,len(pos0)))
        list_neg_shuffle.append(neg0) 
    return list_neg_pro,list_neg_shuffle        

def print_top(name0,list0,p_list):
    for p0 in p_list:
        val0 = str(np.percentile(list0, p0))
        print(name0+',Top'+str(p0)+','+str(val0))
    

def process_data(list0,model_rnn):
    pos_data = get_v_pos_list(list0)
    pos_data = list(set(pos_data))
    print('Positive Ig Constant peptides recovered from MHCII='+str(len(pos_data)))
    print(pos_data[0:10])
    [list_neg_pro,list_neg_shuffle] = get_neg_from_pos(pos_data)
    #dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
    #initiate the training set
    #set_train = set()
    #dict_pos = defaultdict(defaultdict)
    #dict_neg = defaultdict(defaultdict)
    #write training data into a txt file
    #patient_target = []
    #done_list = ['MCL019']
    #patient_target = ['MCL001']
    #patient_target = []
    #output_string = ''
    #if len(patient_target)<1:
    #    patient_target = MCL_data['pid']['pid']
    #get mhc_set   
    #create 10000 random peptide sequence lenth 
    #dict_random = pickle.load(open('/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/random_pep_by_mhc.dict','r'))
    print('start prediction')
    [val_pos,val_neg_pro,val_neg_shuffle] = predict_with_rnn(model_rnn,pos_data,list_neg_pro,list_neg_shuffle)
    auc_raw = cal_auc_from2lists(val_pos,val_neg_pro,val_neg_shuffle)
    print('Model_used='+model_name0)
    print('Weight_used='+weight_name0)
    print('AUC='+str(auc_raw))
    test_data = pos_data*100
    [list_neg_pro,list_neg_shuffle] = get_neg_from_pos(test_data)
    val_pro = get_val_with_rnn(model_rnn, list_neg_pro)
    val_shuffle = get_val_with_rnn(model_rnn, list_neg_shuffle)
    print('Random peptides='+str(len(val_shuffle)))
    print_top('Random peptides',val_pro,[90,95,99])
    print_top('Shuffled peptides',val_shuffle,[90,95,99])



def main(list0):
    model_rnn = import_model(path_model, model_name0, weight_name0)
    process_data(list0, model_rnn)

#file_names_v = [pathig+'MCL_all_V_heavy.list',pathig+'MCL_all_V_light.list']
file_names_v = [pathig+'MCL_all_IG_constant.list']
main(file_names_v)
      