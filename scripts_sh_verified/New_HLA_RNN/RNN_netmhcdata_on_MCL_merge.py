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
from scipy.stats import pearsonr

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
mhc_dict_file = 'DRB1_34_encoding.dict'
#path where the model is saved
path_model = '/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/regression_model/'
model_name0 = 'train1.tsv_s_only.txtn64_hnn128_l20.1_d0.3_test1_model.json'
weight_name0 = 'train1.tsv_s_only.txtn64_hnn128_l20.1_d0.3_test1_weight.h5'

#patients excluded
#length_max
#max0 = 74
max0 = 71

#aa encodin
dict_name = 'aa_21_sparse_encoding.dict'
dict_aa = pickle.load(open(path_encoding+dict_name,'r'))
dict_aa['_'] = np.zeros(21)
###determine the encoding size
chars = dict_aa['A']
##batch_size
b_size = 128
#mixed length size
len0_hla = 34

#ouptu file
file_name_out = path_pep+'rnn_merged_nemhc_data_on_mcl_v2_71_sparse.csv'

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
    #print(mhc_seq0)
    for x in list0:
        #print(mhc_seq0+x)
        list_out.append(mhc_seq0+x)
    return list_out

def split_x(x0,n0):
    x_fixed = x0[:,:n0,:].reshape((x0.shape[0],n0*len(chars)))
    x_variable = x0[:,n0:,:]
    return [x_fixed,x_variable]

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
    [list_pos_1_fixed, list_pos_1_var] = split_x(encoding_data(list_pos_1, max0),len0_hla)
    [list_pos_2_fixed, list_pos_2_var] = split_x(encoding_data(list_pos_2, max0),len0_hla)
    [list_neg_1_fixed, list_neg_1_var] = split_x(encoding_data(list_neg_1, max0),len0_hla)
    [list_neg_2_fixed, list_neg_2_var] = split_x(encoding_data(list_neg_2, max0),len0_hla)
    val_pos_1 = list(model0.predict_proba([list_pos_1_fixed,list_pos_1_var],batch_size=b_size).flat)
    val_pos_2 = list(model0.predict_proba([list_pos_2_fixed,list_pos_2_var],batch_size=b_size).flat)
    val_neg_1 = list(model0.predict_proba([list_neg_1_fixed,list_neg_1_var],batch_size=b_size).flat)
    val_neg_2 = list(model0.predict_proba([list_neg_2_fixed,list_neg_2_var],batch_size=b_size).flat)
    val_pos = get_max_list_from2lists(val_pos_1, val_pos_2)
    val_neg = get_max_list_from2lists(val_neg_1, val_neg_2)
    return val_pos,val_neg,val_pos_1,val_pos_2,val_neg_1,val_neg_2

def get_dict_nth_val(dict0,n0):
    list_out = []
    for _,val0 in dict0.iteritems():
        list_out.append(val0[n0])
    return list_out

def get_len(list0):
    list_out = []
    for x in list0:
        list_out.append(len(x))
    return list_out

def clean_list(list0):
    list_out = []
    for x in list0:
        if not 'o' in x and len(x)<=max0-len0_hla:
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
    done_list = ['MCL019','MCL001']
    patient_target = ['MCL001']
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
            print(pid0)
            mhc1 = MCL_data[pid0]['HLA_typing'][-1]
            mhc2 = MCL_data[pid0]['HLA_typing'][-2]
            dict_pos = pickle.load(open(path_pep+'netmhc_predict_'+pid0+'.pos.dict','r'))
            pos_reference = get_dict_nth_val(dict_pos, 0)
            dict_neg = pickle.load(open(path_pep+'netmhc_predict_'+pid0+'.neg.dict','r'))
            neg_reference = get_dict_nth_val(dict_neg, 0)
            if len(dict_pos)>1:
                list_pos = get_key_list_from_dict(dict_pos)
                list_pos = clean_list(list_pos)
                list_neg = get_key_list_from_dict(dict_neg)
                list_neg = clean_list(list_neg)
                [val_pos,val_neg,val_pos_1,val_pos_2,val_neg_1,val_neg_2] = predict_with_rnn(model_rnn,list_pos,list_neg,mhc1,mhc2)
                print('Positive predicted by RNN')
                print(val_pos[0:20])
                print('Positive predicted by NetMHCpanII')
                print(pos_reference[0:20])
                print('Pearson_regression for positive='+str(pearsonr(val_pos, pos_reference)[0]))
                print('Negative predicted by RNN')
                print(val_neg[0:20])
                print('Negative predicted by NetMHCpanII')
                print(neg_reference[0:20]) 
                print('Pearson_regression for negative='+str(pearsonr(val_neg, neg_reference)[0]))                               
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