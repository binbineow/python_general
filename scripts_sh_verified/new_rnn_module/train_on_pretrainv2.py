# -*- coding: utf-8 -*-
#
###with regulizer and drop out
# how to use the model elsewhere...
#model = model_from_json(open('my_model_architecture.json').read())
#model.load_weights('my_model_weights.h5')

#####################import#################################
from __future__ import print_function
import numpy as np
import cPickle as pickle
import fileinput 
from keras.models import Sequential, Model
from keras.layers.core import Activation, Masking, Dropout, Dense
from keras.layers import recurrent, Input, Merge
from keras.layers import Bidirectional
from keras.layers.merge import concatenate
from keras.callbacks import ModelCheckpoint
#import keras.kernel_constraint
from keras.models import model_from_json
from scipy.stats import pearsonr
from keras.constraints import maxnorm
from keras.regularizers import l2
import keras
print(keras.__version__)
# %matplotlib inline
# import matplotlib
# import matplotlib.pyplot as plt

#encoding will take a string or char, string=sequence and to return a matrix of encoded peptide sequence
#char = class, '0' = non-binding (0,1), '1' = binding (1,0)
def encoding_line(str0, max_len, char_len, dict_aa, classn = 2):

    coded0 = np.zeros((max_len,char_len))
    for i,char0 in enumerate(str0):
        coded0[i,:] = dict_aa[char0] 
    #print(str0)
    #print(coded0)
    return coded0

def encoding(matrix0, input0, len0,dict_aa, char_len):
    for i, sentence in enumerate(input0):
        matrix0[i] = encoding_line(sentence, len0, char_len,dict_aa)
    return matrix0

def output_perf2(list0):
    touch_file(path_save+performance_file_name+'.txt')     
    file_out = open(path_save+performance_file_name+'.txt','a')
    for x in list0:
        file_out.write(str(x)+'\t')
    file_out.write('\n')
    file_out.close()

def shuffle_train(list1,list2):
    from random import shuffle
    # Given list1 and list2
    list1_shuf = []
    list2_shuf = []
    index_shuf = range(len(list1))
    shuffle(index_shuf)
    for i in index_shuf:
        list1_shuf.append(list1[i])
        list2_shuf.append(list2[i])
    return [list1_shuf,list2_shuf]

#you probably don't need this
#AUC evaluation is potentially dsired
def calf1(str1,str2):
    pre0 = float(str1)
    recall0 = float(str2)
    f1_out = 2.0*pre0*recall0/(pre0+recall0)
    return str(f1_out)

# def read_data(path_file0):
#     #read
#     X_0 = []
#     y_0 = []
#     max_len0 = 0
#     for line0 in open(path_file0,'r'):
#         line0 = line0.rstrip().split('\t')
#         X_0.append(line0[0]+line0[1])
#         y_0.append(float(line0[2]))
#         max_len0 = max(max_len0,len(line0[0]+line0[1]))
#     return [X_0,y_0,max_len0]

def encoding_data(list0,MAXLEN,dict_aa,char_len ):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0, MAXLEN, dict_aa, char_len)
    return X_encoded

def encoding_fixed(list0,MAXLEN,dict_aa,char_len,add_placeh=0):
    X_0_m = np.zeros((len(list0), MAXLEN, char_len))
    X_encoded = encoding(X_0_m,list0,MAXLEN, dict_aa,char_len)
    X_encoded = X_encoded.reshape((X_encoded.shape[0],MAXLEN*char_len))
    if add_placeh > 0:
        add_block = np.ones((X_encoded.shape[0],add_placeh))
        #print(X_encoded.shape)
        X_encoded = np.concatenate((X_encoded,add_block),axis=1)
        print(X_encoded.shape)
    return np.array(X_encoded)
    
    
    
def encoding_y(list0,class0=2):
    list_out = []
    for x in list0:
        vec0 = [0]*class0
        vec0[int(x)] = 1
        list_out.append(vec0)
    return list_out


#positive int0 -> 1
#negative int0 -> 0
def make_and_encode_label(len0,int0,class0=2):
    list_temp = [int0]*len0
    list_out = encoding_y(list_temp,class0=class0)
    return list_out

def encode_apair(data_pair,dict_aa,len_mhc,MAXLEN,len_char):
    x_mhc = encoding_fixed(data_pair[0],len_mhc,dict_aa,len_char,add_placeh=3)
    x_seq = encoding_data(data_pair[1],MAXLEN,dict_aa,len_char)
    return [x_mhc, x_seq]


###########get data from my MCL patient package
#get data
path_store = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/'
file_store = 'mcl_patient_mhc2.training_list1_v2'
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
dict_name = 'aa_21_sparse_encoding.dict'
MAXLEN = 26
dict_aa = pickle.load(open(path_encoding+dict_name,'r'))
[train_pos0,train_pos1,val_pos0,val_pos1,train_neg,val_neg, 
  list_random,list_ms_random,list_shuffle_random,
  list_anti_patient_train,list_anti_patient_val, iedb_train ,
  iedb_train_y ,iedb_val ,iedb_val_y] = pickle.load(open(path_store+file_store,'r'))

#get mhc length and dictionary length
len_mhc = len(train_pos0[0][0])
len_char = len(dict_aa['A'])

#encoding training data
y_train_neg = make_and_encode_label(len(train_neg[0]),0)
y_train_pos = make_and_encode_label(len(train_pos0[0]),1)
y_train_iedb = list(encoding_y(iedb_train_y))
y_train = y_train_pos + y_train_neg + y_train_iedb
x_train_pos0 = encode_apair(train_pos0,dict_aa,len_mhc,MAXLEN,len_char)
x_train_pos1 = encode_apair(train_pos1,dict_aa,len_mhc,MAXLEN,len_char)
x_train_neg = encode_apair(train_neg,dict_aa,len_mhc,MAXLEN,len_char)
x_train_iedb = encode_apair(iedb_train,dict_aa,len_mhc,MAXLEN,len_char)

#encoding validation data
y_val_neg = make_and_encode_label(len(val_neg),0)
y_val_pos = make_and_encode_label(len(val_pos0),1)
#y_val = y_val_pos + y_val_neg
x_val_pos0 = encode_apair(val_pos0,dict_aa,len_mhc,MAXLEN,len_char)
x_val_pos1 = encode_apair(val_pos1,dict_aa,len_mhc,MAXLEN,len_char)
x_val_neg = encode_apair(val_neg,dict_aa,len_mhc,MAXLEN,len_char)

#x_val_iedb = encode_apair(iedb_val,dict_aa,len_mhc,MAXLEN,len_char


#encoding other lists
x_list_random = encode_apair(list_random,dict_aa,len_mhc,MAXLEN,len_char)
x_list_ms_random = encode_apair(list_ms_random,dict_aa,len_mhc,MAXLEN,len_char)
x_list_shuffle_random = encode_apair(list_shuffle_random,dict_aa,len_mhc,MAXLEN,len_char)
x_list_anti_patient_train = encode_apair(list_anti_patient_train,dict_aa,len_mhc,MAXLEN,len_char)
x_list_anti_patient_val = encode_apair(list_anti_patient_val,dict_aa,len_mhc,MAXLEN,len_char)

############get model


path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'

def import_model(path_model, model_name0,weight_name0):
    model_name0 = path_model+ model_name0
    weight_name0 = path_model + weight_name0
    model0 = model_from_json(open(model_name0).read())
    #model0.load_weights(weight_name0)
    return model0
#this model has training AUC 0.95 and validation AUC 0.85
#model1 = 'mhc2_iedb_binding_training.list_iedb_pretrain_v1n64_f64_h64_d0.43_l20_layer2_sparse_masking_model.json'
#weight1 = 'mhc2_iedb_binding_training.list_iedb_pretrain_v1n64_f64_h64_d0.43_l20_layer2_sparse_maskinglstm2_weight.h5'

#acutally all neuron numbers = 64 the same neuron connections with dropout = .4 and AUC 0.837 on validation
model1 = 'mhc2_iedb_binding_training.list_iedb_pretrain_v1n128_f128_h128_d0.3_l20.1_layer2_sparse_masking_v1_model.json'
model1 = 'rnn_combine_train_modelv2.2_d0.3.json'
#weight1 = 'mhc2_iedb_binding_training.list_iedb_pretrain_v1n128_f128_h128_d0.3_l20.1_layer2_sparse_masking_v1lstm_0.837_weight.h5'
weight1 = 'rnn_combine_train_modelv1.weight2.2'

model_merge = import_model(path_save,model1,weight1)
rmsprop = keras.optimizers.RMSprop(lr=0.001)
model_merge.compile(loss='categorical_crossentropy', optimizer=rmsprop)
#model_merge.summary()

#supporting funcitons
#for merge model
def cal_performance(model0,x,y,batch0=1024):
    from sklearn.metrics import roc_auc_score
    list_predicted = model0.predict(x,verbose=1,batch_size=batch0)[:,1]
    auc_val = roc_auc_score(y, list_predicted)
    print('AUC='+str(auc_val))
    return auc_val

def cal_auc_pos_neg(pre_pos,pre_neg):
    from sklearn.metrics import roc_auc_score
    list_scores = list(pre_pos) + list(pre_neg)
    list_true = [1]*len(pre_pos) + [0]*len(pre_neg)
    auc_val = roc_auc_score(list_true, list_scores)
    print('AUC='+str(auc_val))
    return auc_val

def get_percent_pos(model0,x,cut_off=0.5,batch0=1024):
    list_scores = model0.predict(x,verbose=1,batch_size=batch0)[:,1]
    percent0 = np.sum(list_scores >= 0.5)/float(len(x[0]))
    print(percent0)
    return percent0
    
##checking the model on training data
def check_model(model_merge, x_train_pos0,x_train_pos1,x_train_neg,batch0 = 1024,record_file=1):
    list_predicted_pos0 = model_merge.predict(x_train_pos0,verbose=1,batch_size=batch0)[:,1]
    list_predicted_pos1 = model_merge.predict(x_train_pos1,verbose=1,batch_size=batch0)[:,1]
    list_predicted_neg = model_merge.predict(x_train_neg,verbose=1,batch_size=batch0)[:,1]
    #np.sum(list_predicted_neg < 0.5)/float(len(list_predicted_neg))
    list_predict_pos = [max(x1,x2) for x1,x2 in zip(list_predicted_pos0,list_predicted_pos1)]
    #list_predict_scores = list_predict_pos + list_predicted_neg
    #list_true = [1]*len(list_predict_pos) + [0] * list_predicted_neg
    #print(np.sum(list_predict_pos >= 0.5)/float(len(list_predict_pos)))
    print('Sensitivity:')
    sens0 = np.sum(np.array(list_predict_pos) >= 0.5)/float(len(list_predict_pos))
    print(sens0)
    print('Specificity:')
    spec0 = 1-np.sum(np.array(list_predicted_neg) >= 0.5)/float(len(list_predicted_neg))
    print(spec0)
    if record_file != 1:
        record_file.write('Sensitivity='+str(sens0)+' Specificity='+str(spec0)+' ')
    #plt.scatter(list_predicted_pos0,list_predicted_pos1)
    auc_model = cal_auc_pos_neg(list_predict_pos,list_predicted_neg)
    return auc_model
    
#check_model(model_merge, x_train_pos0,x_train_pos1,x_train_neg)

#####functions supporting fitting regarding two alleles
def get_pos_for_fit(model_merge,x_train_pos0,x_train_pos1,batch0= 1024):
    list_predicted_pos0 = model_merge.predict(x_train_pos0,verbose=1,batch_size=batch0)[:,1]
    list_predicted_pos1 = model_merge.predict(x_train_pos1,verbose=1,batch_size=batch0)[:,1]
    list_out = [[],[]]
    for i in range(0,len(list_predicted_pos0)):
        if list_predicted_pos0[i] > list_predicted_pos1[i]:
            list_out[0].append(list(x_train_pos0[0][i]))
            list_out[1].append(list(x_train_pos0[1][i]))
        else:
            list_out[0].append(list(x_train_pos1[0][i]))
            list_out[1].append(list(x_train_pos1[1][i]))
    return list_out

#############training######################
#train model
#pos = 31021
#neg = 286582
#ratio = 0.108
#neg_weight = len(x_train_pos0[0])/float(len(x_train_neg[1]))
#class weight
def make_w_list(w_rati0,list_list):
    list_out = []
    for i,list0 in enumerate(list_list):
        len0 = len(list0)
        weight0 = w_rati0[i]/float(len0)
        list_out.extend([weight0]*len0)
    return list_out

path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'

weight_list = make_w_list([80000,100000,800000],[x_train_pos0[0],x_train_neg[0],x_train_iedb[0]])
print(weight_list[0])
print(weight_list[-1])
#to save records
weight_name = 'rnn_wopretrain_weightv2_d0.3.h5'
record_file = path_save + 'record_trainingv2_d0.3.txt'
#parameters
n_iteration = 50
nb0 = 1
vb0 = 0
#auc_best = 0.80
batch0 = 128*2
for i in range(0,n_iteration):
    if isfile(record_file):
        file_write = open(record_file,'a')
    else:
        file_write = open(record_file,'a')
    file_write.write(str(i)+'\n')
    #perparing data
    x_pos_cycle = get_pos_for_fit(model_merge,x_train_pos0,x_train_pos1)
    x_train = [[],[]]
    x_train[0] = list(x_pos_cycle[0]) + list(x_train_neg[0]) + list(x_train_iedb[0])
    x_train[1] = list(x_pos_cycle[1]) + list(x_train_neg[1]) + list(x_train_iedb[1])
    print(len(x_train[0]),len(y_train))
    #fitting
    model_merge.fit([np.array(x_train[0]),np.array(x_train[1])],
                    np.array(y_train),
                    batch_size=batch0,
                    verbose=vb0,
                    epochs=nb0,
                    sample_weight = np.array(weight_list)
                    ) #class_weight={1:1,0:neg_weight}
    
    #calculate performance on training
    file_write.write('Training: ')
    auc_train = check_model(model_merge, x_train_pos0,x_train_pos1,x_train_neg,batch0 = 1024*8,record_file=file_write)
    file_write.write('AUC='+str(auc_train)+'\n')
    #calculate performance on validation
    file_write.write('Validation: ')
    auc_val = check_model(model_merge, x_val_pos0,x_val_pos1,x_val_neg,batch0 = 1024*8,record_file=file_write)
    file_write.write('AUC='+str(auc_val)+'\n')
    file_write.close()
    #save if better
    if auc_val > auc_best+0.01:
        auc_best = auc_val
        model_merge.save_weights(path_save+ weight_name,
                         overwrite=True)
        print('best AUC='+str(auc_best))




