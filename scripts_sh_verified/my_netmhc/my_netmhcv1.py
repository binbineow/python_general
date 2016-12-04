# -*- coding: utf-8 -*-
#Compared to v2 in the same folder
#this version only 

#####This learning script is built upon the previous HLA_RNN v13
#####Each iteration, the scirpt will calculate precision and recall of the training and whole validation set (also F1)
#####and non_identical peptide recall and non_substring peptide recall in the validation set (compared to training)
#####The performance writes into a file each iteration
#####

# how to use the model elsewhere...
#model = model_from_json(open('my_model_architecture.json').read())
#model.load_weights('my_model_weights.h5')

#####################import#################################
from __future__ import print_function
from keras.models import Sequential, Model
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from keras.callbacks import ModelCheckpoint
#from utilities import *
from keras.models import model_from_json
from sklearn.metrics import roc_auc_score,roc_curve
import cPickle as pickle
import numpy as np
import fileinput
#from keras.regularizers import l1,activity_l1

######Path for data as well as performance output are read in from fileinput ###
#no space
'''
path_data='/home/stanford/rbaltman/users/bchen45/data/MCL_data/'
path_save='/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/'
data_file_name='hla_ii_train_val_'
performance_file_name='model_performance'
version='_generalv1_x'
shuffle=False
path_data=/home/stanford/rbaltman/users/bchen45/data/MCL_data/
path_save=/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/
data_file_name=hla_ii_train_val
performance_file_name=model_performance
version=_generalv1_x_
shuffle=False
'''

############################default value##############################
##################import coding path and dictionaries#####################
path_dict = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#dictionary avaialbe:
#aa_21_sparse_encoding.dict
#Blosum50_sparse.dict
#Blosum50_only.dict
#Sparse_only.dict
dict_name = 'BLOSUM50_20long.dict'
b_shuffle = True
loss_function0 = 'categorical_crossentropy'
vb0 = 1
nb0 = 3
n_iteration = 100
ratio_t = 0.5
help_nn = 60
#input file path and parameters from the setting file
for line0 in fileinput.input():
    line0 = line0.rstrip()
    print(line0)
    #ideally add a line to remove spaces in each line
    part1 = line0.split('=')[0]
    #print(type(part1))
    part2 = line0.split('=')[1]
    #print(part2)
    if 'path_data' in part1:
        path_data = part2
    if 'path_save' in part1:
        path_save = part2
    if 'data_file_name' in part1:
        data_file_name = part2
        #print(data_file_name)
    if 'performance_file_name' in part1:
        performance_file_name = part2
    if 'version' in part1:
        v1 = part2
    if 'shuffle' in part1:
        if 'alse' in part2:
            b_shuffle = False
    if  'neuron' in part1:
        node0 = int(part2)
    if 'out_name' in part1:
        out_name = part2
    if 'nb' == part1:
        nb0 = int(part2)
    if 'vb' == part1:
        vb0 = int(part2)
    if 'layer' in part1:
        LAYERS = int(part2)
    if 'loss_function' in part1:
        loss_function0 = part2
    if 'iteration' in part1:
        n_iteration = int(part2)
    if 'encoding' in part1:
        dict_name = part2
    if 'ratio' == part1:
        ratio_t = float(part2)
        print(str(ratio_t))
    if 'help_nn' in part1:
        help_nn = int(part2)
        
 


# load blossom chemical encoding dictionary
chem = pickle.load(open(path_dict+dict_name,'r'))
# set '-' to all zeros
chem['-'] = np.zeros(len(chem['A']))
chem['X'] = np.zeros(len(chem['A']))


   
##########################construct input file name####################  
file_name0 = data_file_name+v1+'.txt'
performance_file_name= performance_file_name +v1+out_name



##########################Parameters for the model and dataset
#TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN


BATCH_SIZE = 128
#will play with Layers 
###class number = binder or non-binder (1 = binder, 0 = non-binder)
classes = [0,1]
    

#save the model
#json_string = model.to_json()
#open(path_save+file_name0+'_model.json', 'w').write(json_string)

#encoding will take a string or char, string=sequence and to return a matrix of encoded peptide sequence
#char = class, '0' = non-binding (0,1), '1' = binding (1,0)

# encode string of peptides with chemical information
def encode_seq(seq):
  return np.array([chem[k] for k in seq]).flatten()

            
# generate sliding windows over core for netmhc technique
# input: string of peptides
# output: yield each netmhc vector encoding
def make_windows(seq):
  #print(seq)
  first = seq[0]
  last = seq[-1]
  for i in range(0,len(seq)-9+1):
    core = seq[i:i+9]
    num_begin = i
    num_end = len(seq)-(i+9)
    nums = np.array([num_begin, 1 - num_begin, num_end, 1 - num_end, len(seq), 1-len(seq)])
    yield np.concatenate((encode_seq(first+core+last),nums))

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

def calf1(str1,str2):
    pre0 = float(str1)
    recall0 = float(str2)
    f1_out = 2.0*pre0*recall0/(pre0+recall0)
    return str(f1_out)

# select strongest binding core for each example
# input: x data matrix, y data, tracking info
# output: strongest bindings for each x datapoint, corresponding y, corresponding outputs
def select_best(x,y,track,model):
  X_new, y_new, outs_new = [], [], []
  outs = model.predict(x, batch_size=BATCH_SIZE,verbose=vb0)[:,1]
  #print(outs[0:100])
  for k,vs in track.items():
    best_max = -float("inf")
    for j in range(vs[0],vs[1]):
      if best_max < outs[j]:
        best_max = outs[j]
        selected_x = x[j]
        selected_y = y[j]
    X_new.append(selected_x)
    y_new.append(selected_y)
    outs_new.append(best_max)
  return np.array(X_new), np.array(y_new), np.array(outs_new)


# load data from files on netmhc website
# third column, binding strength
# input: file
# output: x_data matrix, y_data matrix, tracking info dictionary
def encoding(data_x,data_y):
    X, y = [], []
    track = {}
    for i,line in enumerate(data_x):
        # track which original datapoint each window belongs to
        begin, end = len(X), len(X)
        for s in make_windows(line):
            X.append(s)
            y.append(data_y[i])
            end += 1
            # data in the range of (begin, end) belongs to datapoint i
        track[i] = (begin, end)
    return np.array(X), np.array(y), track
  
#main function
def main():
    #initiate 
    inputs=[]
    outputs=[]
    char_set = set([' '])
    class_set = set()
    max_len = 0
    X_train = []
    y_val = []
    X_val = []
    y_train = []
    y_val = []
    y_val_linear = []
    
    #file_name0 ='HLADRB10401simplev1_tr_1_val.csv'
    for line in open(path_data+file_name0,'r'):
        in_,out_ = [x.rstrip() for x in line.split("\t")]
        if len(out_) != 1:
            raise Exception("Output should be single characer")
        else:
            if out_ == '0' :
                X_train.append(in_)
                y_train.append([1,0])
            elif out_ == '1':
                X_train.append(in_)
                y_train.append([0,1])
            else:
                out_ = str(int(out_) -2)
                if out_ == '0':
                    X_val.append(in_)
                    y_val.append([1,0])
                    y_val_linear.append(0)
                else:
                    X_val.append(in_)
                    y_val.append([0,1])
                    y_val_linear.append(1)
              
            max_len = max([max_len,len(in_),len(out_)])
            inputs.append(in_)
            outputs.append(out_)
            
            #for c in in_: char_set.add(c)
            class_set.add(out_)
        
    y_val_linear = np.array(y_val_linear)
    print(y_val_linear[0:100])
    
    #create training or validation matrix
    X_train,y_train_encoded,track = encoding(X_train, y_train)
    X_test,y_test_encoded,track_test = encoding(X_val, y_val)

    ##########################start a model
    model = Sequential()
    model.add(Dense(help_nn, input_shape=(X_train.shape[1],)))
    model.add(Activation('relu'))
    classes = [0,1]
    model.add(Dense(len(classes)))
    model.add(Activation('softmax'))
    model.compile(loss=loss_function0, optimizer='adam')
        
    print(X_train[1,:])
    print('Training='+str(len(X_train))+' Validation='+str(len(X_val)))
    print("Input loaded ")

    #training
    for iteration in range(0, n_iteration):
        #iterations.append(str(iteration))
        print()
        print('-' * 50)
        print('Iteration', iteration)
        X_new, y_new, outs_new = select_best(X_train, y_train_encoded, track, model)
        print(len(X_new))
        model.fit(X_new, y_train, batch_size=BATCH_SIZE, nb_epoch=nb0,verbose=vb0)
        _, y_test_select, test_outs = select_best(X_test, y_test_encoded, track_test,model)
        print(test_outs[0:100])        
        auc_val = roc_auc_score(y_val_linear, test_outs)
        #print('Val_Precision='+str(float(tp0)/(tp0+fp0)))
        #print('Val_Recall='+str(float(tp0)/(tp0+fn0)))
        print('\n')
        print([iteration,auc_val])
        #print('[iteration,train_pre,train_recall,train_f1,val_pre,val_recall,val_f1,recall_non_i,recall_non_sub]')
        #print([iteration,train_pre,train_recall,train_f1,val_pre,val_recall,val_f1,recall_non_i,recall_non_sub])
        #output_perf2([iteration,train_pre,train_recall,train_f1,val_pre,val_recall,val_f1,recall_non_i,recall_non_sub])
        model.save_weights(path_save+file_name0+out_name+'_weight.h5',overwrite=True)
    pickle.dump(roc_curve(y_val_linear,test_outs),open(path_save+'roc_curve_data+'+out_name+'.list','w+'))
    
main()
