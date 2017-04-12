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
from fileinput import 
from keras.models import Sequential
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent, merge
from keras.callbacks import ModelCheckpoint
#import keras.kernel_constraint
from keras.models import model_from_json
from scipy.stats import pearsonr
from keras.constraints import maxnorm
from keras.regularizers import l2


############################default value##############################
##################import coding path and dictionaries#####################
#dictionary avaialbe:
#aa_21_sparse_encoding.dict
#Blosum50_sparse.dict
#Blosum50_only.dict
#Sparse_only.dict

class RNN_master():

    def __init__(self):
        self.path_dict = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
        self.dict_name = 'Blosum50_only.dict'
        self.b_shuffle = True
        self.loss_function0 = 'mse'
        self.vb0 = 1
        self.nb = 1
        self.n_iteration = 30
        self.l2_value = 0
        self.drop_out_c = 0
        self.help_nn = 32
        self.act_fun = 'relu'
        self.help_layer0 = 1
        self.input_info = ''
        self.mask0 = True
            #maximum length of peptide 
        self.MAXLEN = 26
        self.BATCH_SIZE = 1024
        self.train_file0 = ''
        self.val_file0 = ''
        self.node0 = 32
        self.v1 = ''
        self.out_name = ''
        self.fix_node = 32

    #input file path and parameters from the setting file
    def get_input(self,parameter_file):
        for line0 in fileinput.input(parameter_file):
            self.input_info = input_info + line0
            line0 = line0.rstrip()
            print(line0)
            #ideally add a line to remove spaces in each line
            part1 = line0.split('=')[0]
            #print(type(part1))
            part2 = line0.split('=')[1]
            #print(part2)
            if 'path_data' in part1:
                self.path_data = part2
            if 'path_save' in part1:
                self.path_save = part2
            if 'train_file_name' in part1:
                self.train_file0 = part2
                #print(data_file_name)
            if 'val_file_name' in part1:
                self.val_file0 = part2
            if 'performance_file_name' in part1:
                self.performance_file_name = part2
            if 'version' in part1:
                self.v1 = part2
            if 'shuffle' in part1:
                if 'alse' in part2:
                    self.b_shuffle = False
            if  'neuron' in part1:
                self.node0 = int(part2)
            if 'out_name' in part1:
                self.out_name = part2
            if 'nb' == part1:
                self.nb0 = int(part2)
            if 'vb' == part1:
                self.vb0 = int(part2)
            if 'layer' in part1:
                self.LAYERS = int(part2)
            if 'loss_function' in part1:
                self.loss_function0 = part2
            if 'iteration' in part1:
                self.n_iteration = int(part2)
            if 'encoding' in part1:
                self.dict_name = part2
            if 'fix_node' in part1:
                self.fix_node = int(part2)
            if 'help_nn' in part1:
                self.help_nn = int(part2)
            if 'l2_value' in part1:
                self.l2_c = float(part2)
            if 'drop_out' in part1:
                self.drop_out_c = float(part2)
            if 'activation' in part1:
                self.act_fun = part2
            if 'help_layer' in part1:
                self.help_layer0 = int(part2)
            if 'batch' in part1:
                self.BATCH_SIZE=int(part2)
            if 'masking' in part1:
                self.mask0 = 'rue' in part2.lower()
                print('Masking='+str(mask0)) 
            if 'feature_dict' in part1:
                self.feature_dict0 = part2   
            if 'max_norm' in part1:
                self.constrain_max = float(part2)
        self.dict_aa = pickle.load(open(path_dict+dict_name,'r'))
        self.chars = dict_aa['A']
        self.dict_aa['-'] = np.zeros(len(chars))
    ###determine the encoding size
    ##########################construct input file name################ 
    ##training and validation file name
    print(self.train_file0)
    print(self.val_file0)
    #note_label = 'val_note.txt'
    #note_file0 = data_file_name+v1+note_label
    self.performance_file_name= self.performance_file_name +self.v1+self.out_name
    self.file_name0 = self.train_file0+self.v1


    def make_model(self,len_feature):
        ##########################Parameters for the model and dataset
        #TRAINING_SIZE = len(inputs)
        # Try replacing JZS1 with LSTM, GRU, or SimpleRNN
        HIDDEN_SIZE = self.node0
        RNN = recurrent.LSTM(HIDDEN_SIZE, input_shape=(None, len(self.chars)),
                              return_sequences=False,kernel_regularizer=l2(self.l2_c),
                              bias_regularizer=l2(self.l2_c),recurrent_dropout=self.drop_out_c,
                              dropout=self.drop_out_c)
        #len0_hla = 34
        
        #ratio_t = 1
        ###class number = binder or non-binder (1 = binder, 0 = non-binder)
        #classes = [0,1]
        
        ##########################start a model##########################
        ##########fixed part
        model_fixed = Sequential()
        model_fixed.add(Dense(help_nn,input_dim=len_feature,
                              activation=act_fun, kernel_constraint=maxnorm(constrain_max)))
        model_fixed.add(Dropout(drop_out_c))
        
        ##########recurrent part
        model_r = Sequential()
        if mask0:
            model_r.add(Masking(mask_value=0., input_shape=(MAXLEN, len(dict_aa['A']))))
        model_r.add(RNN)
        
        ####merge
        merged = Merge([model_fixed, model_r],mode='concat') 
        ###final
        final_model = Sequential()
        final_model.add(merged)
        for _ in range(0,help_layer0):
            final_model.add(Dense(help_nn, kernel_constraint=maxnorm(constrain_max)))
            final_model.add(Activation(act_fun))
            final_model.add(Dropout(drop_out_c))
        final_model.add(Dense(1))
        final_model.compile(loss=loss_function0, optimizer="adam")
        model = final_model
        json_string = model.to_json()
        open(path_save+file_name0+out_name+'_model.json', 'w').write(json_string)
        return model

#encoding will take a string or char, string=sequence and to return a matrix of encoded peptide sequence
#char = class, '0' = non-binding (0,1), '1' = binding (1,0)
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

def read_data(path_file0):
    #read
    X_0 = []
    y_0 = []
    max_len0 = 0
    for line0 in open(path_file0,'r'):
        line0 = line0.rstrip().split('\t')
        X_0.append(line0[0]+line0[1])
        y_0.append(float(line0[2]))
        max_len0 = max(max_len0,len(line0[0]+line0[1]))
    return [X_0,y_0,max_len0]

def encoding_data(list0,MAXLEN):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, len(chars)))
    X_encoded = encoding(X_0_m,list0,MAXLEN)
    return X_encoded

def split_x(x0,n0):
    x_fixed = x0[:,:n0,:].reshape((x0.shape[0],n0*len(chars)))
    x_variable = x0[:,n0:,:]
    return [x_fixed,x_variable]
 