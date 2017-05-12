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
from os.path import isfile
print(keras.__version__)


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
        self.path_save = 'path_save=/home/stanford/rbaltman/users/bchen45/results/HLA_pred_general_model/'
        self.b_shuffle = True
        self.loss_function0 = 'categorical_crossentropy'
        self.vb0 = 1
        self.nb0 = 1
        self.fix_len = 0
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
        self.drop_r = 0

        
    def process_line(self,line0):
        self.input_info = self.input_info + line0
        if len(line0) > 1:
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
                #print('Masking='+str(self.mask0)) 
            if 'feature_dict' in part1:
                self.feature_dict0 = part2   
            if 'max_norm' in part1:
                self.constrain_max = float(part2)
            if 'drop_r' in part1:
                self.drop_r=float(part2)
        

    #input file path and parameters from the setting file
    def get_input(self,parameter_file):
        if parameter_file != '':
            for line0 in open(parameter_file,'rU'):
                self.process_line(line0)
        else:
            for line0 in fileinput.input():
                self.process_line(line0)
        self.dict_aa = pickle.load(open(self.path_dict+self.dict_name,'r'))
        self.chars = self.dict_aa['A']
        self.dict_aa['B'] = self.dict_aa['D'] 
        self.dict_aa['Z'] = self.dict_aa['Q']
        self.dict_aa['-'] = np.zeros(len(self.chars))
    
    def create_out_file_names(self):
        ##training and validation file name
        print(self.train_file0)
        print(self.val_file0)
        #note_label = 'val_note.txt'
        #note_file0 = data_file_name+v1+note_label
        self.performance_file_name= self.performance_file_name +self.v1+self.out_name
        self.file_name0 = self.train_file0+self.v1


    def make_model(self,model_type='class2'):
        ##########################Parameters for the model and dataset
        #input layers
        seq_input = Input(shape=(self.MAXLEN,len(self.dict_aa['A'])))
        fix_input = Input(shape=(self.fix_len,))
        #set RNN layer
        if self.mask0:
            seq_input0 = Masking(mask_value=0.0)(seq_input)
        else:
            seq_input0 = seq_input
        rnn0 = recurrent.LSTM(self.node0, activation=self.act_fun, #recurrent_activation=self.act_fun,
                                    use_bias=True, kernel_initializer='glorot_uniform',
                                    recurrent_initializer='orthogonal', bias_initializer='zeros', 
                                    unit_forget_bias=True, kernel_regularizer=None, 
                                    recurrent_regularizer=None, bias_regularizer=None, 
                                    activity_regularizer=None, kernel_constraint=None, 
                                    recurrent_constraint=None, bias_constraint=None, dropout=self.drop_out_c, 
                                    recurrent_dropout=self.drop_r)(seq_input0)
        fix0 = Dense(self.fix_node,activation=self.act_fun)(fix_input)
        fix0_d = Dropout(self.drop_out_c)(fix0)
        merge_layer = concatenate([rnn0,fix0_d])
        combine1 = Dense(self.help_nn,activation=self.act_fun)(merge_layer)
        combine1_d = Dropout(self.drop_out_c)(combine1)
        combine2 = Dense(self.help_nn,activation=self.act_fun)(combine1_d)
        combine2_d = Dropout(self.drop_out_c)(combine2)
        if model_type == 'regression':
            dense0 = Dense(1)(combine2_d)
            final_model = Model(inputs = [seq_input,fix_input],outputs = [dense0])
            final_model.compile(loss=self.loss_function0,optimizer='adam')
        elif 'class' in model_type:
            class_n = int(model_type[-1])
            dense0 = Dense(class_n,activation='softmax')(combine2_d)
            final_model = Model(inputs = [fix_input,seq_input],outputs = [dense0])
            final_model.compile(loss=self.loss_function0, optimizer="RMSprop")

        json_string = final_model.to_json()
        open(self.path_save+self.file_name0+self.out_name+'_model.json', 'w').write(json_string)
        return final_model
    


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

# def split_x(x0,n0):
#     x_fixed = x0[:,:n0,:].reshape((x0.shape[0],n0*len(chars)))
#     x_variable = x0[:,n0:,:]
#     return [x_fixed,x_variable]



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


#testing
path_para = '/home/stanford/rbaltman/users/bchen45/code/slurm/'
para_file = 'run_now_merge_class.txt'
rnn_master = RNN_master()
rnn_master.get_input(path_para+para_file)
rnn_master.create_out_file_names()
rnn_master.mask0
rnn_master.drop_out_c = 0.45
# [mhc_train,seq_train,label_train] = pickle.load(open(rnn_master.path_data+rnn_master.train_file0,'r'))
# [mhc_val,seq_val,label_val] = pickle.load(open(rnn_master.path_data+rnn_master.val_file0,'r'))
# #get max length
# rnn_master.MAXLEN = max([len(seq0) for seq0 in seq_train]+[len(seq0) for seq0 in seq_val])
# print(rnn_master.MAXLEN)
# #encode x 
# len_hla = len(mhc_train[0])
# x_mhc_train = encoding_fixed(mhc_train,len_hla,rnn_master.dict_aa,len(rnn_master.chars),add_placeh=3)
# x_seq_train = encoding_data(seq_train,rnn_master.MAXLEN,rnn_master.dict_aa,len(rnn_master.chars))
# x_mhc_val = encoding_fixed(mhc_val,len_hla,rnn_master.dict_aa,len(rnn_master.chars),add_placeh=3)
# x_seq_val = encoding_data(seq_val,rnn_master.MAXLEN,rnn_master.dict_aa,len(rnn_master.chars))
# #process y
# y_train_iedb = encoding_y(label_train)
# y_val_iedb = encoding_y(label_val)
# label_train = np.array([int(y0) for y0 in label_train])
# label_val = np.array([int(y0) for y0 in label_val])
# #generate model
# fix_len = len(x_mhc_train[0])
# rnn_master.fix_len = fix_len
rnn_master.fix_len = 801
rnn_master.MAXLEN = 26
rnn_master.self_node = 64
model_merge = rnn_master.make_model('class2')

##load weight
#model_merge.load_weights(path_save+ 'rnn_wopretrain_weightv5_d0.4.h5_8_12_8')
#rmsprop = keras.optimizers.RMSprop(lr=0.001)
#model_merge.compile(loss='categorical_crossentropy', optimizer=rmsprop)


###########get data from my MCL patient package
#get data
path_store = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/'
file_store = 'mcl_patient_mhc2.training_list1_v5'
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
dict_name = 'aa_21_sparse_encoding.dict'
path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'
MAXLEN = 26
dict_aa = pickle.load(open(path_encoding+dict_name,'r'))
[train_pos0,train_pos1,val_pos0,val_pos1,train_neg,val_neg, 
  list_random,list_ms_random,list_shuffle_random,
  list_anti_patient_train,list_anti_patient_val, iedb_train ,
  iedb_train_y ,iedb_val ,iedb_val_y] = pickle.load(open(path_store+file_store,'r'))

#get mhc length and dictionary length
len_mhc = len(train_pos0[0][0])
len_char = len(dict_aa['A'])

#make train_pos
train_pos = [[],[]]
train_pos[0] = train_pos0[0] + train_pos1[0]
train_pos[1] = train_pos1[1] + train_pos1[1]
#make val_pos
val_pos = [[],[]]
val_pos[0] = val_pos0[0] + val_pos1[0]
val_pos[1] = val_pos1[1] + val_pos1[1]

#including 04:01 data
iedb_train[0] = iedb_train[0] + iedb_val[0]
iedb_train[1] = iedb_train[1] + iedb_val[1]
iedb_train_y = iedb_train_y + iedb_val_y

#encoding training data
y_train_neg = make_and_encode_label(len(train_neg[0]),0)
y_train_pos = make_and_encode_label(len(train_pos[0]),1)
y_train_iedb = list(encoding_y(iedb_train_y))
y_train = y_train_pos + y_train_neg + y_train_iedb
x_train_pos = encode_apair(train_pos,dict_aa,len_mhc,MAXLEN,len_char)
#x_train_pos1 = encode_apair(train_pos1,dict_aa,len_mhc,MAXLEN,len_char)
x_train_neg = encode_apair(train_neg,dict_aa,len_mhc,MAXLEN,len_char)
x_train_iedb = encode_apair(iedb_train,dict_aa,len_mhc,MAXLEN,len_char)

#encoding validation data
y_val_neg = make_and_encode_label(len(val_neg),0)
y_val_pos = make_and_encode_label(len(val_pos),1)
#y_val = y_val_pos + y_val_neg
x_val_pos = encode_apair(val_pos,dict_aa,len_mhc,MAXLEN,len_char)
#x_val_pos1 = encode_apair(val_pos1,dict_aa,len_mhc,MAXLEN,len_char)
x_val_neg = encode_apair(val_neg,dict_aa,len_mhc,MAXLEN,len_char)
x_val_iedb = encode_apair(iedb_val,dict_aa,len_mhc,MAXLEN,len_char)

#x_val_iedb = encode_apair(iedb_val,dict_aa,len_mhc,MAXLEN,len_char


#encoding other lists
x_list_random = encode_apair(list_random,dict_aa,len_mhc,MAXLEN,len_char)
x_list_ms_random = encode_apair(list_ms_random,dict_aa,len_mhc,MAXLEN,len_char)
x_list_shuffle_random = encode_apair(list_shuffle_random,dict_aa,len_mhc,MAXLEN,len_char)
x_list_anti_patient_train = encode_apair(list_anti_patient_train,dict_aa,len_mhc,MAXLEN,len_char)
x_list_anti_patient_val = encode_apair(list_anti_patient_val,dict_aa,len_mhc,MAXLEN,len_char)


#supporting funcitons
#for merge model
def cal_performance(model0,x,y,batch0=1024,average0='macro'):
    from sklearn.metrics import roc_auc_score
    list_predicted = model0.predict(x,verbose=1,batch_size=batch0)[:,1]
    auc_val = roc_auc_score(y, list_predicted,average=average0)
    print('AUC='+str(auc_val))
    return auc_val

def cal_auc_pos_neg(pre_pos,pre_neg,average0='macro'):
    from sklearn.metrics import roc_auc_score
    list_scores = list(pre_pos) + list(pre_neg)
    list_true = [1]*len(pre_pos) + [0]*len(pre_neg)
    auc_val = roc_auc_score(list_true, list_scores,average=average0)
    print('AUC='+str(auc_val))
    return auc_val

def get_percent_pos(model0,x,cut_off=0.5,batch0=1024):
    list_scores = model0.predict(x,verbose=1,batch_size=batch0)[:,1]
    percent0 = np.sum(list_scores >= cut_off)/float(len(x[0]))
    print(percent0)
    return percent0
    
##checking the model on training data
def check_model(model_merge, x_train_pos, x_train_neg,batch0 = 1024,
                record_file=1,cut_off0=0.5,average0='macro'):
    list_predicted_pos = model_merge.predict(x_train_pos,verbose=1,batch_size=batch0)[:,1]
    #list_predicted_pos1 = model_merge.predict(x_train_pos1,verbose=1,batch_size=batch0)[:,1]
    list_predicted_neg = model_merge.predict(x_train_neg,verbose=1,batch_size=batch0)[:,1]
    #np.sum(list_predicted_neg < 0.5)/float(len(list_predicted_neg))
    #list_predict_pos = [max(x1,x2) for x1,x2 in zip(list_predicted_pos0,list_predicted_pos1)]
    #list_predict_scores = list_predict_pos + list_predicted_neg
    #list_true = [1]*len(list_predict_pos) + [0] * list_predicted_neg
    #print(np.sum(list_predict_pos >= 0.5)/float(len(list_predict_pos)))
    print('Sensitivity:')
    sens0 = np.sum(np.array(list_predicted_pos) >= cut_off0)/float(len(list_predicted_pos))
    print(sens0)
    print('Specificity:')
    spec0 = 1-np.sum(np.array(list_predicted_neg) >= cut_off0)/float(len(list_predicted_neg))
    print(spec0)
    if record_file != 1:
        record_file.write('Sensitivity='+str(sens0)+' Specificity='+str(spec0)+' ')
    #plt.scatter(list_predicted_pos0,list_predicted_pos1)
    auc_model = cal_auc_pos_neg(list_predicted_pos,list_predicted_neg,average0=average0)
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

def make_w_list(w_rati0,list_list):
    list_out = []
    for i,list0 in enumerate(list_list):
        len0 = len(list0)
        weight0 = w_rati0[i]/float(len0)
        list_out.extend([weight0]*len0)
    return list_out

#def generate_iedb_weight(iedb_mhc_list,iedb_label,total_val):

def make_w_list(w_rati0,list_list):
    list_out = []
    for i,list0 in enumerate(list_list):
        len0 = len(list0)
        weight0 = w_rati0[i]/float(len0)
        list_out.extend([weight0]*len0)

    return list_out



#check_model(model_merge,x_val_pos, x_val_neg,cut_off0=0.5,average0='macro')

#############training######################
#train model
#pos = 31021
#neg = 286582
#ratio = 0.108
#neg_weight = len(x_train_pos0[0])/float(len(x_train_neg[1]))
#class weight
path_save = '/cstor/stanford/rbaltman/users/bchen45/mcl_data/model_weight/'
rmsprop = keras.optimizers.RMSprop(lr=0.001)
model_merge.compile(loss='categorical_crossentropy', optimizer=rmsprop)

weight_list = make_w_list([100000,160000,80000],[x_train_pos[0],x_train_neg[0],x_train_iedb[0]])
#weight_list = make_w_list([80000,100000],[x_train_pos[0],x_train_neg[0]])
print(weight_list[0])
print(weight_list[-1])
#to save records
weight_name = 'rnn_wopretrain_weightv5_d0.45.h5.10_16_8'
record_file = path_save + 'record_trainingv5_d0.45.10_16_8.txt'
#parameters
n_iteration = 100
nb0 = 1
vb0 = 1
auc_best = 0.84
batch0 = 128*4
x_train = [[],[]]
x_train[0] = list(x_train_pos[0]) + list(x_train_neg[0]) + list(x_train_iedb[0])
x_train[1] = list(x_train_pos[1]) + list(x_train_neg[1]) + list(x_train_iedb[1])

#training
for i in range(0,n_iteration):
    if isfile(record_file):
        file_write = open(record_file,'a')
    else:
        file_write = open(record_file,'w+')
    file_write.write(str(i)+'\n')
    #perparing data
    #x_pos_cycle = x_pos
    #

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
    auc_train = check_model(model_merge, x_train_pos,x_train_neg,batch0 = 1024*8,record_file=file_write)
    file_write.write('AUC='+str(auc_train)+'\n')
    #calculate performance on validation
    file_write.write('Validation: ')
    auc_val = check_model(model_merge, x_val_pos,x_val_neg,batch0 = 1024*8,record_file=file_write)
    file_write.write('AUC='+str(auc_val)+'\n')
    file_write.close()
    #save if better
    if auc_val > auc_best+0.01:
        auc_best = auc_val
        model_merge.save_weights(path_save+ weight_name,
                         overwrite=True)
        print('best AUC='+str(auc_best))


