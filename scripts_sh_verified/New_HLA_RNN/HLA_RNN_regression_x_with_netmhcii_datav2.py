# -*- coding: utf-8 -*-
#

# how to use the model elsewhere...
#model = model_from_json(open('my_model_architecture.json').read())
#model.load_weights('my_model_weights.h5')

#####################import#################################
from __future__ import print_function
from keras.models import Sequential
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent, Merge
from keras.callbacks import ModelCheckpoint
from utilities import *
from keras.models import model_from_json
from scipy.stats import pearsonr
from keras.regularizers import l2, activity_l2
#from keras.regularizers import l1,activity_l1

############################default value##############################
##################import coding path and dictionaries#####################
path_dict = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#dictionary avaialbe:
#aa_21_sparse_encoding.dict
#Blosum50_sparse.dict
#Blosum50_only.dict
#Sparse_only.dict
dict_name = 'aa_21_sparse_encoding.dict'
b_shuffle = True
loss_function0 = 'mse'
vb0 = 0
nb = 3
n_iteration = 30
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
    if 'train_file_name' in part1:
        train_file0 = part2
        #print(data_file_name)
    if 'val_file_name' in part1:
        val_file0 = part2
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
        
        
 

dict_aa = pickle.load(open(path_dict+dict_name,'r'))
###determine the encoding size
chars = dict_aa['A']
   
##########################construct input file name################ 
file_name0 = train_file0+v1+'.txt'
#note_label = 'val_note.txt'
#note_file0 = data_file_name+v1+note_label
performance_file_name= performance_file_name +v1+out_name

##########################Parameters for the model and dataset
#TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN
RNN = recurrent.LSTM
HIDDEN_SIZE = node0
BATCH_SIZE = 128

#ratio_t = 1
###class number = binder or non-binder (1 = binder, 0 = non-binder)
#classes = [0,1]
    

##########################start a model##########################
'Need to fix the model for regression here'

model_fixed = Sequential()
model_fixed.add(Dense(HIDDEN_SIZE,input_dim=19*len(chars), W_regularizer=l2(0.01), activity_regularizer=activity_l2(0.01)))
model_fixed.add(Dropout(0.5))

model = Sequential()
# "Encode" the input sequence using an RNN, producing an output of HIDDEN_SIZE
#model.add(Masking())
model.add(RNN(HIDDEN_SIZE, input_shape=(None, len(chars)), return_sequences=True))
model.add(Dropout(0.5))
# for _ in xrange(LAYERS-1):
#     model.add(RNN(HIDDEN_SIZE, return_sequences=True))
#    #model.add(Dropout(0.5))
model.add(RNN(HIDDEN_SIZE, return_sequences=False))
model.add(Dropout(0.5))
merged = Merge([model_fixed, model], mode='concat')

final_model = Sequential()
final_model.add(merged)

final_model.add(Dense(1,W_regularizer=l2(0.01), activity_regularizer=activity_l2(0.01)))
final_model.compile(loss=loss_function0, optimizer="adam")
# model.add(Dense(len(classes)))
# model.add(Activation('softmax'))
# model.compile(loss=loss_function0, optimizer='adam')
#save the model
#json_string = model.to_json()
#open(path_save+file_name0+'_model.json', 'w').write(json_string)

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
        X_0.append(line0[0])
        y_0.append(float(line0[1]))
        max_len0 = max(max_len0,len(line0[0]))
    return [X_0,y_0,max_len0]

def encoding_data(list0,MAXLEN):
    #encoding   
    X_0_m = np.zeros((len(list0), MAXLEN, len(chars)))
    X_encoded = encoding(X_0_m,list0,MAXLEN)
    return X_encoded

        
    

def main():
    MAXLEN = 0
    [X_train, y_train,maxlen0] = read_data(path_data+train_file0)
    
    y_train =np.array(y_train)#.reshape((len(y_train),1))
    
    MAXLEN = max(MAXLEN,maxlen0)
        #shuffle if indicated
    if b_shuffle:
        [X_train,y_train] = shuffle_train(X_train, y_train)
        print('after shuffling, len(x)='+str(len(X_train))) 
    else:
        print('without shuffling, len(x)='+str(len(X_train)))
    [X_val, y_val,maxlen0] = read_data(path_data+val_file0)
    
    y_val = np.array(y_val)
    MAXLEN = max(MAXLEN,maxlen0)
    
    ########encoding
    X_train = encoding_data(X_train,MAXLEN)
    #print(X_train)
    X_val = encoding_data(X_val,MAXLEN)
    y_train = np.array(y_train)
    
    x_train_fixed = X_train[:,:19,:].reshape((X_train.shape[0],19*len(chars)))
    x_train_variable = X_train[:,19:,:]
    
    x_val_fixed = X_val[:,:19,:].reshape((X_val.shape[0],19*len(chars)))
    x_val_variable = X_val[:,19:,:]
    
    print(X_train.shape)
    print(x_train_fixed.shape, x_train_variable.shape)
    
    for n0 in range(0,n_iteration+1):
        #fit    
        print(y_train)
        final_model.fit([x_train_fixed, x_train_variable], y_train, batch_size=BATCH_SIZE, verbose=vb0, nb_epoch=nb0,validation_data=([x_val_fixed,x_val_variable], y_val))      
        #calculate the performance
        #calculate Pearson Correltion Coeficient 
        y_predicted = final_model.predict([x_val_fixed,x_val_variable])
        y_predicted = y_predicted.reshape(y_predicted.shape[0])
        print(y_predicted)
        [r0, pval0] = pearsonr(y_predicted,y_val)
        #save performance
        output_perf2([n0,r0,pval0])
        #print performance
        print([n0,r0,pval0])
        print('Predicted binding aff')
        print(y_predicted[0:100])
        print('Measured binding aff')
        print(y_val[0:100])
        #save the model
        model.save_weights(path_save+file_name0+out_name+'_weight.h5',overwrite=True)

main()
