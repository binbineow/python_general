# -*- coding: utf-8 -*-
#
###with regulizer and drop out
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
#from keras.regularizers import l2,activity_l2

############################default value##############################
##################import coding path and dictionaries#####################
path_dict = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#dictionary avaialbe:
#aa_21_sparse_encoding.dict
#Blosum50_sparse.dict
#Blosum50_only.dict
#Sparse_only.dict
dict_name = 'Blosum50_only.dict'
b_shuffle = True
loss_function0 = 'mse'
vb0 = 0
nb = 3
n_iteration = 30
l2_value = 0
drop_out_c = 0
help_nn = 0
act_fun = 'tanh'
help_layer0 = 1
input_info = ''

#input file path and parameters from the setting file
for line0 in fileinput.input():
    input_info = input_info + line0
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
    if 'help_nn' in part1:
        help_nn = int(part2)
    if 'l2_value' in part1:
        l2_c = float(part2)
    if 'drop_out' in part1:
        drop_out_c = float(part2)
    if 'activation' in part1:
        act_fun = part2
    if 'help_layer' in part1:
        help_layer0 = int(part2)
    if 'batch' in part1:
        BATCH_SIZE=int(part2)
        
        


dict_aa = pickle.load(open(path_dict+dict_name,'r'))
chars = dict_aa['A']
dict_aa['-'] = np.zeros(len(chars))
###determine the encoding size

   
##########################construct input file name################ 
##training and validation file name
print(train_file0)
print(val_file0)
#note_label = 'val_note.txt'
#note_file0 = data_file_name+v1+note_label
performance_file_name= performance_file_name +v1+out_name

##########################Parameters for the model and dataset
#TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN
HIDDEN_SIZE = node0
BATCH_SIZE = 1024
RNN = recurrent.LSTM(HIDDEN_SIZE, input_shape=(None, len(chars)), return_sequences=False,W_regularizer=l2(l2_c),b_regularizer=l2(l2_c),dropout_W=drop_out_c,dropout_U=drop_out_c)
len0_hla = 34

#ratio_t = 1
###class number = binder or non-binder (1 = binder, 0 = non-binder)
#classes = [0,1]


##########################start a model##########################
##########fixed part
model_fixed = Sequential()
model_fixed.add(Dense(help_nn,input_dim=len0_hla*len(chars),activation=act_fun))

##########recurrent part
model_r = Sequential()
model_r.add(RNN)
      
####merge
merged = Merge([model_fixed, model_r],mode='concat') 
###final
final_model = Sequential()
final_model.add(merged)
for _ in range(0,help_layer0):
    final_model.add(Dense(help_nn))
    final_model.add(Activation(act_fun))
final_model.add(Dense(1))
final_model.compile(loss=loss_function0, optimizer="adam")
model = final_model
json_string = model.to_json()
open(path_save+file_name0+out_name+'_model.json', 'w').write(json_string)

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
    

def main():
    MAXLEN = 26
    [X_train_fixed,train_seq,y_train] = pickle.load(open(path_data+train_file0,'r'))
    [X_val_fixed,val_seq,y_val] = pickle.load(open(path_data+tval_file0,'r'))
    ########encoding
    X_train_variable = encoding_data(train_seq,MAXLEN)
    X_val_variable = encoding_data(val_seq,MAXLEN)
    r_best = 0
    output_perf2(['Iteration','Training PCC','Training p-val','Val PCC','Val p-val'])
    print('start training')
    for n0 in range(0,n_iteration+1):
        #fit    
        #print(y_train)
        model.fit([X_train_fixed,X_train_variable], y_train, batch_size=BATCH_SIZE, verbose=vb0, nb_epoch=nb0)      
        #calculate the performance
        #calculate Pearson Correltion Coeficient 
        y_train_pred = model.predict([X_train_fixed,X_train_variable],batch_size=BATCH_SIZE)
        y_train_pred = y_train_pred.reshape(y_train_pred.shape[0])
        [r0_train,pval0_train] = pearsonr(y_train_pred,y_train)
        
        y_predicted = model.predict([X_val_fixed,X_val_variable],batch_size=BATCH_SIZE)
        y_predicted = y_predicted.reshape(y_predicted.shape[0])
        print(y_predicted)
        [r0, pval0] = pearsonr(y_predicted,y_val)
        #print('PCC in validation'+str(r0))
        #save performance
        output_perf2([n0,r0_train,pval0_train,r0,pval0])
        #print performance
        #print([n0,r0_train,pval0_train,r0,pval0])
        #print('Predicted binding aff')
        #print(y_predicted[0:10])
        #print('Measured binding aff')
        #print(y_val[0:10])
        #save the model
        if r0 > r_best+0.005:
            model.save_weights(path_save+file_name0+out_name+'_weight.h5',overwrite=True)
            r_best = r0
#print out parameteres
output_perf2([input_info]) 
main()
