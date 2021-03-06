# -*- coding: utf-8 -*-
#does not include chemistry information

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
from keras.models import Sequential
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from keras.callbacks import ModelCheckpoint
from utilities import *
from keras.models import model_from_json
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

#default value
b_shuffle = True
loss_function0 = 'categorical_crossentropy'
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
        
 
##################import coding path and dictionaries#####################
path_dict = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
#dictionary avaialbe:
#aa_21_sparse_encoding.dict
#Blosum50_sparse.dict
#Blosum50_only.dict
#Sparse_only.dict
note_label = 'val_note.txt'
dict_name = 'aa_21_sparse_encoding.dict'
dict_aa = pickle.load(open(path_dict+dict_name,'r'))
   
##########################construct input file name####################  
file_name0 = data_file_name+v1+'.txt'
note_file0 = data_file_name+v1+note_label
performance_file_name= performance_file_name +v1+out_name

#########################construct note label############################
list0 = []
for line0 in open(path_data+note_file0,'r'):
    num0 = float(line0)
    list0.append(num0)
#this array has 1,2,3 to distinghish three types of positive peptides 
list_val0 = np.array(list0)
mask_non_i = list_val0 >= 2
# non_i includes non_sub
len_non_i = sum(mask_non_i)
mask_non_sub = list_val0 == 3
len_non_sub = sum(mask_non_sub)





##########################Parameters for the model and dataset
#TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN
RNN = recurrent.LSTM
HIDDEN_SIZE = node0
BATCH_SIZE = 128
#will play with Layers 
ratio_t = 1
chars = 'ACDEFGHIKLMNPQRSTVWXY'#'0123456789+ '
if dict_name == 'Blosum50_sparse.dict':
    chars = chars + chars
classes = [0,1]
    

##########################start a model
model = Sequential()
# "Encode" the input sequence using an RNN, producing an output of HIDDEN_SIZE
#model.add(Masking())
model.add(RNN(HIDDEN_SIZE, input_shape=(None, len(chars)), return_sequences=True))
for _ in xrange(LAYERS-1):
    model.add(RNN(HIDDEN_SIZE, return_sequences=True))
#    #model.add(Dropout(0.5))
model.add(RNN(HIDDEN_SIZE, return_sequences=False))
model.add(Dense(len(classes)))
model.add(Activation('softmax'))
model.compile(loss=loss_function0, optimizer='adam')
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

def calf1(str1,str2):
    pre0 = float(str1)
    recall0 = float(str2)
    f1_out = 2.0*pre0*recall0/(pre0+recall0)
    return str(f1_out)

  
#main function
for _ in range(0,1):
    #initiate 
    inputs=[]
    outputs=[]
    char_set = set([' '])
    class_set = set()
    max_len = 0
    X_train = []
    X_train_p = []
    X_train_n = []
    X_val_p = []
    X_val_n = []
    y_train = []
    y_val_p = []
    y_val_n = []
    
    #file_name0 ='HLADRB10401simplev1_tr_1_val.csv'
    for line in open(path_data+file_name0,'r'):
        in_,out_ = [x.rstrip() for x in line.split("\t")]
        if len(out_) != 1:
            raise Exception("Output should be single characer")
        else:
            if out_ == '0' :
                X_train_n.append(in_)
                X_train.append(in_)
                y_train.append(out_)
            elif out_ == '1':
                X_train_p.append(in_)
                X_train.append(in_)
                y_train.append(out_)
            else:
                out_ = str(int(out_) -2)
                if out_ == '0':
                    X_val_n.append(in_)
                    y_val_n.append(out_)
                else:
                    X_val_p.append(in_)
                    y_val_p.append(out_)
              
            max_len = max([max_len,len(in_),len(out_)])
            inputs.append(in_)
            outputs.append(out_)
            
            #for c in in_: char_set.add(c)
            class_set.add(out_)
    #
    
    #shuffle if indicated
    if b_shuffle:
        [X_train,y_train] = shuffle_train(X_train, y_train)
        print('after shuffling, len(x)='+str(len(X_train)))      
     
    #creating encoding table
    print(class_set)
    classes = ''.join(class_set)
    #ctable = CharacterTable(chars, MAXLEN)
    #classtable = CharacterTable(classes, 1)
    MAXLEN = max_len #DIGITS + 1 + DIGITS

    #create training or validation matrix
    X_train_m = np.zeros((len(X_train), MAXLEN, len(chars)), dtype=np.bool)
    X_val_p_m = np.zeros((len(X_val_p), MAXLEN, len(chars)), dtype=np.bool)
    X_val_n_m = np.zeros((len(X_val_n), MAXLEN, len(chars)), dtype=np.bool)
    X_train_p_m = np.zeros((len(X_train_p), MAXLEN, len(chars)), dtype=np.bool)
    X_train_n_m = np.zeros((len(X_train_n), MAXLEN, len(chars)), dtype=np.bool)
    y_train_m = np.zeros((len(y_train), len(classes)), dtype=np.bool)
    y_val_p_m = np.zeros((len(y_val_p), len(classes)), dtype=np.bool)
    y_val_n_m = np.zeros((len(y_val_n), len(classes)), dtype=np.bool)
    
    X_train = encoding(X_train_m, X_train,MAXLEN)
    X_train_p = encoding(X_train_p_m, X_train_p,MAXLEN)
    X_train_n = encoding(X_train_n_m, X_train_n,MAXLEN)
    X_val_p = encoding(X_val_p_m, X_val_p,MAXLEN)
    X_val_n = encoding(X_val_n_m, X_val_n,MAXLEN)
    y_train = encoding(y_train_m, y_train,1)
    y_val_p = encoding(y_val_p_m, y_val_p,1)
    y_val_n = encoding(y_val_n_m, y_val_n,1)
    
    X_val = np.concatenate((X_val_n,X_val_p))
    y_val = np.concatenate((y_val_n,y_val_p))
    print('Training='+str(len(X_train))+' Validation='+str(len(X_val)))
    print("Input loaded ")
    

    
    #Create checkpoint
    #checkpointer = ModelCheckpoint(filepath=model_name+'.weight', verbose=1, save_best_only=True)
    # Train the model each generation and show predictions against the validation dataset
    '''
    if os.path.isfile(path_save+performance_file_name+v1+out_name+'.csv'):     
        file_out = open(path_save+performance_file_name+v1+out_name+'.csv','a')
    else:
        file_out = open(path_save+performance_file_name+v1+out_name+'.csv','w+')
    '''
    #iterations = []
    f2_val_best = []
    n_best = []
    ptotal0 = len(X_train_p)
    ntotal0 = len(X_train_n)
    training_n = str(ptotal0+ntotal0)
    for iteration in range(1, n_iteration):
        #iterations.append(str(iteration))
        print()
        print('-' * 50)
        print('Iteration', iteration)
        
        model.fit(X_train, y_train, batch_size=BATCH_SIZE, verbose=vb0, nb_epoch=nb0, class_weight={1:1,0:1.0/ratio_t/2})      
        #####predicting training
        ptotal0 = len(X_train_p)
        print('p_training='+str(ptotal0))
        ntotal0 = len(X_train_n)
        #print('Train_Postive')
        #print(model.predict_classes(X_val_p)) 
        tp0 = sum(model.predict_classes(X_train_p,verbose=vb0))+0.1
        #print('Train_Negative')
        #print(model.predict_classes(X_val_n)) 
        fp0 = sum(model.predict_classes(X_train_n,verbose=vb0))
        tn0 = ntotal0 - fp0
        fn0 = ptotal0 - tp0
        train_pre = str(float(tp0)/(tp0+fp0))
        train_recall = str(float(tp0)/(tp0+fn0))
        train_f1 = calf1(train_pre, train_recall)
        #print('Train_Precision='+str(float(tp0)/(tp0+fp0)))
        #print('Train_Recall='+str(float(tp0)/(tp0+fn0)))
        
        ######predicting validation
        #print('Val_Postive')
        #print(model.predict_classes(X_val_p)) 
        ptotal0 = len(X_val_p)
        ntotal0 = len(X_val_n)
        #predict
        p_predicted = model.predict_classes(X_val_p,verbose=vb0)
        #overall true positive
        tp0 = sum(p_predicted)+0.1
        #recall = tp/(total positive by gold standard)
        #print('\n')
        #print('X_train')
        #print('p_predicted='+str(len(p_predicted)))
        #print('mask_non_i='+str(len(mask_non_i)))
        recall_non_i = sum(p_predicted[mask_non_i])/float(len_non_i)
        recall_non_sub = sum(p_predicted[mask_non_sub])/float(len_non_sub)
        #print('Val_Negative')
        #print(model.predict_classes(X_val_n)) 
        fp0 = sum(model.predict_classes(X_val_n,verbose=vb0))
        tn0 = ntotal0 - fp0
        fn0 = ptotal0 - tp0
        val_pre=str(float(tp0)/(tp0+fp0))
        val_recall=str(float(tp0)/(tp0+fn0))
        val_f1 = calf1(val_pre,val_recall)
        #print('Val_Precision='+str(float(tp0)/(tp0+fp0)))
        #print('Val_Recall='+str(float(tp0)/(tp0+fn0)))
        print('\n')
        print('[iteration,train_pre,train_recall,train_f1,val_pre,val_recall,val_f1,recall_non_i,recall_non_sub]')
        print([iteration,train_pre,train_recall,train_f1,val_pre,val_recall,val_f1,recall_non_i,recall_non_sub])
        output_perf2([iteration,train_pre,train_recall,train_f1,val_pre,val_recall,val_f1,recall_non_i,recall_non_sub])
        model.save_weights(path_save+file_name0+out_name+'_weight.h5',overwrite=True)
    #save weights and performance info
    #output_perf(file_out,file_name0,iterations,training_n, train_pre,train_recall,val_pre,val_recall)
    #model.save_weights(path_save+file_name0+v1+'_weight.h5',overwrite=True)