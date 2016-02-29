# -*- coding: utf-8 -*-
from __future__ import print_function
from keras.models import Sequential, slice_X
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from keras.callbacks import ModelCheckpoint
from utilities import *

# model reconstruction from JSON:
from keras.models import model_from_json
path_save = '/scratch/users/bchen45/HLA_prediction/RNN_data/'
# how to use the model elsewhere...
#model = model_from_json(open('my_model_architecture.json').read())
#model.load_weights('my_model_weights.h5')



class CharacterTable(object):
    def __init__(self, chars, maxlen):
        self.chars = sorted(set(chars))
        self.char_indices = dict((c, i) for i, c in enumerate(self.chars))
        self.indices_char = dict((i, c) for i, c in enumerate(self.chars))
        self.maxlen = maxlen

    def encode(self, C, maxlen=None):
        maxlen = maxlen if maxlen else self.maxlen
        if maxlen == 1:
            X = np.zeros(len(self.chars))
            X[self.char_indices[C]] = 1
            return X
        else:
            X = np.zeros((maxlen, len(self.chars)))
            for i, c in enumerate(list(C)):
                if i >= maxlen: break
                X[i, self.char_indices[c]] = 1
        return X

    def decode(self, X, calc_argmax=True, class0=True):
        if class0:            
            idxs = [i for i,x in enumerate(X) if x == 1][0]
            return self.indices_char[idxs]
        else:
            if calc_argmax:
                X = X.argmax(axis=-1)
            return ''.join(self.indices_char[x] for x in X)

def encoding(matrix0,input0, ctable0,len0):
    for i, sentence in enumerate(input0):
        matrix0[i] = ctable0.encode(sentence, maxlen=len0)
    return matrix0

def output_perf(file_out, file_name0, iteraions,train_pre,train_recall,val_pre,val_recall):
    file_out.write(file_name+'_'+'iterations'+'\t')
    for x0 in iterations:
        file_out.write(x0+'\t')
    file_out.write('n')
    file_out.write(file_name+'_'+'Training_precision'+'\t')
    for x0 in train_pre:
        file_out.write(x0+'\t')
    file_out.write('n')
    file_out.write(file_name+'_'+'Training_recall'+'\t')
    for x0 in train_recall:
        file_out.write(x0+'\t')
    file_out.write('n')
    file_out.write(file_name+'_'+'Validation_precision'+'\t')
    for x0 in val_pre:
        file_out.write(x0+'\t')
    file_out.write('n')
    file_out.write(file_name+'_'+'Validation_recall'+'\t')
    for x0 in val_recall:
        file_out.write(x0+'\t')
    file_out.write('n')
    file_out.close()

for file_name0 in fileinput.input():
    
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
    
    
    for line in fileinput.input(path_save+file_name0):
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
            
            for c in in_: char_set.add(c)
            class_set.add(out_)
    file_name0 = file_name0.split('.')[0]      
    # Parameters for the model and dataset
    #TRAINING_SIZE = len(inputs)
    # Try replacing JZS1 with LSTM, GRU, or SimpleRNN
    RNN = recurrent.JZS1
    n_iteration = 30
    HIDDEN_SIZE = 28
    BATCH_SIZE = 20
    LAYERS = 2
    ratio_t = 1
    MAXLEN = max_len #DIGITS + 1 + DIGITS
         
    #creating encoding table
    print(class_set)
    chars = ''.join(char_set) #'0123456789+ '
    classes = ''.join(class_set)
    ctable = CharacterTable(chars, MAXLEN)
    classtable = CharacterTable(classes, 1)
    
    
    
    
    #create training or validation matrix
    X_train_m = np.zeros((len(X_train), MAXLEN, len(chars)), dtype=np.bool)
    X_val_p_m = np.zeros((len(X_val_p), MAXLEN, len(chars)), dtype=np.bool)
    X_val_n_m = np.zeros((len(X_val_n), MAXLEN, len(chars)), dtype=np.bool)
    X_train_p_m = np.zeros((len(X_train_p), MAXLEN, len(chars)), dtype=np.bool)
    X_train_n_m = np.zeros((len(X_train_n), MAXLEN, len(chars)), dtype=np.bool)
    y_train_m = np.zeros((len(y_train), len(classes)), dtype=np.bool)
    y_val_p_m = np.zeros((len(y_val_p), len(classes)), dtype=np.bool)
    y_val_n_m = np.zeros((len(y_val_n), len(classes)), dtype=np.bool)
    
    X_train = encoding(X_train_m, X_train,ctable,MAXLEN)
    X_train_p = encoding(X_train_p_m, X_train_p,ctable,MAXLEN)
    X_train_n = encoding(X_train_n_m, X_train_n,ctable,MAXLEN)
    X_val_p = encoding(X_val_p_m, X_val_p,ctable,MAXLEN)
    X_val_n = encoding(X_val_n_m, X_val_n,ctable,MAXLEN)
    y_train = encoding(y_train_m, y_train,classtable,1)
    y_val_p = encoding(y_val_p_m, y_val_p,classtable,1)
    y_val_n = encoding(y_val_n_m, y_val_n,classtable,1)
    
    X_val = np.concatenate((X_val_n,X_val_p))
    y_val = np.concatenate((y_val_n,y_val_p))
    print(len(X_train),len(X_val))
    print("loaded input")
    
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
    model.compile(loss='categorical_crossentropy', optimizer='adam')
    #save the model
    json_string = model.to_json()
    open(path_save+file_name0+'_model.json', 'w').write(json_string)
    
    #Create checkpoint
    #checkpointer = ModelCheckpoint(filepath=model_name+'.weight', verbose=1, save_best_only=True)
    # Train the model each generation and show predictions against the validation dataset
    file_out = open(path_save+'model_performance_simplev1.csv')
    iterations = []
    train_pre = []
    train_recall = []
    val_pre = []
    val_recall = []
    for iteration in range(1, n_iteration):
        iterations.append(iteration)
        print()
        print('-' * 50)
        print('Iteration', iteration)
        #to save weight callbacks=[checkpointer]
        model.fit(X_train, y_train, batch_size=BATCH_SIZE, nb_epoch=1, class_weight={1:1,0:1.0/ratio_t/2},validation_data=(X_val, y_val),show_accuracy=True)
        ptotal0 = len(X_train_p)
        ntotal0 = len(X_train_n)
        #print('Train_Postive')
        #print(model.predict_classes(X_val_p)) 
        tp0 = sum(model.predict_classes(X_train_p))
        #print('Train_Negative')
        #print(model.predict_classes(X_val_n)) 
        fp0 = sum(model.predict_classes(X_train_n))
        tn0 = ntotal0 - fp0
        fn0 = ptotal0 - tp0
        train_pre.append(str(float(tp0)/(tp0+fp0)))
        train_recall.append(str(float(tp0)/(tp0+fn0)))
        print('Train_Precision='+str(float(tp0)/(tp0+fp0)))
        print('Train_Recall='+str(float(tp0)/(tp0+fn0)))
        
        ######predicting validation
        #print('Val_Postive')
        #print(model.predict_classes(X_val_p)) 
        ptotal0 = len(X_val_p)
        ntotal0 = len(X_val_n)
        tp0 = sum(model.predict_classes(X_val_p))
        #print('Val_Negative')
        #print(model.predict_classes(X_val_n)) 
        fp0 = sum(model.predict_classes(X_val_n))
        tn0 = ntotal0 - fp0
        fn0 = ptotal0 - tp0
        val_pre.append(str(float(tp0)/(tp0+fp0)))
        val_recall.append(str(float(tp0)/(tp0+fn0)))
        print('Train_Precision='+str(float(tp0)/(tp0+fp0)))
        print('Train_Recall='+str(float(tp0)/(tp0+fn0)))
    #save weights and performance info
    output_perf(file_out,file_name0,iteraions,train_pre,train_recall,val_pre,val_recall)
    model.save_weights(path_save+file_name0+'_weight.h5')