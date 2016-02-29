# -*- coding: utf-8 -*-
from __future__ import print_function
from keras.models import Sequential, slice_X
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from keras.callbacks import ModelCheckpoint
from sklearn.utils import shuffle
import numpy as np
import fileinput
import Levenshtein
import cPickle as pickle
n_ratio = 4 #Negative to positive ratio
n_split = 50 #x out of 100
n_iteration = 100
n_goal_n0 = 300
model_name = 'MCLX001'
n_version = 'v1'
p_Leven = 0.8

dict_chem = dict()

  

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
        

class colors:
    ok = '\033[92m'
    fail = '\033[91m'
    close = '\033[0m'

inputs=[]
outputs=[]
char_set = set([' '])
class_set = set()
max_len = 0

def make_train_val(X,y,tbl,tblseq):
    #print(X,y)
    # Shuffle (X, y) in unison as the later parts of X will almost all be larger digits
    X_, y_ = shuffle(X, y)
    X_train, X_val_n, X_val_p, y_train,y_val_n,y_val_p = [],[],[],[],[],[]
    n_goal_p = round(len(X_)/(float(n_ratio+1))*n_split/100)
    n_goal_n = n_goal_n0
    n_set =round(len(X_)/(n_ratio+1))-n_goal_p
    set_peptide = set()
    for i in range(len(X_)-1,0,-1):
        if tbl.decode(y_[i]) == '1':
            set_peptide.add(tblseq.decode(X_[i],class0=False))
            n_set -= 1
        if n_set < 1:
            break
    print(len(set_peptide))
    # Explicitly set apart 10% for validation data that we never train over
    for i in range(0,len(X_)):
        #print(X_[i])
        #print(len(X_val),len(X_train))
        if tbl.decode(y_[i]) == '1':
            if n_goal_p>0:
                include0 = True
                for test0 in set_peptide:
                    if Levenshtein.ratio(test0,tblseq.decode(X_[i],class0=False)) > p_Leven:
                        include0 = False
                        break
                if include0:
                    X_val_p.append(X_[i])
                    y_val_p.append(y_[i])
                else:
                    X_train.append(X_[i])
                    y_train.append(y_[i])                    
                n_goal_p -= 1
            else:
                X_train.append(X_[i])
                y_train.append(y_[i])
        else:
            if n_goal_n>0:
                X_val_n.append(X_[i])
                y_val_n.append(y_[i])
                n_goal_n -= 1
            else:
                X_train.append(X_[i])
                y_train.append(y_[i])
                
    print('len of X_val_n='+str(len(X_val_n)))
    print('len of X_val_p='+str(len(X_val_p)))
    #this section is to predict with Levenhstein
    
    
    
    return [np.array(X_train),np.array(X_val_p),np.array(X_val_n),np.array(y_train),np.array(y_val_p),np.array(y_val_n)]



tests = []
for line in fileinput.input():
  in_,out_ = [x.rstrip() for x in line.split("\t")]
  if out_ == '2':
      tests.append(in_)
  else:
      max_len = max([max_len,len(in_),len(out_)])
      inputs.append(in_)
      outputs.append(out_)
      if len(out_) != 1:
        raise Exception("Output should be single characer")
      for c in in_: char_set.add(c)
      class_set.add(out_)

print(class_set)

# Parameters for the model and dataset
TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN
RNN = recurrent.JZS1
HIDDEN_SIZE = 28
BATCH_SIZE = 20
LAYERS = 2
MAXLEN = max_len #DIGITS + 1 + DIGITS


chars = ''.join(char_set) #'0123456789+ '
classes = ''.join(class_set)
ctable = CharacterTable(chars, MAXLEN)
classtable = CharacterTable(classes, 1)
X = np.zeros((len(inputs), MAXLEN, len(chars)), dtype=np.bool)
#X_test = np.zeros((len(tests), MAXLEN, len(chars)), dtype=np.bool)
y = np.zeros((len(inputs), len(classes)), dtype=np.bool)
for i, sentence in enumerate(inputs):
    X[i] = ctable.encode(sentence, maxlen=MAXLEN)
for i, c in enumerate(outputs):
    y[i] = classtable.encode(c, maxlen=1)#MAXLEN)
#for i, c in enumerate(tests):
#    X_test[i] = ctable.encode(c, maxlen=MAXLEN)#MAXLEN)

print("loaded input")

[X_train,X_val_p,X_val_n,y_train,y_val_p,y_val_n] = make_train_val(X,y,classtable,ctable)
#print(X_train,X_val,y_train,y_val)
#print(X_val_p)
#print(X_val_n)
#X_val = [X_val_n,X_val_p]
#X_val = np.array(X_val)
#y_val = [y_val_n,y_val_p]
#y_val = np.array(y_val)
#print(X_val)
X_val = np.concatenate((X_val_n,X_val_p))
y_val = np.concatenate((y_val_n,y_val_p))
print(len(X_train),len(X_val))

model = Sequential()
# "Encode" the input sequence using an RNN, producing an output of HIDDEN_SIZE
#model.add(Masking())
model.add(RNN(HIDDEN_SIZE, input_shape=(None, len(chars)), return_sequences=True))
#for _ in xrange(LAYERS-1):
#    model.add(RNN(HIDDEN_SIZE, HIDDEN_SIZE, return_sequences=True))
#    #model.add(Dropout(0.5))
model.add(RNN(HIDDEN_SIZE, return_sequences=False))
model.add(Dense(len(classes)))
model.add(Activation('softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adam')
#Create checkpoint
#checkpointer = ModelCheckpoint(filepath=model_name+'.weight', verbose=1, save_best_only=True)
# Train the model each generation and show predictions against the validation dataset
for iteration in range(1, n_iteration):
    print()
    print('-' * 50)
    print('Iteration', iteration)
    #to save weight callbacks=[checkpointer]
    model.fit(X_train, y_train, batch_size=BATCH_SIZE, nb_epoch=1, class_weight={1:1,0:1.0/n_ratio},validation_data=(X_val, y_val),show_accuracy=True)
    #pickle(model,open(model_name+'.model'+n_version,'w+'))
    print('Postive')
    print(model.predict_classes(X_val_p)) 
    ptotal0 = len(X_val_p)
    ntotal0 = len(X_val_n)
    tp0 = sum(model.predict_classes(X_val_p))
    print('Negative')
    print(model.predict_classes(X_val_n)) 
    fp0 = sum(model.predict_classes(X_val_n))
    tn0 = ntotal0 - fp0
    fn0 = ptotal0 - tp0
    print('Precision='+str(float(tp0)/(tp0+fp0)))
    print('Recall='+str(float(tp0)/(tp0+fn0)))
    #for testing on 049
    #tp0 = sum(model.predict_classes(X_test))
    #print('Recall on Patient MCL049='+str(float(tp0)/(len(X_test))))