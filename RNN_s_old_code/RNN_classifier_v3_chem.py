# -*- coding: utf-8 -*-
from __future__ import print_function
from keras.models import Sequential, slice_X
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from sklearn.utils import shuffle
import numpy as np
import fileinput
import Levenshtein
n_ratio = 3 #Negative to positive ratio
n_split = 10 #x out of 100

# create a python dictionary, map letter (string) to a numpy vector (row/column)
# dict["A"] -> numpy vector

class CharacterTable(object):
    def __init__(self, chars, maxlen):
        self.chars = sorted(set(chars))
        self.char_indices = dict((c, i) for i, c in enumerate(self.chars))
        self.indices_char = dict((i, c) for i, c in enumerate(self.chars))
        self.maxlen = maxlen

    def encode(self, C, maxlen=None,chemical=False):
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
              X[i] = dict[c]
          return X

    def decode(self, X, calc_argmax=True):
      idxs = [i for i,x in enumerate(X) if x == 1][0]
      return self.indices_char[idxs]

class colors:
    ok = '\033[92m'
    fail = '\033[91m'
    close = '\033[0m'

inputs=[]
outputs=[]
char_set = set([' '])
class_set = set()
max_len = 0

def make_train_val(X,y,tbl):
    print(X,y)
    # Shuffle (X, y) in unison as the later parts of X will almost all be larger digits
    X_, y_ = shuffle(X, y)
    X_train, X_val, y_train, y_val = [],[],[],[]
    n_goal_p = round(len(X_)/(float(n_ratio+1))*n_split/100)
    n_goal_n = n_goal_p
    # Explicitly set apart 10% for validation data that we never train over
    for i in range(0,len(X_)):
        #print(X_[i])
        #print(len(X_val),len(X_train))
        if tbl.decode(y_[i]) == '1':
            if n_goal_p>0:
                X_val.append(X_[i])
                y_val.append(y_[i])
                n_goal_p -= 1
            else:
                X_train.append(X_[i])
                y_train.append(y_[i])
        else:
            if n_goal_n>0:
                X_val.append(X_[i])
                y_val.append(y_[i])
                n_goal_n -= 1
            else:
                X_train.append(X_[i])
                y_train.append(y_[i])
    return [np.array(X_train),np.array(X_val),np.array(y_train),np.array(y_val)]
                


for line in fileinput.input():
  in_,out_ = [x.rstrip() for x in line.split("\t")]
  max_len = max([max_len,len(in_),len(out_)])
  inputs.append(in_)
  outputs.append(out_)
  if len(out_) != 1:
    raise Exception("Output should be single characer")
  for c in in_: char_set.add(c)
  class_set.add(out_)

# Parameters for the model and dataset
TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN
RNN = recurrent.JZS1
HIDDEN_SIZE = 28
BATCH_SIZE = 5
LAYERS = 2
MAXLEN = max_len #DIGITS + 1 + DIGITS

chars = ''.join(char_set) #'0123456789+ '
classes = ''.join(class_set)
ctable = CharacterTable(chars, MAXLEN)
classtable = CharacterTable(classes, 1)
X = np.zeros((len(inputs), MAXLEN, len(chars)), dtype=np.bool)
y = np.zeros((len(inputs), len(classes)), dtype=np.bool)
for i, sentence in enumerate(inputs):
    X[i] = ctable.encode(sentence, maxlen=MAXLEN)
for i, c in enumerate(outputs):
    y[i] = classtable.encode(c, maxlen=1)#MAXLEN)

print("loaded input")

X_train,X_val,y_train,y_val = make_train_val(X,y,classtable)
#print(X_train,X_val,y_train,y_val)
#print(len(X_train),len(X_val))

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


# Train the model each generation and show predictions against the validation dataset
for iteration in range(1, 30):
    print()
    print('-' * 50)
    print('Iteration', iteration)
    model.fit(X_train, y_train, batch_size=BATCH_SIZE, nb_epoch=1, class_weight={1:1,0:1/n_ratio},validation_data=(X_val, y_val), show_accuracy=True)

