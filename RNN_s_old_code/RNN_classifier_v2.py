# -*- coding: utf-8 -*-
from __future__ import print_function
from keras.models import Sequential, slice_X
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from sklearn.utils import shuffle
import numpy as np
import fileinput


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

    def decode(self, X, calc_argmax=True):
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


print(X,y)

print("loaded input")


# Shuffle (X, y) in unison as the later parts of X will almost all be larger digits
X, y = shuffle(X, y)
# Explicitly set apart 10% for validation data that we never train over
split_at = len(X) - len(X) / 10
(X_train, X_val) = (slice_X(X, 0, split_at), slice_X(X, split_at))
(y_train, y_val) = (y[:split_at], y[split_at:])

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
    model.fit(X_train, y_train, batch_size=BATCH_SIZE, nb_epoch=1, class_weight={1:0.3,0:1},validation_data=(X_val, y_val), show_accuracy=True)

