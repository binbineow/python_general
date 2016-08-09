# -*- coding: utf-8 -*-
#including chemistry information
#for v13, iterations will be reduced to 20 and the best non-overfit performance will be recorded
#non-overfit is defiend as iteration > 10, F1(training)-F1(Validation)<5% (5% is working approximation)

from __future__ import print_function
from keras.models import Sequential, slice_X
from keras.layers.core import Activation, Masking, Dropout, Dense, RepeatVector
from keras.layers import recurrent
from keras.callbacks import ModelCheckpoint
#from utilities import *

# model reconstruction from JSON:
from keras.models import model_from_json
#path_save = '/scratch/users/bchen45/HLA_prediction/RNN_data/'
# how to use the model elsewhere...
#model = model_from_json(open('my_model_architecture.json').read())
#model.load_weights('my_model_weights.h5')

##import coding dictionary
#path_dict = '/scratch/users/bchen45/code/python_general/python_general/encoding_dict/'
#Blosum50_sparse.dict
#Blosum50_only.dict
#Sparse_only.dict


# Parameters for the model and dataset
#TRAINING_SIZE = len(inputs)
# Try replacing JZS1 with LSTM, GRU, or SimpleRNN
RNN = recurrent.JZS1
n_iteration = 20
HIDDEN_SIZE = 30
BATCH_SIZE = 20
LAYERS = 2
ratio_t = 1

    

#start a model
model = Sequential()
# "Encode" the input sequence using an RNN, producing an output of HIDDEN_SIZE
#model.add(Masking())
model.add(RNN(HIDDEN_SIZE, input_shape=(None, len(chars)), return_sequences=True))
for _ in xrange(LAYERS):
    model.add(RNN(HIDDEN_SIZE, return_sequences=True))
#    #model.add(Dropout(0.5))
model.add(RNN(HIDDEN_SIZE, return_sequences=False))
model.add(Dense(len(classes)))
model.add(Activation('softmax'))
model.compile(loss='categorical_crossentropy', optimizer='adam')
model1 = model

path_save = '/home/stanford/rbaltman/users/bchen45/code/test/'

from keras.utils.visualize_util import plot
plot(model, to_file=path_save+'hla_rnn_model_nodropout_layer_real2.png')

