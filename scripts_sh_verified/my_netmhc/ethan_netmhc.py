from __future__ import print_function
import numpy as np
np.random.seed(1337)  # for reproducibility
import sys
import math
from scipy.stats import pearsonr
from keras.preprocessing import sequence
from keras.models import Model
from keras.layers import merge, GRU, Dense, Dropout, Embedding, LSTM, Input, Bidirectional
from keras.layers.core import Flatten
from keras.datasets import imdb
from keras.regularizers import l2, activity_l2
import pickle

# load blossom chemical encoding dictionary
chem = pickle.load(open(sys.argv[1],"rb"))
# set '-' to all zeros
chem['-'] = np.zeros(len(chem['A']))

# encode string of peptides with chemical information
def encode_seq(seq):
  return np.array([chem[k] for k in seq]).flatten()

# generate sliding windows over core for netmhc technique
# input: string of peptides
# output: yield each netmhc vector encoding
def make_windows(seq):
  first = seq[0]
  last = seq[-1]
  for i in range(0,len(seq)-9+1):
    core = seq[i:i+9]
    num_begin = i
    num_end = len(seq)-(i+9)
    nums = np.array([num_begin, 1 - num_begin, num_end, 1 - num_end, len(seq), 1-len(seq)])
    yield np.concatenate((encode_seq(first+core+last),nums))

# load data from files on netmhc website
# first column, fixed length peptide string
# second column, variable length peptide string
# third column, binding strength
# input: file
# output: x_data matrix, y_data matrix, tracking info dictionary
def load_data(filename):
    X, y = [], []
    track = {}
    with open(filename,"r") as f:
        for i,line in enumerate(f):
            cols = line.rstrip().split()
            # vector for fixed peptide
            fixed_rep = encode_seq(cols[0])
            # track which original datapoint each window belongs to
            begin, end = len(X), len(X)
            for s in make_windows(cols[1]):
              X.append(np.concatenate((fixed_rep,s)))
              y.append(float(cols[2]))
              end += 1
            # data in the range of (begin, end) belongs to datapoint i
            track[i] = (begin, end)
    return np.array(X), np.array(y), track

# load training data file
X_train, y_train, track = load_data(sys.argv[2])
# load testing data file
X_test, y_test, track_test = load_data(sys.argv[3])

# neural network
input_layer = Input(shape=(X_train.shape[1],), name='main_input')
dense_layer = Dense(64,activation="relu")(input_layer)
output = Dense(1,activation="linear")(dense_layer)
model = Model(input=[input_layer], output=[output])
model.compile('adam', 'mse')

# select strongest binding core for each example
# input: x data matrix, y data, tracking info
# output: strongest bindings for each x datapoint, corresponding y, corresponding outputs
def select_best(x,y,track):
  X_new, y_new, outs_new = [], [], []
  outs = model.predict(x, batch_size=256)
  for k,vs in track.items():
    best_max = -float("inf")
    for j in range(vs[0],vs[1]):
      if best_max < outs[j]:
        best_max = outs[j]
        selected_x = x[j]
        selected_y = y[j]
    X_new.append(selected_x)
    y_new.append(selected_y)
    outs_new.append(best_max)
  return np.array(X_new), np.array(y_new), np.array(outs_new)

# save a model to file
def save_model(model,i):
  json_string = model.to_json()
  with open('model{}.json'.format(i), 'w+') as f:
    f.write(json_string)
  model.save_weights('model{}_weight.h5'.format(i),overwrite=True)

best_score = -float("inf")
scores = []
for i in range(0,300):
  print(i)
  # select strongest binders for each datapoint
  X_new, y_new, outs_new = select_best(X_train, y_train, track)
  # train on selected binders
  model.fit(X_new, y_new, batch_size=256, nb_epoch=1)
  # use model to predict binding scores on test data
  _, y_test_select, test_outs = select_best(X_test, y_test, track_test)
  # compute pearson on test data
  score, pval = pearsonr(y_test_select, test_outs.flatten())
  print(score)
  scores.append(score)
  with open("eval-scores.pkl","wb") as fff:
    pickle.dump(scores,fff)
  if score > best_score:
    best_score = score
    save_model(model,i)
