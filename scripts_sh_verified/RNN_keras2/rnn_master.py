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
print(keras.__version__)