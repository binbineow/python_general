'''
Created on Aug 25, 2016

@author: binbineow2
'''

import cPickle as pickle
import fileinput
import numpy as np


protvec_dict = dict()

def make_np_array(list0):
    list1 = []
    for x in list0:
        list1.append(float(x))
    return np.array(list1)
        

for line0 in fileinput.input():
    line0 = line0.rstrip()
    line0 = line0[1:-1]
    line0 = line0.split('\t')
    protvec_dict[line0[0]] = make_np_array(line0[1:])

pickle.dump(protvec_dict,open('protvec.dict','w+'))
    

