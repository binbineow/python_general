#this dictionary will take in a string of nonredundant letters
#and sparse-encoodes each character and return a library

#string for MCL peptides: ACDEFGHIKLMNPQRSTVWXY
#path_file0 = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/aa_21_sparse_encoding.dict

def make_char_one_hot(str0,path_file0):
    import cPickle as pickle 
    import numpy as np
    dict0 = dict()
    len0 = len(str0)
    for n0 in range(0,len0):
        vec0 = np.zeros(len0)
        vec0[n0] = 1
        dict0[str0[n0]] = vec0
    pickle.dump(dict0,open(path_file0,'w+'))    
    
    
    
     
    
    