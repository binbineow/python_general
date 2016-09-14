from sklearn.metrics import roc_auc_score
import numpy as np
import cPickle as pickle



def get_val_list(dict0,n0):
    list_out = []
    for key0,val0 in dict0.iteritems():
        if n0 == 0:
            list_out.append(max(val0[n0],val0[n0+3]))
        else:
            list_out.append(min(val0[n0],val0[n0+3]))
    return list_out

file_name_pos = 'netmhc_predict_MCLJeko.pos.dict'
file_name_neg = 'netmhc_predict_MCLJeko.neg.dict'

dict_pos = pickle.load(open(file_name_pos,'r'))
dict_neg = pickle.load(open(file_name_neg,'r'))

for n0 in range(0,3):
    print('using the first value '+str(n0)+' at '+file_name_pos)
    list_val_p =get_val_list(dict_pos,n0)
    list_val_n =get_val_list(dict_neg,n0)
    list_values = np.concatenate((list_val_p,list_val_n))
    if n0 == 0:
        list_true = np.concatenate((np.ones(len(list_val_p)),np.zeros(len(list_val_n))))
    else:
        list_true = np.concatenate((np.zeros(len(list_val_p)),np.ones(len(list_val_n))))
    auc_val = roc_auc_score(list_true, list_values)
    print(str(auc_val))
    
