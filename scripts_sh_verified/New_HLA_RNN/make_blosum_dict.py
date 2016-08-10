import fileinput
import cPickle as pickle 
import numpy as np
n_aa = 23
order0 = ''
dict_aa = dict()
for line0 in fileinput.input():
    np0 = np.zeros(n_aa*2)
    line0 = line0.rstrip()
    line0 = line0.split(' ')
    line0 = filter(None,line0)
    n_current = n_aa
    #print line0
    for n0,x0 in enumerate(line0):
        if n0 == 0:
            order0 = order0 + x0
            np0[len(order0)-1] = 1
        else:
            np0[n_current] = int(x0)
            #print n_current
            n_current += 1
    dict_aa[line0[0]] = np0
print dict_aa['P']
print order0
dict_aa['order'] = order0
pickle.dump(dict_aa,open('Blosum50_sparse.dict','w+'))
dict_chemonly = {k:v[23:] for k,v in dict_aa.iteritems()}
dict_sparseonly = {k:v[0:23] for k,v in dict_aa.iteritems()}
# for item0 in dict_chemonly['order']:
#     dict_chemonly[item0] = dict_aa[item0][23:]
#     dict_sparseonly[item0] = dict_aa[item0][0:23]
print dict_aa['P']
print dict_chemonly['P']
print dict_sparseonly['P']
pickle.dump(dict_chemonly,open('Blosum50_only.dict','w+'))
pickle.dump(dict_sparseonly,open('Sparse_only.dict','w+'))


 
    
    