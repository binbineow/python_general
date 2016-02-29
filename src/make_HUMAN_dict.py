from utilities import *
import cPickle as pickle
import fileinput
dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/UniprotID.dict'

def make_any_dict(dictName,pickle_s,header_s):
    #dict0 = pickle.load(open(dict_file,'r'))
    for line0 in fileinput.input():
        if header_s:
            header_s = False
        else:
            line0 = line0.rstrip('\t')
            dict0[line0[0][:-]] = line0[1]
    if pickle_s:
        pickle.dump(dict0,open(dictName,'w+'))
    else:
        return dict0