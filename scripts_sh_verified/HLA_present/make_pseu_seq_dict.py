#This function makes a dictonary based on NetMHCIIpan's pseduosequence
#list but only grabbing the beta sequence (after 15AA, the alpha part)
#and format to my current HLAII format,e.g. HLA-DRB1*11:04 (from DRB1_1156).
import cPickle as pickle
path_in = '/scratch/users/bchen45/HLA_prediction/IEDB/netMHCIIpan-3.1/data/'
file_in = 'pseudosequences.dat'
path_out = '/scratch/users/bchen45/code/python_general/python_general/encoding_dict/'
file_out = 'DRB1_pseudo_seq.dict'

def format_hla(str0):
    str_out = 'HLA-DRB1*'+str0[-4:-2]+':'+str0[-2:]
    return str_out

def get_seq(str0):
    return str0[15:]
    
#main
dict0 = dict()

for line0 in open(path_in+file_in):
    line0 = line0.rstrip().split('\t')
    if 'DRB1_' in line0[0]:
        line0[0] = format_hla(line0[0])
        line0[1] = get_seq(line0[1])
        dict0[line0[0]] = line0[1]
        
#output
pickle.dump(dict0,open(path_out+file_out,'w+'))
''' 
dict1 = dict()
for key0,val0 in dict0.iteritems():
    if 'DRB1_' in key0:
        key1 = format_hla(key0)
        dict1[key1] = val0
pickle.dump(dict1,open('DRB1_34_encoding.dict','w+'))

'''  
