#taking the training data 
#sequence \t binding affinity \t DR type (e.g. DRB1_0101)
#code to convert a list into cluster 
from utilities import *
import fileinput
import math
import random
import re
import subprocess

    
def format_hla(str0):
    str_out = 'HLA-DRB1*'+str0[-4:-2]+':'+str0[-2:]
    return str_out

########set path and output files
###output path
path_save = '/home/stanford/rbaltman/users/bchen45/data/HLA_pred_data/'
#folder for encoding dictionary
path_encoding = '/home/stanford/rbaltman/users/bchen45/code/python_general/encoding_dict/'
hla_dict_file = 'DRB1_pseudo_seq.dict'
###initiate the output file
file_name0 = 'netMHCIIpan_train1'
touch_file(path_save+file_name0)
file_out = open(path_save+file_name0,'a')

############import dictionary for hla -> pseudosequence
dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))

##############main###################################
###########use command line file input###############
for line0 in fileinput.input():
    [seq0,aff0,hla0] = line0.split('\t')
    if 'DRB1' in hla0:
        #only use DRB1 for this model
        hla0 = format_hla(hla0)
        hla_seq0 = dict_hla[hla0]
        write_list(file_out,[hla_seq0+seq0,aff0],'\t')
    
    
        

