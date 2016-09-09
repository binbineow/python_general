#make a dictionary = HLA-type to list of peptides
#code to convert a list into cluster 
from utilities import *
from predict_netmhciipan import *
import math
import random
import re
import subprocess
from collections import defaultdict
import os
from flake8.formatting.default import Default

path_save = '/share/PI/rbaltman/bchen45/software/IEDB/MCL_netmhc_predict_results'
path0 = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
one_gene_path = '/share/PI/rbaltman/bchen45/software/IEDB/test0/human_proteinome_oneline.str'
'''
####iedb_path should be in the bash path file, so the program can call NetMHCIIpan directly

#convert HLA-DRB1*04:01 into DRB1_0401 format
def convert_hla(str0):
    str1 = str0.split('*')[1].split(':')
    num1 = str1[0]
    num2 = str1[1]
    return 'DRB1_'+num1+num2

#given a list of strings, return a list of their lengths in order
def get_len_list(list1):
    list0 = []
    for x in list1:
        list0.append(len(x))
    return list0

def make_one_line(list0):
    str_out = ''
    for x in list0:
        str_out = str_out + x
    return str_out

def read_netmhc_xls(file0,list_run):
    file_out = open(file0,'r')
    dict0 = dict()
    for line0 in file_out:
        line0 = line0.rstrip()
        if 'Sequence' in line0:
            line0 = line0.split('\t')
        if line0[1] in list_run:
            dict0[line0[1]] = [float(line0[3]),float(line0[4]),float(line0[5])]
    return dict0

def remove_file(file0):
    cmd_line = 'rm '+file0
    cmd0 = subprocess.Popen(cmd_line, shell=True)
    cmd0.wait()     

def run_netmhciipan(hla_type_run,list_run,len_run):
    set_run = set(list_run)
    list_run = list(set_run)
    str_long = make_one_line(list_run)
    #create a pep file
    file_name_in = hla_type_run+str(len_run)+'.pep'
    file_pep = open(file_name_in,'w+')
    file_pep.write('>temp0\n')
    file_pep.write(str_long)
    file_pep.close()    
    #run
    #cmd_line = 'netMHCIIpan -f '+file_name_in+ ' -inptype 0 -a '+ hla_type_run+ \
    # ' -length '+str(len_run)+' > '+file_name_in+'.temp' + ' -xls -xlsfile '+file_name_in+'.xls ' + \
    #'-tdir /home/stanford/rbaltman/users/bchen45/software/netMHCIIpan-3.1/tmp'
    cmd_line = 'netMHCIIpan -f '+file_name_in+ ' -inptype 0 -a '+ hla_type_run+ ' >'+file_name_in+'.temp'\
    ' -length '+str(len_run)+ ' -xls -xlsfile '+file_name_in+'.xls ' + \
    '-tdir /share/PI/rbaltman/bchen45/software/IEDB/netMHCIIpan-3.1/tmp'
    print cmd_line
    #cmd_line_list = cmd_line.split(' ')
    #print cmd_line_list
    #subprocess.call(cmd_line_list)
    cmd0 = subprocess.Popen(cmd_line,shell=True)      
    cmd0.wait() 
    #get data
    dict_out = read_netmhc_xls(file_name_in+'.xls', list_run)
    #remove temp file
    if os.path.isfile(file_name_in):
        remove_file(file_name_in)
        remove_file(file_name_in+'.temp')
        remove_file(file_name_in+'.xls')
    return dict_out

#main function
#give the hla_type in HLA_DRB1*04:01 format, and a list of sequences, return a dictionary
#The dictionary maps individual sequences -> [binding affinity, ranking score 
def predict_netmhciipan(hla_type0,list_seq):
    dict0 = dict()
    list_len0 = get_len_list(list_seq)
    set_len0 = set(list_len0)
    if '*' in hla_type0:
        hla_type_run = convert_hla(hla_type0)
    else:
        hla_type_run = hla_type0
    for len_run in set_len0:
        list_run = []
        for x in range(0,len(list_len0)):
            if list_len0[x] == len_run:
                list_run.append(list_seq[x])
        dict_run = run_netmhciipan(hla_type_run,list_run,len_run)
        print(dict_run)
        dict0.update(dict_run)
    return dict0
'''
def make_neg(pep_list):
    list_out = []
    for x in pep_list:
        neg0 = ''.join(random.sample(x,len(x)))
        rand0 = random.randint(0,len_one)
        neg1= onegenestr[rand0:rand0+len(x)]
        list_out.append(neg0)
        list_out.append(neg1)
    return list_out

#requires every elements in the list >= 9
def keep_long(list0):
    list_out = []
    for x in list0:
        if len(x)>8:
            list_out.append(x)
    return list_out

def update_dict(dict1,dict2):
    for key0, value0 in dict1.iteritems():
        dict1[key0] = dict1[key0] + dict2[key0]
    return dict1
        

def make_predict_dict(pid0,hla1,hla2):
    dict_run_pos = dict()
    dict_run_neg = dict()
    pep_list = keep_long(MCL_data[pid0]['MHC2_frag'])
    pep_list_neg = make_neg(pep_list)
    dict_run_pos = predict_netmhciipan(hla1,pep_list)
    dict_second = predict_netmhciipan(hla2,pep_list)
    dict_run_pos = update_dict(dict_run_pos,dict_second)
    dict_run_neg = predict_netmhciipan(hla1,pep_list_neg)
    dict_second = predict_netmhciipan(hla2,pep_list_neg)
    dict_run_neg = update_dict(dict_run_neg,dict_second)
    return [dict_run_pos,dict_run_neg]

#generating the random peptide sequence
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
#read in data
#patient_val = ['MCL041','MCL128','MCL019']
patient_target = ['MCL128','MCL041','MCL019']
patient_target = ['MCL019']
MCL_data = pickle.load(open(path0+'MCL_data11_18_2015v1.1.dict','r'))
#dict_hla = pickle.load(open(path_encoding+hla_dict_file,'r'))
#initiate the training set
set_train = set()
#dict_pos = defaultdict(defaultdict)
#dict_neg = defaultdict(defaultdict)
dict_pos = defaultdict()
dict_neg = Defaultdict()
#write training data into a txt file
if len(patient_target)<1:
    patient_target = MCL_data['pid']['pid']
for pid0 in patient_target:
    hla1 = MCL_data[pid0]['HLA_typing'][-1]
    hla2 = MCL_data[pid0]['HLA_typing'][-2]
    [dict_pos,dict_neg] = make_predict_dict(pid0,hla1,hla2)

#save
pickle.dump(dict_pos,open(path_save+'netmhc_predict_mcl019.pos.dict'))
pickle.dump(dict_neg,open(path_save+'netmhc_predict_mcl019.neg.dict'))   
