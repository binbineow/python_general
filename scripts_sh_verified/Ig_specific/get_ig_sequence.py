#this script reads in MCL dictionary and return
#a dictionary of patients to two sets of peptides (constant region and variable region)
#dict_mhc_v_h[pid] -> list of mhc2 peptides
#dict_mhc_v_l[pid] -> list of mhc2 peptides
#a list of variable region peptides (all patients)
#for variable region pepetides returning with a list with the same lenght indicating the pep length
#a set peptides not in IG
#print out unique matched IgV heavy chain and light chain recovered from the patient
#format
#patient id, IgV heavy chian number, IgV light chain

#importing
from utilities import *
import math
import random
import re
import subprocess



#paths and files
#folder for the main peptide data
path0 = '/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/'
mcl_file = 'MCL_data11_18_2015v1.1.dict'
file_v_region = 'MCL_IG_pep_seqs.csv'
file_c_region = 'ig_constant_region.txt'
file_pid = 'target_patient_test.txt'

#helping Function
def get_list_from_file(file_name0):
    file0 = open(file_name0,'r')
    list_out = []
    for x in file0:
        x = x.rstrip()
        list_out.append(x)
    return list_out

#read in constant region and return two lists of gene name -> sequence
def get_c_region(file_name0):
    list_name = []
    list_seq = []
    for line0 in open(file_name0):
        line0 = line0.rstrip()
        if '#' in line0:
            list_name.append(line0[1:])
        else:
            list_seq.append(line0)
    return [list_name,list_seq]

#read in variable region files and make the following two items
#1 two dictionaries h and l [pid] -> sequence
#2 a list of all sequences (h & l)
def get_v_region(file_name0):
    list_v = []
    dict_h = dict()
    dict_l = dict()
    for line0 in open(file_name0,'rU'):
        line0 = line0.rstrip()
        if not 'Patient' in line0:
            line0 = line0.split(',')
            pid0 = line0[0]
            h0 = line0[1]
            l0 = line0[2]
            list_v.append(h0)
            list_v.append(l0)
            dict_h[pid0] = h0
            dict_l[pid0] = l0
    return [list_v,dict_h,dict_l]

def is_in_list(str0,list0):
    out0 = False
    for x in list0:
        if str0 in x:
            out0 = True
            break
    return out0

def get_ig_nonig(pid0,mhc2_frag,mhc2_gene,con_seq_list,list_var,v_h0,v_l0):
    list_nonig = []
    list_nonig_gene = []
    list_l = []
    list_h = []
    list_c = []
    set_any = set()
    num_h = 0
    num_l = 0
    for n0 in range(0,len(mhc2_frag)):
        frag0 = mhc2_frag[n0]
        gene0 = mhc2_gene[n0]
        if not frag0 in set_any:
            set_any.add(frag0)
            is_con = False
            is_ig = False
            for con0 in con_seq_list:
                if frag0 in con0:
                    list_l.append(frag0)
                    is_con = True
                    is_ig = True
                    break
            if not is_con:
                if frag0 in v_h0:
                    list_h.append(frag0)
                    num_h +=1
                    is_ig = True
                if frag0 in v_l0:
                    list_l.append(frag0)
                    num_l +=1
                    is_ig = True
            if not is_ig:
                list_nonig.append(frag0)
                list_nonig_gene.append(gene0)
    #print out results
    print(pid0+','+str(num_h)+','+str(num_l))
    return [list_c,list_h,list_l,list_nonig,list_nonig_gene]

def get_val_from_dict(dict0,make_set=True):
    list_val = []
    for _,val0 in dict0.iteritems():
        list_val0.append(val0)
    if make_set:
        list_val = list(set(list_val))
    return list_val
            
def get_nonredun_frag(long_list,short_list,short_list_gene):
    list_out = []
    list_out_gene = []
    for n0 in range(0,len(short_list)):
        if not short_list[n0] in long_list:
            list_out.append(short_list[n0])
            list_out_gene.append(short_list_gene[n0])
    return [list_out,list_out_gene]

# get dictionaries and lists of variable regions and non-Ig peptides
def get_mhc_seq(dict_mcl,pid_list,con_seq_list,list_var,dict_h,dict_l):
    list_allnonig = []
    list_allnonig_gene = []
    list_all_h = []
    list_all_l = []
    list_all_c = []
    #dict_mhc_v = dict()
    #dict_mhc_c = dict()
    #dict_mhc_nonig = dict()
    for pid0 in pid_list:
        mhc2_frag = dict_mcl[pid0]['MHC2_frag']
        mhc2_gene = dict_mcl[pid0]['MHC2_gene']
        [list_c,list_h,list_l,list_nonig,list_nonig_gene] = get_ig_nonig(pid0,mhc2_frag,mhc2_gene,con_seq_list,list_var,dict_h[pid0],dict_l[pid0])
        dict_mcl[pid0]['Variable_h'] = list_h
        dict_mcl[pid0]['Variable_l'] = list_l
        dict_mcl[pid0]['Constant'] = list_c
        dict_mcl[pid0]['Non_ig'] = list_nonig
        dict_mcl[pid0]['Non_ig_gene'] = list_nonig_gene
        dict_mcl[pid0]['Variable_h_seq'] = dict_h[pid0]
        dict_mcl[pid0]['Variable_l_seq'] = dict_l[pid0]
        list_all_h.extend(list_h)
        list_all_l.extend(list_l)
        list_all_c.extend(list_c)
        [list_add,list_add_gene] = get_nonredun_frag(list_allnonig, list_nonig, list_nonig_gene)
        list_allnonig.extend(list_add)
        list_allnonig_gene.extend(list_add_gene)
        
    return [list_allnonig,list_allnonig_gene,list_all_h,list_all_l,list_all_c]
        
        
#main Function
def main(path0,mcl_file,file_pid,file_v_region,file_c_region):
    #load input

    pid_list = get_list_from_file(path0+file_pid)
    [con_gene_list,con_seq_list] = get_c_region(path0+file_c_region)
    #print('DNA Sequence data loaded')
    [list_var,dict_h,dict_l] = get_v_region(path0+file_v_region)
    for n0 in range(0,len(con_gene_list)):
        dict_mcl['Constant'][con_gene_list[n0]]=con_seq_list[n0]
    print('DNA Sequence data loaded')
    #print(list_var)
    #print(dict_h)
    #get two pid -> mhc2 peptide dictionaries, two lists of variable regions, 
    #a list of constant region peptides,
    #a list non-Ig peptides 
    [list_allnonig,list_allnonig_gene,list_all_h,list_all_l,list_all_c] = get_mhc_seq(dict_mcl,pid_list,con_seq_list,list_var,dict_h,dict_l)
    #output
    return [list_allnonig,list_allnonig_gene,list_all_h,list_all_l,list_all_c]

def save_files(list0,list_names):
    for n0 in range(0,len(list0)):
        pickle.dump(list0[n0],open(list_names[n0],'w+'))
    
    
dict_mcl = pickle.load(open(path0+mcl_file))
[list_allnonig,list_allnonig_gene,list_all_h,list_all_l,list_all_c] = main(path0,file_pid,file_v_region,file_c_region)
list_files = ['MCL_all_nonIg.list','MCL_all_nonIG_gene.list','MCL_all_V_heavy.list','MCL_all_V_light.list','MCL_all_IG_constant.list']
list0 = [list_allnonig,list_allnonig_gene,list_all_h,list_all_l,list_all_c]
save_files(list0,list_files)
#save the new MCL data
pickle.dump(dict_mcl,open('MCL_data11_18_2015_10_9_2016v1.1.dict','w+'))

