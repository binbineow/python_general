import fileinput
from collections import defaultdict
import cPickle as pickle

#>sp|Q6GZX2|003R_FRG3G Uncharacterized protein 3R OS=Frog virus 3 (isolate Goorha) GN=FV3-003R PE=3 SV=1
def check_letter(str0):
    pass0 = True
    for item0 in 'BJOUXZ':
        if item0 in str0:
            pass0 = False
            break
    #check if BJOUXZ is in the string
    return pass0

def make_uniprot_dict(file_raw,file_id_dict):
    #this function takes in raw uniprot file and a dictionary converting 
    #uniprot ids to gene ids
    dict_id = pickle.load(open(file_id_dict,'r'))
    dict_out = defaultdict(list)
    seq_run = ''
    read_seq = False
    for line0 in open(file_raw):
        line0 = line0.rstrip()
        #first part read in line0 sequence info if in read in mode
        #write dictionary when the program hits anotehr gene entry
        if read_seq and not '>' in line0:
            seq_run = seq_run + line0
        elif '>' in line0 and len(seq_run)>1:
            if name_run == '' or not check_letter(seq_run):
                print seq_run
                read_seq = False
                seq_run = ''
                name2 = ''
                name_run = ''
            else:
                dict_out[name_run].append(seq_run)
                if len(name2)>1:
                    dict_out[name2].append(seq_run)
                read_seq = False
                seq_run = ''
                name2 = ''
                name_run = ''
            
        if 'OS=Homo sapiens' in line0 and len(line0.rstrip().split(' ')[0].split('|'))>2:
            name_run = line0.rstrip().split(' ')[0].split('|')[1]
            if name_run in dict_id:
                name2 = dict_id[name_run]
            else:
                name2 = ''
            read_seq = True
    #check if any sequence to be written into the dictionary
    if len(seq_run)>1 and check_letter(seq_run):
        dict_out[name_run].append(seq_run)
        if len(name2)>1:
            dict_out[name2].append(seq_run)
    
    return dict_out        
            
#put the first string of each entry in a default dict together and make a long string     
#separted by sep0 between each string
def one_string_dict(dict0,sep0):
    set0 = set()
    for _, value0 in dict0.iteritems():
        set0.add(value0[0])
    line_out = ''
    for x in set0:
        line_out = line_out+sep0+x
    return line_out

#stick a list of strings together separated by sep0
def one_string_list(list0,sep0):
    str0 = ''
    for x in list0:
        str0 = str0+sep0+x
    return str0

def dumb(): return defaultdict(list)
    
def count_exclude(pep0,exclude0):
    import re
    num_out = 0
    for x in exclude0:
        num_out += len(re.findall('N'+x+'S',pep0))
        num_out += len(re.findall('N'+x+'T',pep0))
    return num_out
        
###this script gather all strings from MHC class and patient ids interested
##into a string separated by sep0
def get_mhc_frag(file0,class_name,patient_ids):
    import cPickle as pickle
    list_out = []
    dict0 = pickle.load(open(file0,'r'))
    if patient_ids == ['all']:
        patient_ids = dict0['pid']['pid']
    for id0 in patient_ids:
        list_out.extend(dict0[id0][class_name+'_frag'])
    return one_string_list(list_out,'XXXXX')
    
###this script gather all strings from MHC class and patient ids interested
##into a string separated by sep0
def get_mhc_frag_cal(file0,class_name,patient_ids,exclude0):
    list_out = []
    total_aa = 0
    motif_aa = 0
    dict0 = pickle.load(open(file0,'r'))
    if patient_ids == ['all']:
        patient_ids = dict0['pid']['pid']
    for id0 in patient_ids:
        list_out.extend(dict0[id0][class_name+'_frag'])
    for pep0 in list_out:
        motif_run = re.findall('N.T',pep0)+re.findall('N.S',pep0) 
        motif_p_run = count_exclude(pep0,exclude0)
        motif_num = len(motif_run) - motif_p_run
        motif_aa += motif_num
        total_aa += len(pep0)
    return [float(motif_aa)/total_aa]   