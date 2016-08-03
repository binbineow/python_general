#this script maps desired short peptides to a protein of interest given a cell line name
#output is a heatmap representing the most common region the peptides are recovered

path_in = '/Users/binbineow2/Documents/Machine_Learning/MCL_data/Cell_line/'
path_out = '/Users/binbineow2/Documents/Machine_Learning/MCL_data/Cell_line/'
file_seq = 'sequence.txt'


import pandas as pd
import Levenshtein
from collections import defaultdict

cell_name = ['JEKO','L128']
lev_cutoff = 0.9

def get_seq(name0):
    #get sequence info from the sequence file
    found0 = False
    for line0 in open(path_in+file_seq,'r'):
        line0 = line0.rstrip()
        if found0:
            return line0
            break
        if line0  == name0:
            found0 = True
def check_Lev(str_long, str_short):
    #output potential similar peptides with the detected peptide, similar peptide, and Levenhstein ratio
    len0 = len(str_short)
    for n0 in range(0,len(str_long)-len0+2):
        ratio0 = Levenshtein.ratio(str_long[n0:n0+len0],str_short)
        if ratio0 > lev_cutoff:
            return ([str_short, str_long[n0:n0+len0],str(ratio0)])
            break
    return 'not matched'
        
def update_dic(dict0,dict_area,n0,len0,area0):
    #update dictionary based on binary or area
    for x in range(n0,n0+len0):
        #print x
        dict0[x] += 1
        dict_area[x] += area0
    return [dict0,dict_area]
    
def print_header(file_out,seq0,note0):
    #global file_out
    if not note0 == '':
        file_out.write(note0+',')
    for x in seq0:
        file_out.write(str(x)+',')
    file_out.write('\n')

def create_dict_n(n0):
    dict0 = defaultdict()
    for x in range(0,n0):
        dict0[x] = 0
    return dict0

def dict_to_list(dict0):
    len0  = len(dict0)
    list0 = []
    for x in range(0,len0):
        list0.append(dict0[x])
    return list0

def main(name0):
    #initiate output
    file_out = open(path_out+name0+'_out.csv','w+')
    file_lev = open(path_out+name0+'_lev.csv','w+')
    #read in data
    seq0 = get_seq(name0)
    dict0 = create_dict_n(len(seq0))
    dict_area = create_dict_n(len(seq0))
    print_header(file_out, seq0, 'Protein_sequence')
    file_lev.write('Detected peptide sequence,Potential matching peptide,Similarity score')
    pep_df = pd.read_csv(path_in+name0+'_pep_chart.csv',',')
    pep_df = pep_df.fillna(0)
    pep_seq = list(pep_df['Sequence'])
    pep_area = list(pep_df['Area: F1: Sample'])
    #checking
    #print dict0
    #print pep_area
    #matching strings
    count0 = 0
    for pep0 in pep_seq:
        len0 = len(pep0)
        match0 = seq0.find(pep0)
        if match0 >= 0:
            dict0,dict_area = update_dic(dict0, dict_area, match0, len0 , float(pep_area[count0]))
        else:
            try_lev = check_Lev(seq0,pep0)
            if len(try_lev) == 3:
                print_header(file_lev,try_lev,'')  
        count0 += 1
    #output heatmap results    
    #checking
    #print dict0
    #print dict_area
    list0 = dict_to_list(dict0)
    list_area = dict_to_list(dict_area)
    print_header(file_out,list0,'Binary counting')
    print_header(file_out,list_area,'Area weighted counting')    
        
    
    #close files
    file_lev.close()
    file_out.close()

for name0 in cell_name:
    main(name0) 