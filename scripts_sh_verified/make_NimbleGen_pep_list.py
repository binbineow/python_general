import fileinput
from utilities import *
target_mhc = ['HLA-DRB1*04:01','HLA-DRB1*07:01','HLA-DRB1*15:01']
path_MCL = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
path_IEDB = '/scratch/users/bchen45/HLA_prediction/IEDB/raw_data/'
path_core = '/scratch/users/bchen45/HLA_prediction/RNN_data/core/'

#master format
#in the output csv file, each peptide has the following information
#0Sequence
#1Type   
#2ID    
#3Qualitative Measurement    
#4Quantitative Measurement    
#5Predicted binding core    
#6Note (mannually added)

def discover_core():  

#1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
def get_IEDB_pep():
    list_target = []
    for line0 in fileinput.input():
        line0=line0.rstrip()
        line0=line0.split('"')
        list0 = []
        if len(line0) > 127:
            print line0[127]
            if line0[127] in target_mhc:
                list0 = [line0[1],line0[8],line[23],line[101],line[109],line[111],line[127]]
            if not list0 == []:
                list_target.append(list0)
    for list0 in list_target:
        line0 = ''
        for ele0 in list0:
            line0=line0+','+ele0
        print line0[1:]

def get_MCL_pep():
    

def get_lisa_pep():