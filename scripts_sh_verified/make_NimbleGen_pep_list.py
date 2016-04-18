import fileinput
from utilities import *
import random
import re
target_mhc = ['HLA-DRB1*04:01','HLA-DRB1*07:01','HLA-DRB1*15:01']
path_MCL = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
path_IEDB = '/scratch/users/bchen45/HLA_prediction/IEDB/raw_data/'
path_core = '/scratch/users/bchen45/HLA_prediction/RNN_data/core/'
path_Lisa = '/scratch/users/bchen45/code/python_general/python_general/small_helping_files/'
#########loading core dictionary
list_core = pickle.load(open(path_core+'drb1_MCL_IEDB_core.list'))

#########loading MCL peptides
dict_MCL = pickle.load(open(path_MCL+'MCL_pep_by_DRB1.dict'))




def discover_core(pep0):
    found0 = False
    for x in list_core:
        if x in pep0:
            return x
            found0 = True
    if not found0:
        return pep0

def print_result(file0,line0):
    for element0 in line0:
        file0.write(element+',')
    if (not line0[-1] ==  'N/A') and float(line0[-2])>16 and len(line0[-1])<=16:
        file0.write('Try core sequence since the original sequence is > 16AAs')
        
#1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
def get_IEDB_pep(hla0,file0):
    list_target = []
    for line0 in open(path_IEDB+'mhc_ligand_full.csv','r'):
        line0=line0.rstrip()
        line0=line0.split('"')
        list0 = []
        if len(line0) > 127:
            if 'HLA-DRB1' in line0[127] and not '/' in line0 [127]:
                line0[127] = line0[127][0:len('HLA-DRB1*08:02')+1]
                #print line0[127]
                if line0[127] == hla0:
                    list0=[line0[23],'IEDB',line0[1],line0[109],line0[111],line0[107],str(len(line0[23])),discover_core(line0[23])]
                    print_result(file0,list0)
    

#master format
#in the output csv file, each peptide has the following information
#0Sequence
#1Type   (assigned)
#2ID    (IEDB 1)
#3Qualitative Measurement   (IEDB-109) 
#4Quantitative Measurement    (IEDB-111)
#5Measurement Unit (IEDB-107)
#6peptide length
#7Predicted binding core    (using my list)
#8Note (mannually added)


def get_MCL_pep(file0,hla0):
    list_pep = dict_MCL(hla0)
    if len(list_pep) < 10:
        print('too few peptides in '+hla0)
    for pep0 in list_pep:
        list0 = [pep0,'MCL_patient','MCL_patient','Patient_positive','N/A','N/A',str(len(pep0)),discover_core(pep0)]
        print_result(file0,list0)
        pep0_original = pep0
        pep0 = ''.join(random.sample(pep0,len(pep0)))
        list0 = [pep0,'MCL_shuffled',pep0_original,'Negative','N/A','N/A',str(len(pep0)),'N/A']

def get_lisa_pep(file0,hla0):
    #example line
    #DR4_pp65CMVpep2,FCEDVPSGKLFMHVTLGSDV,DRB1*04:01,y,
    #my hla0 input:HLA-DRB1*07:01'
    for line0  in open(path_Lisa+'Lisa_assay_results.csv','r'):
        line0 = line0.rstrip().slipt(',')
        if not line0[0] == 'peptideName' and line0[2].split('DRB1')[1] == hla0.split('DRB1')[1]:
            if line0[3] == 'y':
                line0[3] = 'Positive'
            if line0[3] == 'n':
                line0[3] = "Negative"
            if line0[3] == 'weak':
                line0[3] = 'Positive-Low'
            list0 = [line0[1],'Tested',line0[0],line0[3],'N/A','N/A',str(len(line0[1])),discover_core(line0[1])]
            print_result(file0,list0)
            line0[1] = ''.join(random.sample(line0[1],len(line0[1])))
            list0 = [line0[1],'Tested_shuffled',line0[0]+'_shuffled','Negative','N/A','N/A',str(len(line0[1])),'N/A']
            
    
    

#open the output file and run grabbing processings
for hla0 in target_mhc:
    hla0_name = re.sub(r'[^\w]', '', hla0)
    file0 = open(path_MCL+hla0_name+'peplist_NimblGen.csv','w+')
    file0.write('Sequence,Type,ID,Qualitative,Measurement,Quantitative measurement,Measurement unit,Peptide length,Predicted binding core,Note')
    get_IEDB_pep(hla0,file0)
    get_MCL_pep(hla0,file0)
    get_lisa_pep(hla0,file0)
    file0.close()
    
    
    