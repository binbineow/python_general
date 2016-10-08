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


#################open the output file and run grabbing processings
# I suggest to add the following columns for annotation: 
# DR [DR4, DR7, DR15]; X
# Data_source [IEDB, MCL, tested];  X
# Qualitative_measurement [Positive, Negative, Shuffled]; X
# Array_type [Substitution, Single], X
# Peptide_type [Full, Core] and what else that you think is important. X
# Information from these columns will be combined in a single column that is
#  included in the data output file, so our standard data analysis scripts can 
#  be used without much modification.
####master output file



#This version will only include peptides or their core shorter than 13AA
len_cut_off = 12
#Label repeat peptides with a repeat
#Label peptides length > 12 but core < 13 with Use_core
pep_set = set()
core_set = set()
core_counter = Counter()
#collect intersection of MCL peptides
non_core_counter = Counter()
IEDB_set = set()



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
        file0.write(element0+',')
#     if (not line0[-1] ==  'N/A') and float(line0[-2])>16 and len(line0[-1])<=16:
#         file0.write('Try core sequence since the original sequence is > 16AAs')
#     elif (not line0[-1] ==  'N/A') and float(line0[-2])>16:
#         file0.write('Peptide sequence and predicted core seuqnces are both > 16AAs')
    file0.write('\n')
        
#1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
def get_IEDB_pep(hla0,file0):
    global IEDB_set
    global core_counter
    global non_core_counter
    list_target = []
    set0 = set()
    for line0 in open(path_IEDB+'mhc_ligand_full.csv','r'):
        line0=line0.rstrip()
        line0=line0.split('"')
        list0 = []
        if len(line0) > 127:
            if not '(' in line0[23]:
                if 'HLA-DRB1' in line0[127] and not '/' in line0 [127]:
                    line0[127] = line0[127][0:len('HLA-DRB1*08:02')+1]
                    #print line0[127]
                    pep0_core= discover_core(line0[23])
                    pep0 = line0[23]
                    if line0[127] == hla0:                        
                        if len(line0[23]) <= len_cut_off:
                            list0=[line0[23],'IEDB',line0[1],line0[109],line0[103],line0[111],line0[105],str(len(line0[23])),pep0_core,test_repeat(line0[23])]
                            print_result(file0,list0)
                            pep_set.add(pep0)
                            core_set.add(pep0_core)
                            set0.add(pep0)
                            set0.add(pep0_core)
                            if  not (('Negative' in line0[109]) or ('Low' in line0[109])):
                                core_counter[pep0_core] +=1
                                if len(pep0) > len(pep0_core):
                                    non_core_counter[pep0] +=1

                        elif len(pep0_core) <= len_cut_off:
                            core_set.add(pep0_core)
                            list0=[line0[23],'IEDB',line0[1],line0[109],line0[103],line0[111],line0[105],str(len(line0[23])),pep0_core,'Use_core']
                            print_result(file0,list0)
                            set0.add(pep0_core)
                            if  not (('Negative' in line0[109]) or ('Low' in line0[109])):
                                core_counter[pep0_core] +=1
    if IEDB_set == set():
        IEDB_set = set0
    else:
        IEDB_set = IEDB_set.intersection(set0)                            

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

def test_repeat(str0):
    if str0 in pep_set:
        return 'Repeat'
    else:
        return ''
        

def get_MCL_pep(hla0,file0):
    list_pep = dict_MCL[hla0]
    #build a set of short peptides for this allele 
    set0 = set()
    #print list_pep
    if len(list_pep) < 10:
        print('too few peptides in '+hla0)
    for pep0 in list_pep:
        pep0_core = discover_core(pep0)
        if len(pep0) <= len_cut_off:
            uniqe0 = test_repeat(pep0)
            list0 = [pep0,'MCL_patient','MCL_patient','Patient_positive','Mass spect','N/A','N/A',str(len(pep0)),pep0_core,uniqe0]
            print_result(file0,list0)
            pep_set.add(pep0)
            core_set.add(pep0_core)
            core_counter[pep0_core] +=1
            if len(pep0) > len(pep0_core):
                non_core_counter[pep0] +=1
            if uniqe0 == '':
                pep0_original = pep0
                pep0 = ''.join(random.sample(pep0,len(pep0)))
                list0 = [pep0,'MCL_shuffled',pep0_original,'Negative','sequence shuffled','N/A','N/A',str(len(pep0)),'N/A']
                print_result(file0,list0)
                pep_set.add(pep0)
        elif len(pep0_core) <= len_cut_off:
            list0 = [pep0,'MCL_patient','MCL_patient','Patient_positive','Mass spect','N/A','N/A',str(len(pep0)),pep0_core,'Use_core']
            print_result(file0,list0)
            core_set.add(pep0_core)
            core_counter[pep0_core] +=1
            


def get_lisa_pep(hla0,file0):
    #example line
    #DR4_pp65CMVpep2,FCEDVPSGKLFMHVTLGSDV,DRB1*04:01,y,
    #my hla0 input:HLA-DRB1*07:01'
    for line0  in open(path_Lisa+'Lisa_pep_testedv2.csv','r'):
        line0 = line0.rstrip().split(',')
        #print line0
        #print(line0[2].split('DRB1')[1])
        #print(hla0.split('DRB1')[1])
        if not line0[0] == 'peptideName' and line0[2].split('DRB1')[1] == hla0.split('DRB1')[1]:
            if line0[3] == 'y':
                line0[3] = 'Positive'
            if line0[3] == 'n':
                line0[3] = "Negative"
            if line0[3] == 'weak':
                line0[3] = 'Positive-Low'
            if len(line0[1]) <= len_cut_off:
                pep0 = line0[1]
                pep0_core = discover_core(pep0)
                list0 = [line0[1],'Tested',line0[0],line0[3],'Tetramer binding','N/A','N/A',str(len(line0[1])),pep0_core,test_repeat(line0[1])]
                print_result(file0,list0)
                line0[1] = ''.join(random.sample(line0[1],len(line0[1])))
                list0 = [line0[1],'Tested_shuffled',line0[0]+'_shuffled','Negative','sequence shuffled','N/A','N/A',str(len(line0[1])),'N/A']
                print_result(file0,list0)
                pep_set.add(pep0)
                core_set.add(pep0_core)
    
    
file0 = open(path_MCL+hla0_name+'peplist_NimblGenv2.csv','w+')
file0.write('Sequence,Data source,ID or unshuffled sequence,Qualitative Measurement,Assay type,Quantitative measurement,Measurement unit,Peptide length,Predicted binding core,MHC_type,Array_type,Peptide_type,Note\n')
    

for hla0 in target_mhc:
    hla0_name = re.sub(r'[^\w]', '', hla0)
    get_IEDB_pep(hla0,file0)
    get_MCL_pep(hla0,file0)
    get_lisa_pep(hla0,file0)
    file0.close()

#####peptides for substitution
output = open(path_MCL+'NimbleGen_substitution_pep.csv','w+')
single_set = IEDB_set
non_c_c =0
output.write('Peptide sequence,Frequency in the database')
for key0, value0 in core_counter.iteritems():
    if value0 > 1:
        single_set.add(key0)
for key0, value0 in non_core_counter.iteritems():
    if value0 > 0:
        single_set.add(key0)
        non_c_c += 1
print('Total substitution candidates= '+str(len(single_set)))
for x in single_set:
    output.write(x+',')
    if x in core_counter:
        output.write(str(core_counter[x])+'\n')
    elif x in non_core_counter:
        output.write(str(non_core_counter[x])+'\n')
    else:
        print(x)
        output.write('Negative_control\n')        
output.close()

#########print peptides for single experiments
set2 = pep_set | core_set
print('Total unique peptides='+str(len(set2)))
output = open(path_MCL+'NimbleGen_single_test.csv','w+')
for x in set2:
    output.write(x+'\n')
output.close()

#######
#print('Total core peptides='+str(c_core))
#print('Total MCL pep overlapping'+str(len(MCL_set)))
#print('Total IEDB pep overlapping'+str(len(IEDB_set)))
#print core_counter
#print non_core_set
print('Total non core peptides='+str(non_c_c))
    
    