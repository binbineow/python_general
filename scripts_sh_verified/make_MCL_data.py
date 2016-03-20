import cPickle as pickle
import fileinput
from utilities import *

dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/UniprotID.dict'
path_MS = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/MS_11_11/'
path_Mut = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/Neoantigen_binding_predictions/'
path_geneFrag = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/FragFromMut/'
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'

from collections import defaultdict
def dumb(): return defaultdict(list)
data_MCL = defaultdict(dumb)
#first patient id
#second data type: MHC1 peptides, MHC1 peptide genes, MHC2 p, MHC2 pn, mutation info, HLA typing
#all lists
dictRNA = pickle.load(open(dictRNA_file,'r'))
dictRNA['idiotype'] = 1000
pickle.dump(dictRNA,open(dictRNA_file,'w'))

def MS_clean(frag_list,gene_list):
    #get rid of any SPA 
    for i in range(len(gene_list)-1,0,-1):
        if '|SPA|' in gene_list[i] or len(gene_list[i])<6:
            del gene_list[i]
            del frag_list[i]
    return [frag_list,gene_list]

#get info from HLA lines
def get_HLA_typing (str0):
    str0 = str0.rstrip()
    if 'Locus' in str0:
        return [False,0,0,0]
    else:
        str0 = str0.split('\t')
        pid = str0[0]
        for i in range(2,4):
            part1 = 'HLA-'+str0[i].split('*')[0]
            part2 = '*'+str0[i].split('*')[1][0:5]
            if i == 2:
                alle1 = part1+part2
            else:
                alle2 = part1+part2
        return [True,pid,alle1,alle2]

def get_gene_from_ms_UP(list0,pid=''):
    list_return = []
    #if idiotype, return IGHM (P01871)
    #extract uniprot ID if appropriate format
    #or return No_gene
    for x in list0:
        if 'idiotype' in x and pid in x:
            list_return.append('P01871')
        elif len(x.split('|'))==3:
            list_return.append(x.split('|')[0])
                #print x
        elif len(x.split('|'))==4:
            list_return.append(x.split('|')[1])
        else:
            list_return.append('No_gene')
    return list_return

                      
#load dictionary
#gene_dict = pickle.load(open(dict_file,'r'))
#dictRNA=pickle.load(open(dictRNA_file,'r'))

for line0 in fileinput.input():
    if len(line0)> 15:
        #get HLA typing: 0,1 HLA-A 2,3 HLA-B 4,5 HLA-C 6,7 HLA-DQA1 8,9 HLA-DQB1 10,11 HLA-DRB1 12,13  
        [info_line,pid,alle1,alle2] = get_HLA_typing(line0)
        if info_line:
            data_MCL[pid]['HLA_typing'].append(alle1)
            data_MCL[pid]['HLA_typing'].append(alle2)
    else:
        pid = 'MCL'+line0.rstrip()
        data_MCL['pid']['pid'].append(pid)
        ###get mutation data
        filename0 = path_Mut+pid+'.table_of_binding_peptides.out'
        #column 2, separated by space
        mut_gene_list = read_col(filename0,' ',2,True)
        mut_gene_list = remove_merged(mut_gene_list)
        data_MCL[pid]['mutation']=mut_gene_list
        data_MCL[pid]['mutation'].append('idiotype')
        #get MHC1 or MHC2 fragments and genes
        for type0 in ['MHC1','MHC2']:
            filename0 = path_MS+pid+'_'+type0+'_20151109_SQLPowerTool.csv'
            #outfile0.write('#MCL'+pid+'\t'+type0+'\n')
            #print '#MCL'+pid+'\t'+type0+'\n' 
            ms_frag_list = read_col(filename0,',',0,True)
            ms_gene_list = read_col(filename0,',',2,True)
            [ms_frag_list,ms_gene_list] = MS_clean(ms_frag_list,ms_gene_list)
            ms_gene_list = get_gene_from_ms(ms_gene_list,gene_dict,pid=pid)
            data_MCL[pid][type0+'_frag'] = ms_frag_list
            data_MCL[pid][type0+'_gene'] = ms_gene_list

pickle.dump(data_MCL,open('MCL_data11_18_2015v1.2_UP.dict','wb+'))


    
