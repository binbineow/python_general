import cPickle as pickle
import fileinput
from utilities import *
import pandas as pd

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

def MS_clean(frag_list,gene_list):
    #get rid of any SPA 
    for i in range(len(gene_list)-1,0,-1):
        if '|SPA|' in gene_list[i] or len(gene_list[i])<6:
            del gene_list[i]
            del frag_list[i]
    return [frag_list,gene_list]

def is_human(gene0):
    if '|SPA|' in gene0 or len(gene0) < 7:
        return False
    else:
        return True

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


#get list from mhc peptide ms data file
def get_list_from_ms(file0):
    list_pep = []
    list_gene = []
    list_gn = []
    i = 0
    for line0 in open(file0,'rU'):
        if i > 0:
            line_ori = line0
            line0 = line0.rstrip().split(',')
            if is_human(line0[2]):
                if 'GN=' in line0:
                    list_gn.append(line_ori.split('GN=')[1].split(' ')[0])
                else:
                    list_gn.append('')
                list_gene.append(line0[2])
                list_pep.append(line0[0])
        else:
            i +=1
    return list_pep,list_gene,list_gn


#convert gene name format to uniprot id and gn
#possible 
#P0CW22|RS17L|1
#tr|D7RIG0|D7RIG0|1
#sp|P04114|APOB|1
#MCL041IgH|idiotype
#GN=MCL041_IgH (related)
def get_gene_from_msv2(list0,list_gn):
    list_sp = []
    list_gn_v2 = []
    for name0,name_gn in zip(list0,list_gn):
        if 'idotype' in name0:
            gene0 = name0.split('|')[0]
            gn_backup = gene0[0:len('MCL041')]+'_'+gene0[-3:]
        name0 = name0.split('|')
        if len(name0) == 4:
            gene0 = name0[1]
            gn_backup =  name0[2]
        elif len(name0) == 3:
            gene0 = name0[0]
            gn_backup = name0[1]
        elif len(name0) == 2:
            gene0 = name0[0]
        else:
            gene0 = name0[0]
            print(gene0)
        list_sp.append(gene0)
        if name_gn == '':
            list_gn_v2.append(gn_backup)
        else:
            list_gn_v2.append(name_gn)
    #print name0
    #print list_sp
    return list_sp, list_gn_v2

                      
#load dictionary
gene_dict = pickle.load(open(dict_file,'r'))
#dictRNA=pickle.load(open(dictRNA_file,'r'))
#initiate 
mcl_mhc1_all = set()
mcl_mhc2_all = set()
mcl_mut_original = set()
mcl_wt_pep = set()
mcl_mut_pep = set()

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
        data_mut = pd.read_csv(filename0,' ')
        mut_gene_list = list(data_mut['Gene'])
        mut_original_list  = list(data_mut['Full_length_peptide']) 
        mut_pep_mhc1 = list(data_mut['Mut_peptide']) 
        wt_pep_mhc1 = list(data_mut['Wt_peptide']) 
        mut_gene_list = remove_merged(mut_gene_list)
        data_MCL[pid]['mutation']=mut_gene_list
        data_MCL[pid]['mutation'].append('IGHM')
        data_MCL[pid]['mutation_pep_original'] = mut_original_list
        data_MCL[pid]['mutation_pep_mhc1'] = mut_pep_mhc1
        data_MCL[pid]['wt_pep_mhc1'] = wt_pep_mhc1
        #update total list
        mcl_mut_original = mcl_mut_original | set(mut_original_list)
        mcl_wt_pep = mcl_wt_pep | set(wt_pep_mhc1)
        mcl_mut_pep = mcl_mut_pep | set(mut_pep_mhc1)
        #get MHC1 or MHC2 fragments and genes
        for type0 in ['MHC1','MHC2']:
            filename0 = path_MS+pid+'_'+type0+'_20151109_SQLPowerTool.csv'
            ms_frag_list, ms_gene_list, ms_gn_list = get_list_from_ms(filename0)
            #print(len(ms_gene_list))
            #ms_frag_list = read_col(filename0,',',0,True)
            #ms_gene_list = read_col(filename0,',',2,True)
            #[ms_frag_list,ms_gene_list] = MS_clean(ms_frag_list,ms_gene_list)
            ms_gene_list,ms_gn_list = get_gene_from_msv2(ms_gene_list,ms_gn_list)
            print(len(ms_gene_list))
            data_MCL[pid][type0+'_frag'] = ms_frag_list
            data_MCL[pid][type0+'_gene'] = ms_gene_list
            data_MCL[pid][type0+'_gn'] = ms_gn_list
            if type0 == 'MHC1':
                mcl_mhc1_all = mcl_mhc1_all | set(ms_frag_list)
            else:
                mcl_mhc2_all = mcl_mhc2_all | set(ms_frag_list)

data_MCL['pooled']['all_mhc1'] = list(mcl_mhc1_all)
data_MCL['pooled']['all_mhc2'] = list(mcl_mhc2_all)
data_MCL['pooled']['all_original_mut'] = list(mcl_mut_original)
data_MCL['pooled']['all_wt_mhc1'] = list(mcl_wt_pep)[1:]
data_MCL['pooled']['all_mut_mhc1'] = list(mcl_mut_pep)
data_MCL['pooled']['key_words'] = ['all_mhc1','all_mhc2','all_original_mut','all_wt_mhc1','all_mut_mhc1']
pickle.dump(data_MCL,open('MCL_2015_data_4_12_2017.dict','wb+'))


    