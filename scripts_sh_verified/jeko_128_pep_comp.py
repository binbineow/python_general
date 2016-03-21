 #importing
import fileinput
from utilities import *

#dict_g_u = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_gene_to_uni.dict'))
#dict_u_g = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_uni_to_gene.dict'))
dict_jeko = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/jeko_pep_countv2.dict'))
#dict_L128 = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/L128_pep_countv2.dict'))



def get_name(str):
    if str in dict_u_g:
        name0 = dict_u_g[str]
    else:
        name0 = str      
    if '_' in name0:
        name0 = name0.split('_')[0]        
    return name0

#load data

MCL_data = pickle.load(open('MCL_data11_18_2015v1.2_UP.dict','r'))
MCL_data_jeko = list(MCL_data['MCLJeko']['MHC1_gene'])


#main
filename = 'Jeko_MHC1_gene_pep'
hist_num = []
for gene0 in MCL_data_jeko:
    if gene0 in dict_jeko:
        num0 = dict_jeko[gene0]
        hist_num.append(num0)
    else:
        hist_num.append(0)
           
plot_counter(hist_num,filename,filename+' Histgram','Gene frequency','peptide abundance in Jeko',True,200)

MCL_data_jeko = list(MCL_data['MCLJeko']['MHC2_gene'])


#main
filename = 'Jeko_MHC2_gene_pep'
hist_num = []
for gene0 in MCL_data_jeko:
    if gene0 in dict_jeko:
        num0 = dict_jeko[gene0]
        hist_num.append(num0)
    else:
        hist_num.append(0)
           
plot_counter(hist_num,filename,filename+' Histgram','Gene frequency','peptide abundance in Jeko',True,200)
            

