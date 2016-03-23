 #importing

from utilities import *
from scipy.stats import pearsonr
path0 = '../images/'
import math 

def square0(list):
    return [i ** 2 for i in list]

def plot_scatter(x,y,x_name,y_name,title,path0=''):
    plt.figure()
    plt.plot(x,y,'b.')
    plt.ylabel(y_name)
    plt.xlabel(x_name)
    plt.ylim([0,max(y)+1])
    plt.title(title)
    save(path0+title)

dict_g_u = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_gene_to_uni.dict'))
#dict_u_g = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_uni_to_gene.dict'))
dict_jeko = pickle.load(open('/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/jeko_pep_countv3.dict'))
#dict_L128 = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/L128_pep_countv2.dict'))

#_0 means adding a zero to pep list if no entry can be found in the jeko_pep_countv3 (assuming zero pepetide)
jeko_rna_0 = []
jeko_rna_no0 = []
jeko_pep_0 = []
jeko_pep_no0 = []

for line0 in fileinput.input():
    line0 = line0.rstrip().split(',')
    if not 'Gene name' in line0[0]:
        if line0[0] in dict_g_u:
            rna_gene0 = dict_g_u[line0[0]]
            rna_value0 = math.log(float(line0[1])+0.0001,10)
            if rna_gene0 in dict_jeko:
                jeko_rna_no0.append(rna_value0)
                jeko_pep_no0.append(dict_jeko[rna_gene0])
                jeko_rna_0.append(rna_value0)
                jeko_pep_0.append(dict_jeko[rna_gene0])
            else:
                jeko_rna_0.append(rna_value0)
                jeko_pep_0.append(0.0)
            
#no0
x_name = 'RNASeq log10(RPKM) of Jeko Gene'
y_name = 'Peptide emPIA of Jeko Gene'
title0 = 'Regression_between_RNA_expression_and_peptide_abundance_\nin_Jeko_cell_line_(gene_emPIA_>0)'
r0 = str(pearsonr(jeko_rna_no0,square0(jeko_pep_no0))[0])
print('Regression R without adding zero peptide='+r0)
plot_scatter(jeko_rna_no0,jeko_pep_no0,x_name+' (R='+r0+')',y_name,title0,path0)   
#0
x_name = 'RNASeq log10(RPKM) of Jeko Gene'
y_name = 'Peptide emPIA of Jeko Gene'
title0 = 'Regression_between_RNA_expression_and_peptide_abundance_\nin_Jeko_cell_line'
r0 = str(pearsonr(jeko_rna_0,square0(jeko_pep_0))[0])
print('Regression R without adding zero peptide='+r0)
plot_scatter(jeko_rna_0,jeko_pep_0,x_name+' (R='+r0+')',y_name,title0,path0)   
