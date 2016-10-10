from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import pandas as pd

def count_shared(listy, listx):
    n_shared = 0
    for frag0 in listy:
        for fragx in listx:
            #if frag0 in fragx or Levenhstein.ratio(frag0,fragx) < 0.9 or fragx in frag0:
            if frag0 in fragx or fragx in frag0:
                n_shared += 1
                break
    return n_shared

mhc1_matrix = []
mhc2_matrix = []
MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict','r'))
pid_list = MCL_data['pid']['pid']
#the format: Row X, and Column Y is percentage of patient X is shared with patient Y
#mhc1 vs. mhc2
#MCL_data[x]['MHC1_frag'][i]
for mhc0 in ['MHC1','MHC2']:
    #row
    matrix0 = []
    for x in pid_list: 
        # column
        n0 = len(MCL_data[x][mhc0+'_frag'])
        row0 = []
        for y in pid_list:
            #y on x, put y first in the function
            n_shared = count_shared(MCL_data[y][mhc0+'_frag'],MCL_data[x][mhc0+'_frag'])
            
            
        matrix0.append(row0)
        
    gene1_set_list.append(set(MCL_data[x]['MC1_gene']))
    gene2_set_list.append(set(MCL_data[x]['MC2_gene']))
counter1 = count_set(gene1_set_list)
#print counter1
counter2 = count_set(gene2_set_list)
RNAdict = pickle.load(open(dictRNA_file,'r'))
#print counter2
x0 = []
y0 = []
for key, value in counter1.iteritems():
    if key in RNAdict: 
        y0.append(value)
        x0.append(math.log(RNAdict[key]))
plot_scatter(x0,y0,'Estimated RNA Expression (log)','Frequency of Gene in MHC1 Ligandome','RNA expression profiles of genes in MHC1 Ligandome')

x0 = []
y0 = []
for key, value in counter2.iteritems():
    if key in RNAdict: 
        y0.append(value)
        x0.append(math.log(RNAdict[key]))
plot_scatter(x0,y0,'Estimated RNA Expression (log)','Frequency of Gene in MHC2 Ligandome','RNA expression profiles of genes in MHC2 Ligandome')
    

'''
set1 = set()
set2 = set()
for x in MCL_data['pid']['pid']:
    set1 = set1 | set(MCL_data[x]['MHC1_gene'])
    set2 = set2 | set(MCL_data[x]['MHC2_gene'])
set1_share = pickle.load(open('MHC1_shared_MCL.set','r'))
set2_share = pickle.load(open('MHC2_shared_MCL.set','r'))
set_mut1 = pickle.load(open('filtered_mutation_MCL.set','r'))
set_mut2 = pickle.load(open('filtered_mutation_MCL.set','r'))
'''

##########heatmap generating


pickle.dump()

    
    
    
    

             