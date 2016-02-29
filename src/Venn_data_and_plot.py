#this script:
#1)make a gene dictonary from the raw MS data (rawer) (once)
#1.5)Load gene dictonary (peptide to gene)
#2)get gene list from predicted list and processed MS data
#3)create three sets of genes, in MS only, in Mutant only, intersection
from utilities import *
import cPickle as pickle
import fileinput
import re
import Levenshtein
from collections import Counter
from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/UniprotID.dict'
path_MS = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/20151002_MCL_MHC_profiling_export_db/20151002_MCL_db_v4_export/'
path_Mut = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/Neoantigen_binding_predictions/'
path_geneFrag = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/FragFromMut/'


def get_frag(str0,n):
    list_out = []
    list_out.append(str0)
    if len(str0) > n:
        for i in range(0,len(str0)-n):
            list_out.append(str0[i:i+n])
    else:
        list_out.append(str0)
    return list_out

def plot_two_set(set1,set2,set1_name,set2_name,title):
    plt.figure()
    venn2([set1,set2],set_labels = (set1_name, set2_name))
    plt.title(title)
    save(title)
          
def main(pid):
    global set1
    global set2
    global set_m
    #load a dictionary
    gene_dict = pickle.load(open(dict_file,'r'))
    #get gene list from files
    #get mutated genes from mut data
    filename0 = path_Mut+'MCL'+pid+'.table_of_binding_peptides.out'
    #column 2, separated by space
    mut_gene_list = read_col(filename0,' ',2,True)
    mut_gene_set = set(mut_gene_list)
    set_m = set_m | mut_gene_set
    for type0 in ['MHC1','MHC2']:
        outfile0 = open(path_geneFrag+'MCL'+pid+'_'+type0+'mutFrag.csv','w+')
        filename0 = path_MS+'MCL'+pid+'_'+type0+'.csv'
        #outfile0.write('#MCL'+pid+'\t'+type0+'\n')
        print '#MCL'+pid+'\t'+type0+'\n' 
        ms_frag_list = read_col(filename0,',',0,True)
        ms_gene_list = read_col(filename0,',',2,True)
        #Convert gene list
        ms_gene_list = get_gene_from_ms(ms_gene_list,gene_dict)
        #delete any genes with 'SPA' or not mappable to the gene_dict (No_genes)
        for i in range(len(ms_frag_list)-1,0,-1):
            if ms_gene_list[i].upper() == 'SPA' or ms_gene_list[i] == 'No_gene':
                del ms_gene_list[i]
                del ms_frag_list[i]
        set_f_total = set()
        set_f_m = set()    
        ms_gene_set = set(ms_gene_list)
        if type0 == 'MHC1':
            set1 = set1 | ms_gene_set
        else:
            set2 = set2 | ms_gene_set
        #calculate patient-specific gene overlap between mutation and total gene detected
        plot_two_set(mut_gene_set,ms_gene_set,'Genes mutated \nin MCL'+pid,'Genes in \nMCL'+pid+' '+type0+' Ligandome','plot/'+type0+'_'+pid+'_Venn')
        #calculate fragement composition from mutated genes between mutation and total gene detected
        
        
        
#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
#count1 = Counter()
#count2 = Counter()
set_m = set()
set1 = set()
set2 = set()
for line0 in fileinput.input():
    line0 = line0.rstrip()
    main(line0) 
plot_two_set(set_m,set1,'Genes with mutations \nin all patients','Genes detected \nin MHCI Ligandome','plot/MHCI_All_patients_Venn')
plot_two_set(set_m,set2,'Genes with mutations \nin all patients','Genes detected \nin MHCII Ligandome','plot/MHCII_All_patients_Venn')
# print 'All IGHM MHC1'
# print set1
# print 'All IGHM MHC2'
# print set2
