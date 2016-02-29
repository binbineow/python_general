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
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
RNA_cut_off = 1000
filter_s = True
def plot_two_set(set1,set2,set1_name,set2_name,title):
    plt.figure()
    venn2([set1,set2],set_labels = (set1_name, set2_name))
    plt.title(title)
    save(title)

def plot_two_set_num(Ab,aB,AB,set1_name,set2_name,title):
    plt.figure()
    venn2(subsets=(Ab,aB,AB),set_labels = (set1_name, set2_name))
    plt.title(title)
    save(title)

def remove_merged(list0):
    for i in range(0,len(list0)):
        if ',' in list0[i]:
            list0[i] = list0[i].split(',')[0]
        if ';' in list0[i]:
            list0[i] = list0[i].split(';')[0]
        if 'merged_' in list0[i]:
            list0[i] = list0[i].split('_')[1]
    return list0

def remove_low_RNA(set0,print_s=False):
    set_remove = set()
    global set_unmapped
    for x in set0:
        if x in dictRNA:
            if float(dictRNA[x]) < RNA_cut_off:
                set_remove.add(x)
                print x
        else:
            if print_s:
                print(x)
            set_unmapped.add(x)
            set_remove.add(x)
    set0 -= set_remove
    return set0
    
          
def main(pid):
    global set1
    global set2
    global set_m
    global set1_share
    global set2_share
    
    #get gene list from files
    #get mutated genes from mut data
    filename0 = path_Mut+'MCL'+pid+'.table_of_binding_peptides.out'
    #column 2, separated by space
    mut_gene_list = read_col(filename0,' ',2,True)
    mut_gene_set = set(remove_merged(mut_gene_list))
    print(pid+str(len(mut_gene_set)))
    if filter_s:
        mut_gene_set = remove_low_RNA(mut_gene_set)
    set_m = set_m | mut_gene_set
    #print set_m
    for type0 in ['MHC1','MHC2']:
        print(type0)
        #outfile0 = open(path_geneFrag+'MCL'+pid+'_'+type0+'mutFrag.csv','w+')
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
        ms_gene_set = set(ms_gene_list)
#        if filter_s:
#            ms_gene_set = remove_low_RNA(ms_gene_set,print_s=True)
        if type0 == 'MHC1':
            set1 = set1 | ms_gene_set
            set1_share = set1_share | (mut_gene_set & ms_gene_set)
        else:
            set2 = set2 | ms_gene_set
            set2_share = set2_share | (mut_gene_set & ms_gene_set)
        #calculate patient-specific gene overlap between mutation and total gene detected
        #plot_two_set(mut_gene_set,ms_gene_set,'Genes mutated \nin MCL'+pid,'Genes in \nMCL'+pid+' '+type0+' Ligandome','plot/'+type0+'_'+pid+'_Vennv2')
        #calculate fragement composition from mutated genes between mutation and total gene detected
        
        
        
#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
#count1 = Counter()
#count2 = Counter()
set_unmapped = set()
set_m = set('IGHM')
set1 = set()
set1_share = set()
set2 = set()
set2_share = set()
#load a dictionary
dictRNA=pickle.load(open(dictRNA_file,'r'))
dictRNA['IGHM'] = 1000
dictRNA['KMT2D'] = 5
dictRNA['KMT2B'] = 10
dictRNA['KMT2A'] = 5
gene_dict = pickle.load(open(dict_file,'r'))
for line0 in fileinput.input():
    line0 = line0.rstrip()
    main(line0) 
print(len(set_m))
print(len(set1))
print(len(set1_share))
print set_unmapped
print(str(len(set_unmapped)))
end0 = 'v2.2'
if filter_s:
    end0 = end0+'_RPKM'+'_'+str(RNA_cut_off) 
plot_two_set_num(len(set_m)-len(set1_share),len(set1)-len(set1_share),len(set1_share),'Genes with mutations \nin all patients','Genes detected \nin MHCI Ligandome','plot/MHCI_All_patients_Venn_'+end0)
plot_two_set_num(len(set_m)-len(set2_share),len(set2)-len(set2_share),len(set2_share),'Genes with mutations \nin all patients','Genes detected \nin MHCII Ligandome','plot/MHCII_All_patients_Venn_'+end0)
# print 'All IGHM MHC1'
# print set1
# print 'All IGHM MHC2'
# print set2
