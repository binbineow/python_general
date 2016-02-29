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
from collections import defaultdict
dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/UniprotID.dict'
path_MS = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/20151002_MCL_MHC_profiling_export_db/20151002_MCL_db_v4_export/'
path_Mut = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/Neoantigen_binding_predictions/'
n1 = 9
n2 = 12
cut_off_high = 0.9
cut_off = 0.7

def get_frag(str0,n):
    list_out = []
    list_out.append(str0)
    if len(str0) > n:
        for i in range(0,len(str0)-n):
            list_out.append(str0[i:i+n])
    else:
        list_out.append(str0)
    return list_out
def get_total_and_common(gene_set_list):
    #get total gene number
    set_total = set()
    for x in gene_set_list:
        set_total = set_total | x
        #get common gene number
    set_common = gene_set_list[0] 
    for x in gene_set_list:
        set_common = set_common & x 
    print set_common
    return [set_total,set_common]

def plot_hist(gene_set_list):
    counter = defaultdict(int)
    for s in gene_set_list:
        for el in s:
            counter[el] += 1
    #print(counter)
    counter2 = defaultdict(int)
    for key,value in counter.iteritems():
        counter2[value] += 1
    for key,value in counter2.iteritems():
        print str(key)
    for key,value in counter2.iteritems():
        print str(value)
             
def main(pid_list):
    #load a dictionary
    gene_dict = pickle.load(open(dict_file,'r'))
    #get gene list from files
    #get mutated genes from mut data
    #get observed genes from MS data
    #column 0, separated by ,
    #total_gene_set = set()
    for type0 in ['MHC1','MHC2']:
        #get genelist
        print type0
        gene_set_list = []
        frag_set_list = []
        for pid in pid_list:
            filename0 = path_MS+'MCL'+pid+'_'+type0+'.csv'
            ms_frag_list = read_col(filename0,',',0,True)
            ms_gene_list = read_col(filename0,',',2,True)      
            ms_gene_list = get_gene_from_ms(ms_gene_list,gene_dict)
            gene_set_list.append(set(ms_gene_list))
            frag_set_list.append(set(ms_frag_list))
        [set_total,set_common] = get_total_and_common(gene_set_list)
        print 'Total gene present in'+type0+': '+str(len(set_total))
        print 'The common gene in '+type0+': '+str(len(set_common))
        print set_common
        [set_total,set_common] = get_total_and_common(frag_set_list)
        print 'Total fragments present in'+type0+': '+str(len(set_total))
        print 'The common fragments in '+type0+': '+str(len(set_common))    
        print set_common    
        # individual set
        #for i in range(0,len(pid_list)):
        #    print 'Not_in_common for '+pid_list[i]+': '+str(len(gene_set_list[i]-set_common))
        plot_hist(gene_set_list)
        plot_hist(frag_set_list)
        
        

#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
pid_list = []
for line0 in fileinput.input():
    line0 = line0.rstrip()
    pid_list.append(line0)
main(pid_list)
