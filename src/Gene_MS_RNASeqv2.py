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
import math
from collections import defaultdict
dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/UniprotID.dict'
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
path_MS = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/20151002_MCL_MHC_profiling_export_db/20151002_MCL_db_v4_export/'
path_Mut = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/Neoantigen_binding_predictions/'
n1 = 9
n2 = 12
cut_off_high = 0.1
cut_off = 0.01

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
    #print set_common
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
    #print gene_dict
    RNA_dict = pickle.load(open(dictRNA_file,'r'))
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
            ms_stuff = zip(ms_frag_list,ms_gene_list)
            #filter out undesired genes (SPA and 'KRT')
            filtered = [x for x in ms_stuff if (not 'KRT' in x[1] ) and (not 'SPA' == x[1].upper())]
            ms_gene_list = [x[1] for x in filtered]
            ms_frag_list = [x[0] for x in filtered]
            gene_set_list.append(set(ms_gene_list))
            frag_set_list.append(set(ms_frag_list))
        [set_total,set_common] = get_total_and_common(gene_set_list)
        print 'Total genes present in '+type0+': '+str(len(set_total))
        n_high = 0
        n_low = 0
        for x in set_total:
            if x in RNA_dict:
                if float(RNA_dict[x]) < cut_off_high:
                    n_high += 1
                    print x
                if float(RNA_dict[x]) < cut_off:
                    n_low +=1
                    
        print 'Total genes<'+str(cut_off_high)+'RPKM: '+str(n_high)
        print  'Percentage='+str(float(n_high)/float(len(set_total)))
        print 'Total genes<'+str(cut_off)+'RPKM: '+str(n_low)
        print  'Percentage='+str(float(n_low)/float(len(set_total)))
        #print 'The common gene in '+type0+': '+str(len(set_common))
        #print set_common
        #print out gene and its associated RNA expression
        '''
        for x in set_common:
            if x in RNA_dict:
                print x+'\t'+str(RNA_dict[x])
        rna_out_total = [math.log(float(RNA_dict[x])+0.0001) for x in set_total if x in RNA_dict ]
        plt.hist(rna_out_total)
        plt.title("RNASeq Expression Histogram")
        plt.xlabel("RPKM")
        plt.ylabel("Frequency")
        save(type0+'RNASeq_outv1','png')
        [set_total,set_common] = get_total_and_common(frag_set_list)
        print 'Total fragments present in'+type0+': '+str(len(set_total))
        print 'The common fragments in '+type0+': '+str(len(set_common))    
        print set_common    
        # individual set
        #for i in range(0,len(pid_list)):
        #    print 'Not_in_common for '+pid_list[i]+': '+str(len(gene_set_list[i]-set_common))
        plot_hist(gene_set_list)
        plot_hist(frag_set_list)
        '''
        
        
        

#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
pid_list = []
for line0 in fileinput.input():
    line0 = line0.rstrip()
    pid_list.append(line0)
main(pid_list)
