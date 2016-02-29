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
                    
def main(pid):
    global set1
    global set2
    #load a dictionary
    gene_dict = pickle.load(open(dict_file,'r'))
    #get gene list from files
    #get mutated genes from mut data
    filename0 = path_Mut+'MCL'+pid+'.table_of_binding_peptides.out'
    #column 2, separated by space
    mut_gene_list = read_col(filename0,' ',2,True)
    mut_gene_set = set(mut_gene_list)
    for type0 in ['MHC1','MHC2']:
        outfile0 = open(path_geneFrag+'MCL'+pid+'_'+type0+'mutFrag.csv','w+')
        filename0 = path_MS+'MCL'+pid+'_'+type0+'.csv'
        #outfile0.write('#MCL'+pid+'\t'+type0+'\n')
        print '#MCL'+pid+'\t'+type0+'\n' 
        ms_frag_list = read_col(filename0,',',0,True)
        ms_gene_list = read_col(filename0,',',2,True)
        ms_gene_list = get_gene_from_ms(ms_gene_list,gene_dict)
        n0 = len(mut_gene_set)
        out_dict0 = dict()
        for i in range(0,len(ms_gene_list)):
            if ms_gene_list[i] == "IGHM":
                print ms_frag_list[i] 
                if type0 == 'MHC1':
                    set1.add(ms_frag_list[i])
                else:
                    set2.add(ms_frag_list[i])
        for mut_gene0 in mut_gene_set:
            out_dict0[mut_gene0] = ''
            for i in range(0,len(ms_gene_list)):
                if mut_gene0 == ms_gene_list[i]:
                    out_dict0[mut_gene0] = out_dict0[mut_gene0]+ms_frag_list[i]+'\n'
            if out_dict0[mut_gene0] == '':
                n0 = n0 -1
        line0 = 'MCL'+pid+'_'+type0+' has fragments from '+str(n0)+' out of '+str(len(mut_gene_set))+' mutated genes\n'
        outfile0.write(line0)
        #print line0[:-1]
        for key,value0 in out_dict0.iteritems():
            if len(value0) > 2:
                outfile0.write(key+'\n')
                outfile0.write(value0)
        outfile0.close()
            
            
        
                
                
        #print 'MCL'+pid+'\t'+type0+'\t'+str(counter)+'\t'+str(len(mut_gene_set & ms_gene_set))+'\t'+str(len(mut_gene_set - (mut_gene_set & ms_gene_set))) + '\t' + str(len(ms_gene_set - (mut_gene_set & ms_gene_set)))

#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
count1 = Counter()
count2 = Counter()
set1 = set()
set2 = set()
for line0 in fileinput.input():
    line0 = line0.rstrip()
    main(line0) 
print 'All IGHM MHC1'
print set1
print 'All IGHM MHC2'
print set2
