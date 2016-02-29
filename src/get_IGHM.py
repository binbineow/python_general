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
                    
def main():
    pid_list = []
    print 'sample\tMHCI_peptides\tMHCI_IGHM_peptides\tprecent\tMHCII_peptides\tMHCII_IGHM_peptides\tpercent'
    for line0 in fileinput.input():
        line0 = line0.rstrip()
        pid_list.append(line0)
    #load a dictionary    
    gene_dict = pickle.load(open(dict_file,'r'))
    count1 = Counter()
    count2 = Counter()
    for pid in pid_list:
        #get gene list from files
        #get mutated genes from mut data
        #filename0 = path_Mut+'MCL'+pid+'.table_of_binding_peptides.out'
        #column 2, separated by space
        #mut_gene_list = read_col(filename0,' ',2,True)
        out_line = pid+'\t'
        for type0 in ['MHC1','MHC2']:
            #outfile0 = open(path_geneFrag+'MCL'+pid+'_'+type0+'mutFrag.csv','w+')
            filename0 = path_MS+'MCL'+pid+'_'+type0+'.csv'
            #outfile0.write('#MCL'+pid+'\t'+type0+'\n')
            #print '#MCL'+pid+'\t'+type0
            ms_frag_list = read_col(filename0,',',0,True)
            ms_gene_list = read_col(filename0,',',2,True)
            ms_gene_list = get_gene_from_ms(ms_gene_list,gene_dict)
            nt = len(ms_frag_list)
            out_line=out_line+str(nt)+'\t'
            #out_dict0 = dict()
            n0 = 0
            for i in range(0,len(ms_gene_list)):
                if ms_gene_list[i] == "IGHM":
                    #print ms_frag_list[i] 
                    if type0 == 'MHC1':
                        count1[ms_frag_list[i]]+=1
                        n0 = n0+1
                    else:
                        count2[ms_frag_list[i]]+=1
                        n0 = n0+1
            out_line=out_line+str(n0)+'\t'+str(float(n0)/float(nt)*100)+'\t'
        print out_line
    count1 = count1.most_common(len(count1))
    count2 = count2.most_common(len(count2))
    print count1
    print count2
    print 'MHC1'
    print len(count1)
    for key in count1:
        print key[0]
    print 'MHC2'
    print len(count2)
    for key in count2:
        print key[0]
          
                
        #print 'MCL'+pid+'\t'+type0+'\t'+str(counter)+'\t'+str(len(mut_gene_set & ms_gene_set))+'\t'+str(len(mut_gene_set - (mut_gene_set & ms_gene_set))) + '\t' + str(len(ms_gene_set - (mut_gene_set & ms_gene_set)))

#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
main()
