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
                    
def main(pid):
    #load a dictionary
    gene_dict = pickle.load(open(dict_file,'r'))
    #get gene list from files
    #get mutated genes from mut data
    filename0 = path_Mut+'MCL'+pid+'.table_of_binding_peptides.out'
    #column 2, separated by space
    mut_gene_list = read_col(filename0,' ',2,True)
    mut_frag_list = read_col(filename0,' ',7,True)
    gene_dict_pred = dict()
    for i in range(0,len(mut_gene_list)):
        gene_dict_pred[mut_frag_list[i]] = mut_gene_list[i]
    mut_frag_set = set(mut_frag_list)
    mut_frag_list1 = []
    mut_frag_list2 = []
    for x in mut_frag_set:
        #from get_frag
        #for each line in mut_frag_list, [0] is the whole string, [1:] is the fragment with length n1 or n2
        #print x
        mut_frag_list1.append(get_frag(x,n1))
        mut_frag_list2.append(get_frag(x,n2))
    #get observed genes from MS data
    #column 0, separated by ,
    for type0 in ['MHC1','MHC2']:
        filename0 = path_MS+'MCL'+pid+'_'+type0+'.csv'
        print 'MCL'+pid+'\t'+type0
        ms_frag_list = read_col(filename0,',',0,True)
        ms_gene_list = read_col(filename0,',',2,True)      
        ms_gene_list = get_gene_from_ms(ms_gene_list,gene_dict)
        #print len(ms_frag_list)
        #print len(ms_gene_list)
        #print ms_gene_list
        #print ms_gene_set
        if type0 == 'MHC1':
            mut_frag_list = mut_frag_list1
        else:
            mut_frag_list = mut_frag_list2
        #print mut_frag_list
        count = 0
        for x in ms_frag_list:
            for y in mut_frag_list:
                for i in range(1,len(y)):
                    similar_score = Levenshtein.ratio(x,y[i])
                    if y[0] in gene_dict_pred:
                        pre_gene0 = gene_dict_pred[y[0]]
                    else:
                        print y[0]
                        pre_gene0 = 'No_gene'  
                    if ms_gene_list[count] == pre_gene0:
                        similar_score_str = str(similar_score)
                    else:
                        #print [ms_gene_list[i],pre_gene0]
                        similar_score_str = str(similar_score)+'*'
                    if similar_score>cut_off_high:
                        print x+'\t'+y[i]+'\t'+y[0]+'\t'+ms_gene_list[count]+'\t'+pre_gene0+'\t'+similar_score_str
                    elif similar_score>cut_off:
                        if not similar_score_str[-1] == '*':
                            print x+'\t'+y[i]+'\t'+y[0]+'\t'+ms_gene_list[count]+'\t'+pre_gene0+'\t'+similar_score_str
                        #else:
                        #    print x+'\t'+y[i]+'\t'+y[0]+'\t'+ms_gene_list[count]+'\t'+pre_gene0+'\t'+similar_score_str
            count = count + 1
        #set operation
        '''
        print 'MCL'+pid+type0
        print 'Not in the dictionary'
        print counter
        print 'Mutated genes detected by MS'
        print len(mut_gene_set & ms_gene_set)
        print 'Mutated genes not detected by MS'
        print len(mut_gene_set - (mut_gene_set & ms_gene_set))
        print 'Non-mutated genes detected by MS'
        print len(ms_gene_set - (mut_gene_set & ms_gene_set))
        '''
        #print 'MCL'+pid+'\t'+type0+'\t'+str(counter)+'\t'+str(len(mut_gene_set & ms_gene_set))+'\t'+str(len(mut_gene_set - (mut_gene_set & ms_gene_set))) + '\t' + str(len(ms_gene_set - (mut_gene_set & ms_gene_set)))
def plot_sets():
    a = 1
    
#makegene_dict()
#main()
print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
for line0 in fileinput.input():
    line0 = line0.rstrip()
    main(line0) 
