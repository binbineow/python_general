#this script:
#1)make a gene dictonary from the raw MS data (rawer) (once)
#1.5)Load gene dictonary (peptide to gene)
#2)get gene list from predicted list and processed MS data
#3)create three sets of genes, in MS only, in Mutant only, intersection
from utilities import *
import cPickle as pickle
import fileinput
from collections import Counter
dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/UniprotID.dict'
path_MS = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/20151002_MCL_MHC_profiling_export_db/20151002_MCL_db_v4_export/'
path_Mut = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/Neoantigen_binding_predictions/'
path_geneFrag = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/FragFromMut/'


def get_six_comp(str0):
    if str0[0] =='"':
        str0 = str0[1:-1]
    list0 = []
    list0.append(str0[0])
    list0.append(str0[2])
    list0.append(str0[0:3])
    list0.append(str0[-3])
    list0.append(str0[-1])
    list0.append(str0[-3:])
    return list0
                    
def main():
    list_data =[]
    print 'A.BXXXXXXC.D 1=A 2=B 3=A.B 4=C 5=D 6=C.D'
    for line0 in fileinput.input():
        #print line0[0:5]
        if not line0[0:5] == 'ScanF':
            line0 = line0.split(',')
            print line0[-3]
            if len(line0[-3])>4:
                list0 = get_six_comp(line0[-3])
                list_data.append(list0)
    #out_count_list = []
    for i in range(0,6):
        print i+1
        count0 = Counter()
        for list0 in list_data:
            count0[list0[i]] +=1
        print count0.most_common(len(count0))
        
        
        

        #print 'MCL'+pid+'\t'+type0+'\t'+str(counter)+'\t'+str(len(mut_gene_set & ms_gene_set))+'\t'+str(len(mut_gene_set - (mut_gene_set & ms_gene_set))) + '\t' + str(len(ms_gene_set - (mut_gene_set & ms_gene_set)))

#makegene_dict()
#main()
#print 'MS\tPredicted\tFrom_mutation\tMS_gene\tPredicted_gene\tSimilarity'
main()
