#this script:
#1)make a gene dictonary from the raw MS data (rawer) (once)
#1.5)Load gene dictonary (peptide to gene)
#2)get gene list from predicted list and processed MS data
#3)create three sets of genes, in MS only, in Mutant only, intersection
from utilities import *
import cPickle as pickle
import fileinput
import re
dict_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/MS_raw/HLA_MS_gene_list.dict'
path_MS = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/20151002_MCL_MHC_profiling_export_db/20151002_MCL_db_v4_export/'
path_Mut = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/Neoantigen_binding_predictions/'

def makegene_dict():
    dict0 = dict()
    for filename0 in fileinput.input():
        file0 = open(filename0.rstrip(),'r')
        head_s = True
        for line0 in file0:
            if head_s:
                head_s=False
            else:
                line0 = line0.rstrip().split(',')
                line0[13] = line0[13].split('.')[1]
                line0[13] = re.sub(r'[^\w]', '', line0[13])
                if '"' in line0[15]:
                    line0[15] = line0[15].split('"')[1]
                #13 pepetide
                #15 gene symbol
                if line0[13] in dict0:
                    if not dict0[line0[13]] == line0[15]:
                        print [filename0,line0[13],line0[15],'new:',dict0[line0[13]]]
                else:
                    dict0[line0[13]]=line0[15]
        file0.close()
    pickle.dump(dict0,open('HLA_MS_gene_list.dict','w+'))
                    
def main(pid):
    #load a dictionary
    #gene_dict = pickle.load(open(dict_file,'r'))
    #get gene list from files
    #get mutated genes from mut data
    filename0 = path_Mut+'MCL'+pid+'.table_of_binding_peptides.out'
    #column 2, separated by space
    mut_gene_list = read_col(filename0,' ',2,True)
    mut_gene_set = set(mut_gene_list)
    #get observed genes from MS data
    #column 0, separated by ,
    for type0 in ['MHC1','MHC2']:
        filename0 = path_MS+'MCL'+pid+'_'+type0+'.csv'
        ms_gene_list = read_col(filename0,',',2,True)
        #print ms_gene_list
        ms_gene_set = set()
        for x in ms_gene_list:
            if len(x.split('|'))>2:
                #print x
                ms_gene_set.add(x.split('|')[2])
                #print ms_gene_list[i]+'is not in the dictionary'
        #set operation
        print ms_gene_set
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
        print 'MCL'+pid+'\t'+type0+'\t'+str(len(mut_gene_set & ms_gene_set))+'\t'+str(len(mut_gene_set - (mut_gene_set & ms_gene_set))) + '\t' + str(len(ms_gene_set - (mut_gene_set & ms_gene_set)))

    
#makegene_dict()
#main()
print 'patientID\tMHCIorII\tOverlapedGene\tGenes_predicted_only\tGenes_MS_only'
for line0 in fileinput.input():
    line0 = line0.rstrip()
    main(line0) 
