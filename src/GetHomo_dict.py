#template for general purposes
import fileinput
from utilities import *
#import subprocess

#variable
pathin=''
pathout = ''
file_in_name = ''
file_out_name = ''
des_list = []
dict_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/uniprot_geneid_homo.dict'
dict_len_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/uniprot_gene_index.dict'
coutner_genebinding_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/gene_binding.counter'

def make_uniprot_list():
    for line in fileinput.input():
        if '>' in line:
            if 'OS=Homo sapiens' in line and len(line.rstrip().split(' ')[0].split('|'))>2:
                print line.rstrip().split(' ')[0].split('|')[1]

def make_uniprot_homo():
    dict0 = make_any_dict('uniprot_geneid_homo.dict',False,'\t',True)
    print dict0



def main(outfilename):
    dict_gene = pickle.load(open(dict_gene_path,'r'))
    reading0 = False
    gene_line = ''
    output0 = ''
    pass0 = True
    count0 = 0
    dict_len_gene = dict()
    #output1 = []
    for line in fileinput.input():
        line = line.rstrip()
        if '>' in line:
            if 'OS=Homo sapiens' in line:
                reading0 = True
                if len(output0)>1 and pass0:
                    if len(gene_line.rstrip().split(' ')[0].split('|'))>2:
                        gene_name = get_gene_from_fasta(line,dict_gene)
                        for i in range(count0,count0+len(output0)+1):
                            dict_len_gene[i] = gene_name 
                    count0 = count0 + len(output0)
                gene_line = line    
                output0 = ''
                pass0 = True
            else:
                reading0 = False
        elif reading0:
            if pass0:
                for item in 'BJOUXZ':
                    if item in line:
                        pass0 = False
                        break
            if pass0:
                output0 = output0+line
    if len(output0)>1 and pass0:
        if len(gene_line.rstrip().split(' ')[0].split('|'))>2:
            gene_name = get_gene_from_fasta(line,dict_gene)
            for i in range(count0,count0+len(output0)+1):
                dict_len_gene[i] = gene_name
    #print dict_len_gene
    pickle.dump(dict_len_gene,open(outfilename,'w+'))

def pep_to_gene(header_s):
    dict_len_to_gene = pickle.load(open(dict_len_path,'r'))
    gene_list = []
    for line in fileinput.input():
        if header_s:
            header_s = False
        else:
            #print line.rstrip().split('\t')[6]
            if float(line.rstrip().split('\t')[6]) < 100:
                begin0 = float(line.rstrip().split('\t')[2])
                end0 = float(line.rstrip().split('\t')[3])
                if begin0 in dict_len_to_gene and end0 in dict_len_to_gene:
                    #print [begin0,dict_len_to_gene[begin0]]
                    if dict_len_to_gene[begin0] == dict_len_to_gene[end0]:
                        gene_list.append(dict_len_to_gene[begin0])  
            else:
                break
    gene_list = [ element0 for element0 in gene_list if not element0 == 'No_gene']
    #print gene_list
    #plot_counter(gene_list,'Hist_of_genefrag_predicted_v2','Histagrm','Frequency','peptide number per gene',True,50)
    counter1 = get_counter(gene_list)
    pickle.dump(counter1,open(coutner_genebinding_path,'w+'))
    
    
def plto_relationship():
    gene_binding_c   #binding peptide per gene
    RNA_dict #RNA expression per gene
    counter_pepgene #peptide counter among 7 patients of interest
for key,value in counter_pepgene.iteritems():
    #key_list.append(key)
    if key in RNA_dict:
        if value > 1 and RNA_dict[key] < 0.01:
            print [key,value]
    else:
        print key
for key,value in counter_pepgene.iteritems():
    #key_list.append(key)
    if key in RNA_dict:
        if value > 1 and gene_binding_c[key] == 0:
            print [key,value]
    else:
        print key
        
pep_to_gene(True)

