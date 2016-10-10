from utilities import *

list_nonig = pickle.load(open('MCL_all_nonIg.list'))
list_nonig_gene = pickle.load(open('MCL_all_nonIG_gene.list'))
list_constant = pickle.load(open('MCL_all_IG_constant.list'))
cut_off = 0.75
gene_filter = ['B1N7B6','DKFZp686C15213']
set_gene = set()

def find_similar_from_a_list(str0,list_nonig,list_nonig_gene):
    out0 = False
    for n0 in range(0,len(list_nonig)):
        r0 = Levenshtein.ratio(str0,list_nonig[n0])
        if  r0 > cut_off and (not list_nonig_gene[n0] in gene_filter):
            set_gene.add(list_nonig_gene[n0])
            print(str0+','+list_nonig[n0]+','+list_nonig_gene[n0]+',Similar_score='+str(r0))
            out0 = True
    return out0
        
n_total = len(list_constant)
n0 = 0
for str0 in list_constant:
    if find_similar_from_a_list(str0, list_nonig, list_nonig_gene):
        n0 +=1
print('Percentage of Ig constant region peptides can be found homologous in Non-Ig peptides=')
print(str(n0))
print(str(float(n0)/n_total*100)+'%')
print(set_gene)
        
        
