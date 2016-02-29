 #importing
import fileinput
import cPickle as pickle

#initial the dictionary
dict_uni_to_gene = dict()
dict_gene_to_uni = dict()

#main
for line in fileinput.input():
    line = line.rstrip()
    line = line.split('\t')
    if line[1] == 'GeneCards':
        dict_uni_to_gene[line[0]] =line[2]
        dict_gene_to_uni[line[2]] =line[0]

#save
pickle.dump(dict_uni_to_gene,open('dict_uni_to_gene.dict','w+'))
pickle.dump(dict_gene_to_uni,open('dict_gene_to_uni.dict','w+')) 
 
