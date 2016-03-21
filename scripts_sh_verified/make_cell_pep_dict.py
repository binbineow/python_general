 #importing
import fileinput

#initial the dictionary
import cPickle as pickle
#dict_g_u =pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_gene_to_uni.dict'))
#dict_u_g = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_uni_to_gene.dict'))
path_gene_analysis = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/'
dict_cell = dict()


def get_name(str0,str1):
    if '_' in str0:
        name0 = str0.split('_')[0]
        return(name0,name0,name0)
    else:
        name0 = str0
        name1 = str1.split('|')[0]
        name2 = str1.split('|')[3].split('_')[0]
        return(name0,name1,name2)

#main
for line in fileinput.input():
    line = line.rstrip()
    line = line.split(',')
    #if length0 == 0:
    #    length0 = len(line)
    if (not line[0] == 'Checked') and (line[1] == 'High'):
        [name0,name1,name2] = get_name(line[3],line[4])
        if float(line[-3]) < 300:
            dict_cell[name0] = float(line[-3])
            dict_cell[name1] = float(line[-3])
            dict_cell[name2] = float(line[-3])


#save
print dict_cell
pickle.dump(dict_cell,open(path_gene_analysis+'jeko_pep_countv3.dict','w+'))
#pickle.dump(dict_gene_to_uni,open('dict_gene_to_uni.dict','w+'))