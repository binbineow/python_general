 #importing
import fileinput

#initial the dictionary
import cPickle as pickle
#dict_g_u =pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_gene_to_uni.dict'))
dict_u_g = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_uni_to_gene.dict'))
dict_cell = dict()
unconverted = 0
total = 0

def get_name(str):
    global unconverted
    if str in dict_u_g:
        name0 = dict_u_g[str]
    else:
        name0 = str
        unconverted += 1
    if '_' in name0:
        name0 = name0.split('_')[0]
        unconverted -= 1
    return name0

#main
length0 = 0
for line in fileinput.input():
    line = line.rstrip()
    line = line.split(',')
    if length0 == 0:
        length0 = len(line)
    if not line[0] == 'Checked':
        name0 = get_name(line[3])
        total += 1
        dict_cell[name0] = int(line[8+(len(line)-length0)])
        print(int(line[8]))

print float(unconvreted)/total

#save
#pickle.dump(dict_uni_to_gene,open('dict_uni_to_gene.dict','w+'))
#pickle.dump(dict_gene_to_uni,open('dict_gene_to_uni.dict','w+'))