import os
import fileinput

#given a dictonary, return set of the keys
def get_key(dict0):
    key_set = set()
    for key, _ in dict0.iteritems():
        key_set.add(key)
    return key_set

#given cell type name, return a class file line
def get_str(name0,header0):
    str0 = name0
    for element0 in header0:
        if element0 == name0:
            str0 = str0 + '\t' + '1'
        else:
            str0 = str0 + '\t' + '2'
    return str0+'\n'

  

#get all tab, csv and txt files in the current directory
files = [f for f in os.listdir('.') if os.path.isfile(f)]
f_use = []
for f in files:
    if '.csv' in f or '.txt' in f or '.tab' in f:
        f_use.append(f)
print f_use
#overall step
####create a dictionary for each file
#####determine what genes to include (interesection of genes in each file)
#####
######all gene names and values should be separatd by a tab aka '\t
######header is for later class file generating

#create dictionary
header0 = []
dict_list = []
n = 0
for f in f_use:
    n += 1
    print n
    dict0 = dict()
    head0 = True
    for line0 in open(f,'rU'):
        line0 = line0.rstrip()
        if head0:
            head0 = False
            #print line0
            print(f+' header length='+str(len(line0)))
            line0 = line0.split('\t')
            header0.extend(line0[1:])
        else:
            index0 = line0.find('\t')
            name0 = line0[0:index0].upper()
            #print name0
            val0 = line0[index0+1:]
            #print val0
            dict0[name0] = val0
    dict_list.append(dict0)
#####determine what genes to include (interesection of genes in each file)
gene_set = set()  
for dict0 in dict_list:
    gene0 = get_key(dict0)
    #print gene0
    if len(gene_set) > 1:
        gene_set = gene_set.intersection(gene0)
    else:
        gene_set = gene0
#print gene_set
###recreate the reference file from all dictionaries 
#initialize file
ref_f = open('reference_file.input','w+')
ref_f.write('Gene_symbols/Cell_types')
#write cell type
for name0 in header0:
    #print(name0)
    ref_f.write('\t'+name0)
ref_f.write('\n')
#write gene values
for gene0 in gene_set:
    ref_f.write(gene0+'\t')
    for dict0 in dict_list:
        if not dict0[gene0][-1] == '\t':
            ref_f.write(dict0[gene0]+'\t')
        else:
            ref_f.write(dict0[gene0])
    ref_f.write('\n')
ref_f.close()
print('Reference file is saved to reference_file.input')
####recreate the class files
class_f = open('class_file.input','w+')
#create class list
name_list = []
#print header0
for name0 in header0:
    if not name0 in name_list:
        name_list.append(name0)
        str0 = get_str(name0,header0)
        class_f.write(str0)
class_f.close()
print('Class file is saved to class_file.input')
            