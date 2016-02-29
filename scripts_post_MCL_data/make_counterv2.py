from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/jeko_pC.dict'
import math

def plot_hist(gene_set_list):
    counter = defaultdict(int)
    for s in gene_set_list:
        for el in s:
            counter[el] += 1
    #print(counter)
    counter2 = defaultdict(int)
    for key,value in counter.iteritems():
        counter2[value] += 1
    for key,value in counter2.iteritems():
        print str(key)
    for key,value in counter2.iteritems():
        print str(value)
        
def count_set(gene_set_list):
    counter = defaultdict(int)
    for s in gene_set_list:
        for el in s:
            counter[el] += 1
    return counter

def remove_element(set1,dict0,cut_off0):
    set_remove = set()
    for x in set1:
        if dict0[x]<cut_off0:
            set_remove.add(x)
    set1 -= set_remove
    return set1

def plot_scatter(x,y,x_name,y_name,title):
    plt.figure()
    plt.plot(x,y,'bo')
    plt.ylabel(y_name)
    plt.xlabel(x_name)
    plt.ylim([0,max(y)+1])
    plt.title(title)
    save(title)

overlap_data = defaultdict(dumb)
gene1_set_list = []
gene2_set_list = []
MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict','r'))
for x in MCL_data['pid']['pid']:
    print x
    print 'MHC1'
    for i in range(0,len(MCL_data[x]['MHC1_gene'])):
        gene0 = MCL_data[x]['MHC1_gene'][i]
        if gene0 == 'idiotype':
            print MCL_data[x]['MHC1_frag'][i]
    print 'MHC2'
    for i in range(0,len(MCL_data[x]['MHC2_gene'])):
        gene0 = MCL_data[x]['MHC2_gene'][i]
        if gene0 == 'idiotype':
            print MCL_data[x]['MHC2_frag'][i]
    gene1_set_list.append(set(MCL_data[x]['MHC1_gene']))
    gene2_set_list.append(set(MCL_data[x]['MHC2_gene']))
counter1 = count_set(gene1_set_list)
#print counter1
counter2 = count_set(gene2_set_list)
RNAdict = pickle.load(open(dictRNA_file,'r'))
pickle.dump(counter1,open('MCL_MHC1_gene_frequency.count','w+'))
pickle.dump(counter2,open('MCL_MHC2_gene_frequency.count','w+'))
#print counter2
#get scatter plot

x0 = []
y0 = []
for key, value in counter1.iteritems():
    if key in RNAdict: 
        y0.append(value)
        x0.append(math.log(float(RNAdict[key])+0.000001))
plot_scatter(x0,y0,'Estimated RNA Expression (log)','Frequency of Gene in MHC1 Ligandome','RNA expression profiles of genes in MHC1 Ligandome')

x0 = []
y0 = []
for key, value in counter2.iteritems():
    if key in RNAdict: 
        y0.append(value)
        x0.append(math.log(float(RNAdict[key])+0.000001))
plot_scatter(x0,y0,'Estimated RNA Expression (log)','Frequency of Gene in MHC2 Ligandome','RNA expression profiles of genes in MHC2 Ligandome')

#pickle.dump(overlap_data,open('overlap_datav0.9.dict','w+'))
 
#get high frequency genes
'''
print('MHC1 >=17')
for key, value in counter1.iteritems():
    if value>= 17:
        print key
print('MHC2 >=16')
for key, value in counter1.iteritems():
    if value>= 16:
        print key
'''

    
    

             