from utilities import *

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
overlap_data = defaultdict(dumb)
gene1_set_list = []
gene2_set_list = []
MCL_data = pickle.load(open('MCL_data11_18_2015v0.9.dict','r'))
for x in MCL_data['pid']['pid']:
    gene1_set_list.append(set(MCL_data[x]['mutation']))
    gene2_set_list.append(set(MCL_data[x]['mutation']))
counter1 = count_set(gene1_set_list)
print counter1
counter2 = count_set(gene2_set_list)
print counter2


set1 = set()
set2 = set()
for x in MCL_data['pid']['pid']:
    set1 = set1 | set(MCL_data[x]['MHC1_gene'])
    set2 = set2 | set(MCL_data[x]['MHC2_gene'])
set1_share = pickle.load(open('MHC1_shared_MCL.set','r'))
set2_share = pickle.load(open('MHC2_shared_MCL.set','r'))
set_mut1 = pickle.load(open('filtered_mutation_MCL.set','r'))
set_mut2 = pickle.load(open('filtered_mutation_MCL.set','r'))

for n0 in range(0,len(MCL_data['pid']['pid'])+1):
    set1 = remove_element(set1,counter1,n0)
    overlap_data[n0]['MHC1_gene'] = list(set1)
    set2 = remove_element(set2,counter2,n0)
    overlap_data[n0]['MHC2_gene'] = list(set2)
    set1_share = remove_element(set1_share,counter1,n0)
    overlap_data[n0]['MHC1_share'] = list(set1_share)
    set2_share = remove_element(set2_share,counter2,n0)
    overlap_data[n0]['MHC2_share'] = list(set2_share)
    set_mut1 = remove_element(set_mut1,counter1,n0)
    overlap_data[n0]['mutation'] = list(set_mut1)
    set_mut2 = remove_element(set_mut2,counter2,n0)
    end0 = 'RPKM>1_Presence>'+str(n0)
    #plot_two_set_num(len(set_mut1)-len(set1_share),len(set1)-len(set1_share),len(set1_share),'Genes with mutations \nin all patients','Genes detected \nin MHCI Ligandome','plot/MHCI_All_patients_Venn_'+end0)
    #plot_two_set_num(len(set_mut2)-len(set2_share),len(set2)-len(set2_share),len(set2_share),'Genes with mutations \nin all patients','Genes detected \nin MHCII Ligandome','plot/MHCII_All_patients_Venn_'+end0)

pickle.dump(overlap_data,open('overlap_datav0.9.dict','w+'))

    
    
    
    

             