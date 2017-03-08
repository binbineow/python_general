#from utilities import *
dictRNA_file = '/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCLRNASeq_ave.dict'
import math
import cPickle as pickle
import numpy as np
import matplotlib.pyplot as plt
def dumb(): return defaultdict(list)

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
    plt.title(title)
    save(title)
    
#cal jaccard 
# only considering substring not additive distances
def count_shared(listx, listy):
    n_shared = 0
    for frag0 in listx:
        for fragy in listy:
            if frag0 in fragy or fragy in frag0:
            #if frag0 in fragy or fragy in frag0:
                n_shared += 1
                break
    return n_shared

# calculate jaccard, union is the traditional union, intersection allows sub string match
def get_jaccard(listx,listy):
    set_union = set(listx) | set(listy)
    val_j = float(count_shared(listx,listy))
    
    

#boostrap jaccard 



overlap_data = defaultdict(dumb)
gene1_set_list = []
gene2_set_list = []
MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict','r'))
for x in MCL_data['pid']['pid']:
    gene1_set_list.append(set(MCL_data[x]['MC1_gene']))
    gene2_set_list.append(set(MCL_data[x]['MC2_gene']))
counter1 = count_set(gene1_set_list)
#print counter1
counter2 = count_set(gene2_set_list)
RNAdict = pickle.load(open(dictRNA_file,'r'))
#print counter2
x0 = []
y0 = []
for key, value in counter1.iteritems():
    if key in RNAdict: 
        y0.append(value)
        x0.append(math.log(RNAdict[key]))
plot_scatter(x0,y0,'Estimated RNA Expression (log)','Frequency of Gene in MHC1 Ligandome','RNA expression profiles of genes in MHC1 Ligandome')

x0 = []
y0 = []
for key, value in counter2.iteritems():
    if key in RNAdict: 
        y0.append(value)
        x0.append(math.log(RNAdict[key]))
plot_scatter(x0,y0,'Estimated RNA Expression (log)','Frequency of Gene in MHC2 Ligandome','RNA expression profiles of genes in MHC2 Ligandome')
    

'''
set1 = set()
set2 = set()
for x in MCL_data['pid']['pid']:
    set1 = set1 | set(MCL_data[x]['MHC1_gene'])
    set2 = set2 | set(MCL_data[x]['MHC2_gene'])
set1_share = pickle.load(open('MHC1_shared_MCL.set','r'))
set2_share = pickle.load(open('MHC2_shared_MCL.set','r'))
set_mut1 = pickle.load(open('filtered_mutation_MCL.set','r'))
set_mut2 = pickle.load(open('filtered_mutation_MCL.set','r'))
'''

##########heatmap generating


pickle.dump()

    
import cPickle as pickle
import plotly.plotly as py
import plotly.graph_objs as go
path0 = '/Users/binbineow2/Documents/Machine_Learning/HLA_prediction/Figure_MCL/heatmap/'
pid_list = pickle.load(open(path0+'heatmap_pid.list','r'))
z = pickle.load(open(path0+'heatmap_mhc2_J.mat','r'))
#z = pickle.load(open(path0+'heatmap_drb.mat','r'))
for n0 in range(0,len(z)):
    z[n0] = list(reversed(z[n0]))
    
x = pid_list
x = list(reversed(pid_list))
y = pid_list



annotations = []
for n, row in enumerate(z):
    for m, val in enumerate(row):
        var = z[m][n]
        annotations.append(
            dict(
                text=str(val),
                x=x[m], y=y[n],
                xref='x1', yref='y1',
                font=dict(color='white' if val > 0.5 else 'black'),
                showarrow=False)
            )

colorscale = [[0, '#CCE5FF'], [1, '#003366']]  # custom colorscale
#trace = go.Heatmap(x=x, y=y, z=z, colorscale=colorscale, showscale=False)
trace = go.Heatmap(x=x, y=y, z=z, colorscale=colorscale)
fig = go.Figure(data=[trace])
fig['layout'].update(
    title="Overlapping of peptides between patients (Jaccard Index)",
    #title="HLA-DRB1 Sequence similarity between patients (Levenshtein Ratio)",
    #annotations=annotations,
    xaxis=dict(ticks='', side='top'),
    # ticksuffix is a workaround to add a bit of padding
    yaxis=dict(ticks='', ticksuffix='  '),
    width=900,
    height=900,
    autosize=True
)

url = py.plot(fig, filename='Peptide overlap Heatmap', height=1000)
  
    
    

             