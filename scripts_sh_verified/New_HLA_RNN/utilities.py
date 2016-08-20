#update 11/23/2015
import Levenshtein
import os
import matplotlib.pyplot as plt
import cPickle as pickle
import fileinput
from collections import Counter
from collections import defaultdict
from matplotlib import pyplot as plt
import numpy as np
import subprocess

def dumb(): return defaultdict(list)

def touch_file(path_file0):
    cmd = 'touch '+path_file0
    cmd0 = subprocess.Popen(cmd,shell=True)      
    cmd0.wait() 

#write_list take in a file variable, list of strings and and separating symbol
#each line contains all strings in the list separated by the symbol, ended with '\n'
def write_list(file_out,list1,symbol0):
    for x in list1:
        file_out.write(x+symbol0)
    file_out.write('\n')

def remove_merged(list0):
    for i in range(0,len(list0)):
        if ',' in list0[i]:
            list0[i] = list0[i].split(',')[0]
        if ';' in list0[i]:
            list0[i] = list0[i].split(';')[0]
        if 'merged_' in list0[i]:
            list0[i] = list0[i].split('_')[1]
    return list0

def remove_low_RNA(set0,dictRNA,RNA_cut_off):
    set_remove = set()
    #global set_unmapped
    for x in set0:
        if x in dictRNA:
            if float(dictRNA[x]) < RNA_cut_off:
                set_remove.add(x)
        else:
            #set_unmapped.add(x)
            set_remove.add(x)
    set0 -= set_remove
    return set0
    
def get_counter(gene_list0):
    counter2 = Counter()
    for x in gene_list0:
        counter2[x] += 1
    return counter2

def plot_counter2(gene_list0,filename0,title0,y0,x0):
    
    plt.hist(gene_list0)
    plt.title(title0)
    plt.xlabel(x0)
    plt.ylabel(y0)
    save(filename0,'png')
  
    
def plot_counter(gene_list0,filename0,title0,y0,x0,print_s,max0):
    counter2 = Counter()
    for x in gene_list0:
        counter2[x] += 1
    #key_list = []
    value_list = []
    for key,value in counter2.iteritems():
        #key_list.append(key)
        if value<max0:
            value_list.append(value)
    plt.hist(value_list)
    plt.title(title0)
    plt.xlabel(x0)
    plt.ylabel(y0)
    save(filename0,'png')
    if print_s:
        print counter2.most_common(10)
    
    

def get_gene_from_fasta(x,dict0):
    if len(x.split('|'))==3:
        if x.split('|')[1] in dict0:
            list_return = dict0[x.split('|')[1]]
        else:
            list_return = ('No_gene')
    else:
        list_return=('No_gene')
    return list_return

def get_gene_from_ms(list0,dict0,pid=''):
    list_return = []
    for x in list0:
        if 'idiotype' in x and pid in x:
            list_return.append('idiotype')
        elif len(x.split('|'))==3:
            list_return.append(x.split('|')[1])
                #print x
        elif len(x.split('|'))==4:
            if x.split('|')[1] in dict0:
                list_return.append(dict0[x.split('|')[1]])
            else:
                list_return.append(x.split('|')[2])
        else:
            list_return.append('No_gene')
    return list_return

def make_any_dict(dictName,pickle_s,label0,header_s):
    dict0 = dict()
    for line0 in fileinput.input():
        if header_s:
            header_s = False
        else:
            line0 = line0.rstrip()
            line0 = line0.split(label0)
            dict0[line0[0]] = line0[1]
    if pickle_s:
        pickle.dump(dict0,open(dictName,'w+'))
    else:
        return dict0

def read_col(filename,deli0,col_num,head_s):
    #file_name should include path
    #deli0 = deliminator
    #col_num starts with 03
    # head_s = True when there is a header (one line only)
    #return a list
    file0 = open(filename,'r')
    list_out = []
    for line0 in file0:
        if head_s:
            head_s = False
        else:
            line0 = line0.rstrip().split(deli0)
            list_out.append(line0[col_num])
    return list_out

def similar_comp(str1,str2):
    str0 = 'different'
    dis0 = 0
    if str1 == '' or str2 == '':
        print 'One of the strings is empty'
    if str1 == str2:
        str0 = 'str1=str2'
    elif str1 in str2:
        str0 = 'str1-in-str2'
    elif str2 in str1:
        str0 = 'str2-in-str1'
    else:
        dis0 = Levenshtein.distance(str1,str2)
    return [str0,dis0]

def save(path, ext='png', close=True, verbose=True):
    """Save a figure from pyplot.

    Parameters
    ----------
    path : string
        The path (and filename, without the extension) to save the
        figure to.

    ext : string (default='png')
        The file extension. This must be supported by the active
        matplotlib backend (see matplotlib.backends module).  Most
        backends support 'png', 'pdf', 'ps', 'eps', and 'svg'.

    close : boolean (default=True)
        Whether to close the figure after saving.  If you want to save
        the figure multiple times (e.g., to multiple formats), you
        should NOT close it in between saves or you will have to
        re-plot it.

    verbose : boolean (default=True)
        Whether to print information about when and where the image
        has been saved.

    """
    
    # Extract the directory and filename from the given path
    directory = os.path.split(path)[0]
    filename = "%s.%s" % (os.path.split(path)[1], ext)
    if directory == '':
        directory = '.'
 
    # If the directory does not exist, create it
    if not os.path.exists(directory):
        os.makedirs(directory)
 
    # The final path to save to
    savepath = os.path.join(directory, filename)
 
    if verbose:
        print("Saving figure to '%s'..." % savepath),
 
    # Actually save the figure
    plt.savefig(savepath)
    
    # Close it
    if close:
        plt.close()
 
    if verbose:
        print("Done")
######################################
        