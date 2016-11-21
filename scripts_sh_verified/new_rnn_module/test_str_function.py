import random
import cPickle as pickle

#working on xstream only if no one_gene_path given
#
def make_random_pep(list0,one_gene_path = '/home/stanford/rbaltman/users/bchen45/data/protein_general/human_proteinome_oneline.str'):
    import random
    list_out = []
    onegenestr = pickle.load(open(one_gene_path,'r'))
    len_one = len(onegenestr)
    for pos0 in list0:
        rand0 = random.randint(0,len_one)
        neg0 = onegenestr[rand0:rand0+len(pos0)]
        list_out.append(neg0)
    return list_out

def make_list_overlap_str(str0,str1):
    list_out = [str0]
    str_last = str0
    for x in str1:
        str_last = str_last + x
        list_out.append(str_last)
    print list_out
    return list_out

def cal_sensi(list0,cut_off):
    n0 = 0
    for x in list0:
        if x>= cut_off:
            n0 +=1
    return float(n0)/len(list0)

#break a X length long string into first half, middle half, and middle half 16
def break2_3halfs(list0,upper0 = 50):
    list_out = []
    n0 = 0
    for x in list0:
        if (not len(x) == 30) and len(x) >= upper0:
            n0 +=1
        if len(x) <= upper0:
            half_len = len(x)/2
            str0 = x[0:half_len+1]
            str1 = x[half_len-1:]
            str2 = x[half_len-8:half_len+8]
            list_out.extend([str0,str1,str2])
    print('Number of strings longer than '+str(upper0)+' = '+str(n0))
    return list_out
            
#break a X length long string into first half, middle half, and middle half 16
def walking_pep(list0,upper0 = 50,walk0=16):
    list_out = []
    n0 = 0
    for x in list0:
        if (not len(x) == 30) and len(x) >= upper0:
            n0 +=1
        if len(x) <= upper0:
            for i in range(0,len(x)-walk0):
                list_out.append(x[i:i+walk0])
    print('Number of strings longer than '+str(upper0)+' = '+str(n0))
    return list_out  
      
######################################

#get a list from a file name###        
def get_list_from_file(file_name0):
    file0 = open(file_name0,'r')
    list_out = []
    for x in file0:
        x = x.rstrip()
        list_out.append(x)
    return list_out

#remove str0
def remove_str0_from_a_list(list0,str0):
    list_out = []
    for x0 in list0:
        if not str0 in x0:
            list_out.append(x0)
    return(list_out)

#limit length
def limit_len_from_a_list(list0,len0):
    list_out = []
    for x0 in list0:
        if len(x0) <= len0:
            list_out.append(x0)
    return(list_out)