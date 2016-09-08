import fileinput
from __builtin__ import True
target_mhc = 'HLA-DRB1*15:01'

#1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
list_target = []
n_positive = 0 #positive but not positive-low
n_negative = 0 #negative
for line0 in fileinput.input():
    line0=line0.rstrip()
    line0=line0.split('"')
    list0 = []
    if len(line0) > 127:
        if len(line0[127])>1 and line0[127][0] != 'H':
            print line0[127]
        if line0[127] == target_mhc:
            list0 = [line0[1],line0[8],line0[23],line0[101],line0[109],line0[111],line0[127]]
            if 'Negative' in line0[109]:
                n_negative += 1
            if ('Positive' in line0[109]) and (not 'Low' in line0[109]):
                n_positive += 1
        if not list0 == []:
            list_target.append(list0)
print('number of'+target_mhc+'='+str(len(list_target)))
print('n_positive='+str(n_positive))
print('n_negative='+str(n_negative))
for list0 in list_target:
    line0 = ''
    for ele0 in list0:
        line0=line0+','+ele0
    #print line0[1:]

#1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type    
#this code will grab peptide sequence based on desired HLA types and bidning properties 
#return list of peptides and AA numbers of potential N-glycosylation motif and total AAs
def check_hla(hla_list,entry0):
    b_out = False
    for x in hla_list:
        if x in entry0:
            b_out = True
            break
    return b_out

def count_exclude(pep0,exclude0):
    import re
    num_out = 0
    for x in exclude0:
        num_out += len(re.findall('N'+x+'S',pep0))
        num_out += len(re.findall('N'+x+'T',pep0))
    return num_out
        

def grab_iedb_list(file0,hla_list,category0,exclude0):
    import re
    list_out = []
    motif_aa = 0
    total_aa = 0
    for line0 in open(file0,'r'):
        line0 = line0.rstrip()
        line0=line0.split('"')
        if len(line0) > 127:
            if check_hla(hla_list,line0[127]) and category0 in line0[109]:
                pep0 = line0[23]
                list_out.append(pep0)
                motifs_run = re.findall('N.T',pep0)+re.findall('N.S',pep0) 
                motif_p_run = count_exclude(pep0,exclude0)
                motif_num = len(motifs_run) - motif_p_run
                motif_aa += motif_num
                total_aa += len(pep0)
    return [list_out,motif_aa, total_aa]
                
        