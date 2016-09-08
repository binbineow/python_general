 #importing
import fileinput
#initial the dictionary
import cPickle as pickle
#dict_g_u
#=pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_gene_to_uni.dict'))
#dict_u_g =a
#pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_uni_to_gene.dict'))
def get_pep_str(file0):
    list_out = []
    for line0 in open(file0,'r'):
        if 'High' in line0:
            pep0 = line0.split(',')[2]
            list_out.append(pep0)
    return list_out

def one_string_list(list0,sep0):
    str0 = ''
    for x in list0:
        str0 = str0+sep0+x
    return str0


def count_exclude(pep0,exclude0):
    import re
    num_out = 0
    for x in exclude0:
        num_out += len(re.findall('N'+x+'S',pep0))
        num_out += len(re.findall('N'+x+'T',pep0))
    return num_out

def get_weighted_ratio(file0,val_def,exclude0):
    import re
    num_total = 0.0
    aa_total = 0
    num_motif = 0.0
    aa_motif = 0
    for line0 in open(file0,'r'):
        if 'High' in line0:
            pep0 = line0.split(',')[2]
            val0 = line0.split(',')[-5]
            if len(val0) > 1:
                val0 = float(val0)
            else:
                val0 = val_def
            motifs_run = re.findall('N.T',pep0)+re.findall('N.S',pep0) 
            motif_p_run = re.findall('NPT',pep0)+re.findall('NPS',pep0)
            motif_num = len(motifs_run) - len(motif_p_run)
            
            ## counting weighted motif frequency
            val_motif = float(motif_num)*val0
            num_motif = num_motif + val_motif
            num_total = num_total + len(pep0)*val0
            ###only counting AA
            aa_total = aa_total + len(pep0)
            aa_motif = motif_num+aa_motif
    return [float(num_motif)/num_total,aa_motif,aa_total]
            
            