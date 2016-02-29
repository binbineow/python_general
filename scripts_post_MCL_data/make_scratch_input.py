from utilities import *
import random
MCL_data = pickle.load(open('MCL_data11_18_2015v0.9.dict'))
pid = 'MCL041'
one_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/human_proteinome_oneline.str'
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
path0 = 'secondary_prediction/'
set_1 = []
set_2 = []
set_r = []
for type0 in ['MHC1','MHC2']:
    file_random = path0+type0+'random_protein.fasta'
    file_random0 = open(file_random,'w+')
    name0 = type0+'_frag'
    file_MHC1_name = path0+ pid+name0+'.fasta'
    file_out = open(file_MHC1_name,'w+')
    set0 = set(MCL_data[pid][name0])
    for x in set0:
        file_out.write('>pid\n'+x+'\n')
        rand0 = random.randint(0,len_one)
        neg0 = onegenestr[rand0:rand0+len(x)]
        file_random0.write('>pid\n'+neg0+'\n')
    file_out.close()
    file_random0.close()


        
    
    
