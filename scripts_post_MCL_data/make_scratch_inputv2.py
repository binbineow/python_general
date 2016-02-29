from utilities import *
import random
MCL_data = pickle.load(open('MCL_data11_18_2015v1.1.dict'))
one_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/human_proteinome_oneline.str'
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
path0 = 'secondary_prediction/'
set_1 = set()
set_2 = set()
set_r = set()
for type0 in ['MHC1','MHC2']:
    for pid0 in MCL_data['pid']['pid']:
        set_1 = set_1 | set(MCL_data[pid0]['MHC1_frag'])
        set_2 = set_2 | set(MCL_data[pid0]['MHC2_frag'])
#remove any substrings
l = list(set_1)
list_1 = [x for i,x in enumerate(l) if all([(i == i2) or (not (x in y)) for i2,y in enumerate(l)])]       
l = list(set_2)
list_2 = [x for i,x in enumerate(l) if all([(i == i2) or (not (x in y)) for i2,y in enumerate(l)])] 
#list_1 = list(set_1)
#list_2 = list(set_2)
print 'MHC1_num='+str(len(list_1))
print 'MHC2_num='+str(len(list_2))

#create a random sequence control
#mhc1
file_out_1 = path0+'mhc_1_input.fasta'
file_mhc1 = open(file_out_1,'w+')
file_out_1r = path0+'random_1_input.fasta'
file_random1 = open(file_out_1r,'w+')
for x0 in list_1:
    file_mhc1.write('>'+x0+'\n'+x0+'\n')
    rand0 = random.randint(0,len_one)
    rstr0 = onegenestr[rand0:rand0+len(x0)]
    file_random1.write('>'+rstr0+'\n'+rstr0+'\n')
file_mhc1.close()
file_random1.close()    
#mhc2
file_out_2 = path0+'mhc_2_input.fasta'
file_mhc2 = open(file_out_2,'w+')
file_out_2r = path0+'random_2_input.fasta'
file_random2 = open(file_out_2r,'w+')
for x0 in list_2:
    file_mhc2.write('>'+x0+'\n'+x0+'\n')
    rand0 = random.randint(0,len_one)
    rstr0 = onegenestr[rand0:rand0+len(x0)]
    file_random2.write('>'+rstr0+'\n'+rstr0+'\n')
file_mhc2.close()
file_random2.close()    
        
        

    
