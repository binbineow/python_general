from utilities import *
import string

def get_HLA_name(str0):
    print str0
    if ':' in str0:
        str1 = str0.split(':')[0]+':'+str0.split(':')[1]
    else:
        str1 = str0
    return str1

def get_89_motif(str0):
    if len(str0) < 91:
        return str0
    else:
        star0 = string.find(str0,'RFL')
        str0 = str0[star0:star0+89]
        return str0

dict0 = dict()
allele_name = ''
for line0 in fileinput.input():
    line0 = line0.rstrip()
    if '>HLA' in line0 in line0:
        if 'DRB1' in allele_name:
            if not allele_name in dict0:
                dict0[allele_name] = get_89_motif(str0_local)
        allele_name = get_HLA_name(line0.split(' ')[1])
        str0_local = ''
    else:
        str0_local = str0_local + line0
print dict0
pickle.dump(dict0,open('DRB1_seq.dict','w+'))