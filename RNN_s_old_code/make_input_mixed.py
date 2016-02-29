import fileinput
import random
#import Levenshtein
import cPickle as pickle

ration = 2
one_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/human_proteinome_oneline.str'
onegenestr = pickle.load(open(one_gene_path,'r'))
len_one = len(onegenestr)
pos = []
neg = []
for i,line in enumerate(fileinput.input()):
    if i == 0: continue
    line= line.rstrip().split(',')
    #'|SPA|' to be excluded
    if not '|SPA|' in line[2]:
        #neg0 = ''.join(random.sample(line[0],len(line[0])))
        for _ in range(0,ration):
	    neg0 = ''.join(random.sample(line[0],len(line[0])))
            if not neg0 == line[0]:
 		neg.append(neg0+'\t'+'0')
            rand0 = random.randint(0,len_one)
            neg0 = onegenestr[rand0:rand0+len(line[0])]
            if not neg0 == line[0]:
                neg.append(neg0+'\t'+'0')
        pos.append(line[0]+'\t'+'1')
    #5 sequence 6value

for p in pos + neg:
    print(p)
