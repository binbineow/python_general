import fileinput
import random
import cPickle as pickle

pos = []
neg = []
for i,line in enumerate(fileinput.input()):
    if i == 0: continue
    line= line.rstrip().split('\t')
    #5 sequence 6value
    if float(line[6])<50:
        pos.append(line[5]+'\t'+'1')
    if float(line[6])>5000 and float(line[6])<7000:
        neg.append(line[5]+'\t'+'0')
    if float(line[6]) >= 7000: break

pickle.dump(neg,open('neg.list','w+'))
pos_sam = random.sample(pos,5000)
neg_sam = random.sample(neg,5000)

for p in pos_sam + neg_sam:
    print(p)