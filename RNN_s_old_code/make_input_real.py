import fileinput
import random
import Levenshtein


pos = []
neg = []
for i,line in enumerate(fileinput.input()):
    if i == 0: continue
    line= line.rstrip().split(',')
    #'|SPA|' to be excluded
    if not '|SPA|' in line[2]:
	for _ in range(0,3):
       		neg0 = ''.join(random.sample(line[0],len(line[0])))
        	if not neg0 == line[0]:
            		neg.append(neg0+'\t'+'0')
        pos.append(line[0]+'\t'+'1')
    #5 sequence 6value

for p in pos + neg:
    print(p)
