import fileinput
#import random
#import Levenshtein
#import cPickle as pickle

ration = 4
one_gene_path = '/scratch/users/bchen45/HLA_prediction/IEDB/test0/human_proteinome_oneline.str'
#onegenestr = pickle.load(open(one_gene_path,'r'))
#len_one = len(onegenestr)
pos = []
for i,line in enumerate(fileinput.input()):
	if i == 0: continue
	line= line.rstrip().split(',')
	if not '|SPA|' in line[2]:
		print(line[0]+'\t2')

