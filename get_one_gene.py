import fileinput
import cPickle as pickle
line = ''
for line0 in fileinput.input():
    if not line0[0] == '<' or not line0[0] == '>':
        line0 = line0.rstrip()
        line = line+line0
pickle.dump(line,open('human_proteinome_oneline.str','w+'))
