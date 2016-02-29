import fileinput 
from collections import Counter

total_C = Counter()
common_each_C = Counter()
h_list = []
e_list = []
c_list = []
for line0 in fileinput.input():
    if not line0[0] == '>':
        local_C = Counter()
        line0 = line0.rstrip()
        list0 = list(line0)
        for x in list0:
            total_C[x] +=1
            local_C[x] +=1
        x0 = local_C.most_common(1)[0][0]
        
        common_each_C[x0] +=1
print total_C
total_f = 0
f_list = []
for _,value in total_C.iteritems():
    total_f = total_f+value
    f_list.append(value)
for value0 in f_list:
    print float(value0)/total_f
print common_each_C
total_f = 0
f_list = []
for _,value in common_each_C.iteritems():
    total_f = total_f+value
    f_list.append(value)
for value0 in f_list:
    print float(value0)/total_f
        
        