import fileinput
target_mhc = 'HLA-DRB1*04:01'

#1 ligand ID; 8 pubmedID; 23 sequence; 101 Assay; 109 result category; 111 EC50; 127 MHC type
list_target = []
for line0 in fileinput.input():
    line0=line0.rstrip()
    line0=line0.split('"')
    list0 = []
    if len(line0) > 127:
        print line0[127]
        if line0[127] == target_mhc:
            list0 = [line0[1],line0[8],line[23],line[101],line[109],line[111],line[127]]
        if not list0 == []:
            list_target.append(list0)
for list0 in list_target:
    line0 = ''
    for ele0 in list0:
        line0=line0+','+ele0
    print line0[1:]