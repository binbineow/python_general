infile_name = 'NimbleGenChart_for_analysis.csv'
outfile_name = 'NimbleGenChart_for_analysisv2.csv'
dict_name = 'NimbleGen_substitution_pep.csv'
path0 = '/Users/binbineow2/Documents/Machine_Learning/RNN_model/NimbleGen/'
file_out = open(path0+outfile_name,'w+')
dict0 = dict()
#make the dictionary 
for line0 in open(path0+dict_name,'rU'):
    line0 = line0.rstrip()
    line0 = line0.split(',')
    dict0[line0[0]] = line0[1]

pep_set = set()
#adding the sub_info
for line0 in open(path0+infile_name,'rU'):
    line0 = line0.rstrip()
    line1 = line0.split(',')
    out_line = ''
    if not line1[0] == 'Sequence':
        pep0 = line1[0]
        if pep0 in dict0:
            if dict0[pep0] == '1':
                out_line = out_line+',Substitution_low_priority'
            elif 'gative' in dict0[pep0]:
                out_line = out_line+',Substitution_negative_control'
            else:
                out_line = out_line+',Substitution'
        else:
            out_line = out_line+',Single'
        if pep0 in pep_set:
            out_line = out_line+',Repeat'
        else:
            pep_set.add(pep0)
            out_line = out_line+',First_presence'
    else:
        out_line = out_line  + ',Array_type,Notes'
    file_out.write(line0+out_line+'\n')
file_out.close()
            