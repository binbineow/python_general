#######this script takes in .codingchange file and output netmhcpan (class I) prediction
import fileinput
import subprocess
import os


####iedb_path should be in the bash path file, so the program can call NetMHCIIpan directly

#convert HLA-DRB1*04:01 into DRB1_0401 format
def convert_hla(str0):
    str1 = str0.split('*')[1].split(':')
    num1 = str1[0]
    num2 = str1[1]
    return 'DRB1_'+num1+num2

#given a list of strings, return a list of their lengths in order
def get_len_list(list1):
    list0 = []
    for x in list1:
        list0.append(len(x))
    return list0

def make_one_line(list0):
    str_out = ''
    for x in list0:
        str_out = str_out + x
    return str_out

def read_netmhc_xls(file0,list_run):
    file_out = open(file0,'r')
    dict0 = dict()
    for line0 in file_out:
        line0 = line0.rstrip()
        if 'Sequence' in line0:
            line0 = line0.split('\t')
        if line0[1] in list_run:
            dict0[line0[1]] = [float(line0[3]),float(line0[4]),float(line0[5])]
    return dict0

def remove_file(file0):
    cmd_line = 'rm '+file0
    cmd0 = subprocess.Popen(cmd_line, shell=True)
    cmd0.wait()     

def run_netmhciipan(hla_type_run,list_run,len_run):
    set_run = set(list_run)
    list_run = list(set_run)
    str_long = make_one_line(list_run)
    #create a pep file
    file_name_in = hla_type_run+str(len_run)+'.pep'
    file_pep = open(file_name_in,'w+')
    file_pep.write('>temp0\n')
    file_pep.write(str_long)
    file_pep.close()    
    #run
    cmd_line = 'netMHCIIpan -f '+file_name_in+ ' -inptype 0 -a '+ hla_type_run+ \
     ' -length '+str(len_run)+' >'+file_name_in+'.temp' + ' -xls -xlsfile '+file_name_in+'.xls ' + \
    '-tdir /home/stanford/rbaltman/users/bchen45/software/netMHCIIpan-3.1/tmp'
    print cmd_line
    cmd0 = subprocess.Popen(cmd_line,shell=True)      
    cmd0.wait() 
    #get data
    dict_out = read_netmhc_xls(file_name_in+'.xls', list_run)
    #remove temp file
    if os.path.isfile(file_name_in):
        remove_file(file_name_in)
        remove_file(file_name_in+'.temp')
        remove_file(file_name_in+'.xls')
    return dict_out

#main function
#give the hla_type in HLA_DRB1*04:01 format, and a list of sequences, return a dictionary
#The dictionary maps individual sequences -> [binding affinity, ranking score 
def predict_netmhciipan(hla_type0,list_seq):
    dict0 = dict()
    list_len0 = get_len_list(list_seq)
    set_len0 = set(list_len0)
    if '*' in hla_type0:
        hla_type_run = convert_hla(hla_type0)
    else:
        hla_type_run = hla_type0
    for len_run in set_len0:
        list_run = []
        for x in range(0,len(list_len0)):
            if list_len0[x] == len_run:
                list_run.append(list_seq[x])
        dict_run = run_netmhciipan(hla_type_run,list_run,len_run)
        print(dict_run)
        dict0.update(dict_run)
    return dict0