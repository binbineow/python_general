#######this script takes in .codingchange file and output netmhcpan (class I) prediction
import fileinput
import subprocess
import os


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
    #cmd_line = 'netMHCIIpan -f '+file_name_in+ ' -inptype 0 -a '+ hla_type_run+ \
    # ' -length '+str(len_run)+' > '+file_name_in+'.temp' + ' -xls -xlsfile '+file_name_in+'.xls ' + \
    #'-tdir /home/stanford/rbaltman/users/bchen45/software/netMHCIIpan-3.1/tmp'
    cmd_line = 'netMHCIIpan -f '+file_name_in+ ' -inptype 0 -a '+ hla_type_run+ ' >'+file_name_in+'.temp ' +\
    ' -length '+str(len_run)+ ' -xls -xlsfile '+file_name_in+'.xls ' + \
    '-tdir /share/PI/rbaltman/bchen45/software/IEDB/netMHCIIpan-3.1/tmp'
    print cmd_line
    #cmd_line_list = cmd_line.split(' ')
    #print cmd_line_list
    #subprocess.call(cmd_line_list)
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