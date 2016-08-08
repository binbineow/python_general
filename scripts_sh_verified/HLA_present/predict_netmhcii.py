#template for general purposes
import fileinput
import subprocess
import os
import string

#variable


def input0():
    ########input file path
    pathin=''
    #######output file path
    pathout = ''
    #########input file (codingchange)
    input_file = 'TCGA-09-2050-01A-01W-0799-08.codingchange'
    ###########IEDB directory
    iedb_path = '/home/alizadehlab/bchen45/IEDB/mhc_i/src/'
    hla_type = 'HLA-A*01:01'
    len_list = [8,9,10]
    outfile_name = input_file[0:16]+hla_type.translate(string.maketrans("",""), string.punctuation)
    cut_off = 200
    return [pathin,pathout,input_file,iedb_path,hla_type,len_list,outfile_name,cut_off]
    
    
def write_seq(line_str,anchor1,anchor2,len0,file_seq):
    #print anchor1
    #print len0
    #print anchor2
    str0 = line_str[max(anchor1-len0+1,0):anchor2+len0]
    file_seq.write(str0)
    return str0

def pep_predict(str,len0):
    set0 = set()
    len_str= len(str)
    for i in range(0,len0):
        if len_str - i < len0:
            break
        set0.add(str[i:i+len0])
    #print set0
    return set0    

def pep_filter(file_name0,cut_off,set_pre):
    #print set_pre
    file_in = open('temp2'+pathout+file_name0,'r')
    file_out = open(pathout+file_name0+'output','w+')
    file_out2 = open(pathout+file_name0+'summary','w+')
    num0 = 0
    for line in file_in:
        if not 'length' in line:
            if 'Consensus' in line:
                c0 = 8
            else:
                c0 = 14
            line_s = line.rstrip().split('\t')
            line_s = filter(None,line_s)
            #print line_s
            if line_s[5] in set_pre:
                file_out.write(line)
                #print line_s[5]
                if float(line_s[c0])<=cut_off:
                    file_out2.write(line_s[5]+'\n')
                    num0 = num0+1
        else:
            file_out.write(line)
    file_out2.write('#'+str(num0)+'\t'+str(cut_off)+'\n')
    file_in.close()
    file_out.close()
    file_out2.close()
            
    
    

def main(pathin,pathout,input_file,iedb_path,hla_type,len_list,outfile_name,cut_off): 
    #process the files and output binding prediction for each length
    file_check = True
    if os.stat(pathin+input_file).st_size < 1:
        file_check = False
        print pathin+input_file+" is empty"
    for len0 in len_list:
        if not file_check:
            break
        #read in the mutation string and create mutation input file
        total_line = ''   
        line_str = ''
        read0 = False
        file_in = open(pathin+input_file,'r')
        file_seq = open(pathout+'temp'+pathout+outfile_name+'_len'+str(len0),'w+')
        file_seq.write('>test\n')        
        set_pre = set()
        for line in file_in:
            line = line.rstrip() 
            if len(line) <1:
                break       
            if line[0] == '>' and not 'WILDTYPE' in line:
                if len(line_str) >1 :
                    line_str = line_str[0:-1]
                    seq0 = write_seq(line_str,anchor1,anchor2,len0,file_seq)
                    total_line = total_line + seq0
                    #predict the final peptide list
                    #print pep_predict(seq0,len0)
                    set_pre = set_pre | pep_predict(seq0,len0)      
                    line_str = ''
                if not ' synonymous' in line:
                    #old version of annovare
                    #anchor1 = line.split('-')[0].split(' ')[-1]
                    #anchor2 = line.split('-')[1].split(' ')[0]   
                    #new version of annovare        
                    anchor1 = int(line.split('position ')[1].split(' changed')[0]) 
                    anchor2 = anchor1   
                    read0 = True
            elif line[0] == '>' and 'WILDTYPE' in line:
                read0 = False
            else:
                if read0:
                    line_str = line_str+line
        if len(line_str) >1 :
            line_str = line_str[0:-1]
            seq0 = write_seq(line_str,anchor1,anchor2,len0,file_seq)
            total_line = total_line + seq0
            #predict the final peptide list
            set_pre = set_pre | pep_predict(seq0,len0)
            
        #predict with IEDB
        #python $IEDBPath/predict_binding.py IEDB_recommended HLA-A*01:01 9 test0.txt
        file_seq.close()
        cmd = 'python '+iedb_path+'predict_binding.py IEDB_recommended '+hla_type+' '+str(len0)+' '+  \
               pathout+'temp'+outfile_name+'_len'+str(len0)+ ' > '+ 'temp2'+pathout+outfile_name+'_len'+str(len0)
        print cmd
        cmd0 = subprocess.Popen(cmd,shell=True)      
        cmd0.wait() 
        #filter with
        pep_filter(pathout+outfile_name+'_len'+str(len0),cut_off,set_pre)
        #close file
        file_in.close()


if __name__ == '__main__':
    pass
[pathin,pathout,input_file,iedb_path,hla_type,len_list,outfile_name,cut_off] = input0()
main(pathin,pathout,input_file,iedb_path,hla_type,len_list,outfile_name,cut_off)

