#template for general purposes
import fileinput
#import subprocess

#variable
pathin=''
pathout = ''
file_in_name = ''
file_out_name = ''
des_list = []
def main():
    #print '>Tophead'
    reading0 = False
    output0 = ''
    pass0 = True
    #output1 = []
    for line in fileinput.input():
        if '>' in line:
            if 'OS=Homo sapiens' in line:
                reading0 = True
                if len(output0)>1 and pass0:
                    print output0[:-1]
                #print line[:-1]
                output0 = line[0:4]+'\n'
                pass0 = True
            else:
                reading0 = False
        elif reading0:
            if pass0:
                for item in 'BJOUXZ':
                    if item in line:
                        pass0 = False
                        break
            if pass0:
                output0 = output0+line 
    if len(output0)>1 and pass0:
        print output0[:-1]
         
        

if __name__ == '__main__':
    pass

main()

