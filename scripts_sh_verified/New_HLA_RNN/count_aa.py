#this function reads in a list of strings and count the frequency of 
#each letters

def count_aa(list0):
    str0 = ''.join(list0)
    alpha_list = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    print len(alpha_list)
    for let0 in alpha_list:
        print(let0+': '+str(str0.count(let0)))