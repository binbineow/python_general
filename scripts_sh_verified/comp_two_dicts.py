 #importing
from utilities import *
from scipy.stats import pearsonr
path0 = '../images/'
def plot_scatter(x,y,x_name,y_name,title,path0=''):
    plt.figure()
    plt.plot(x,y,'bo')
    plt.ylabel(y_name)
    plt.xlabel(x_name)
    plt.ylim([0,max(y)+1])
    plt.title(title)
    save(path0+title)
    
# #reference for plot function
# def plot_counter2(gene_list0,filename0,title0,y0,x0):
#     
#     plt.hist(gene_list0)
#     plt.title(title0)
#     plt.xlabel(x0)
#     plt.ylabel(y0)
#     save(filename0,'png')

#get input
for line0 in fileinput.input():
    line0 = line0.rstrip().split(';')
    if line0[0] == 'dict1':
        dict1 = pickle.load(open(line0[1]))
    if line0[0] == 'dict1_name':
        dict1_name = line0[1]
    if line0[0] == 'dict2_name':
        dict2_name = line0[1]
    if line0[0] == 'title':
        title = line0[1]
    if line0[0] == 'dict2':
        dict2 = pickle.load(open(line0[1]))

#plot correlation scatter
#x is dict1
x0 = []
#y is dict2
y0 = []
for name0,value0 in dict1.iteritems():
    if name0 in dict2:
        x0.append(float(value0))
        y0.append(float(dict2[name0]))
print(pearsonr(x0,y0))

plot_scatter(x0,y0,dict1_name,dict2_name,title,path0)        



        
