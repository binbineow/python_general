 #importing
import fileinput
from utilities import *
path0 = '../images/v2'
label_value = -20
cap_value = 50

#dict_g_u = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_gene_to_uni.dict'))
#dict_u_g = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/dict_uni_to_gene.dict'))
dict_jeko = pickle.load(open('/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/jeko_pep_countv3.dict'))
dict_L128 = pickle.load(open('/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/L128_pep_countv3.dict'))
#dict_L128 = pickle.load(open('/scratch/users/bchen45/HLA_prediction/Uniprot/L128_pep_countv2.dict'))
MCL_data = pickle.load(open('/scratch/users/bchen45/HLA_prediction/MCL_MHC_project/gene_analysis/MCL_data11_18_2015v1.2_UP.dict','r'))


def plot_counter2(gene_list0,title0,x0,y0,path0):
    
    plt.hist(gene_list0)
    plt.title(title0)
    plt.xlabel(x0)
    plt.ylabel(y0)
    save(path0+title0,'png')
    
#pid0 = 128 or Jeko; mhc0 = MHC1 or MHC2
def plot_hist_pep(pid0,mhc0):
    #load data
    MCL_data_sub = list(MCL_data['MCL'+pid0][mhc0+'_gene'])
    if pid0 == 'Jeko':
        dict_sub = dict_jeko
    if pid0 == '128':
        dict_sub = dict_L128
    #main
    filename = 'Histagram_of_MCL'+pid0+'_'+mhc0+'peptide_vs_proteome'
    hist_num = []
    for gene0 in MCL_data_sub:
        if gene0 in dict_sub:
            num0 = dict_sub[gene0]
            if num0 < cap_value:
                hist_num.append(num0)
        #else:
            #hist_num.append(label_value)
               
    plot_counter2(hist_num,filename,'peptide emPAI in Jeko whole cells',mhc0+' Gene frequency',path0)

def plot_hist_pep_norm(pid0):
    #load data
    hist_num = []
    if pid0 == 'Jeko':
        dict_sub = dict_jeko
    if pid0 == '128':
        dict_sub = dict_L128
    for gene0,value0 in dict_sub.iteritems():
        if value0<cap_value:
            hist_num.append(value0)
    filename = 'Histogram_of_MCL'+pid0+'_proteome_distribution'
    plot_counter2(hist_num,filename,'peptide emPAI in Jeko whole cells','Gene frequency',path0)

#main
plot_hist_pep('Jeko','MHC1')
plot_hist_pep('Jeko','MHC2')
plot_hist_pep('128','MHC1')    
plot_hist_pep('128','MHC2')
plot_hist_pep_norm('Jeko')
plot_hist_pep_norm('128')