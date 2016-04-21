import fileinput
import Levenshtein
from utilities import *
#from matplotlib.dviread import fPrev
path_train = '/scratch/users/bchen45/HLA_prediction/RNN_data/training_files/'
path_AUC = '/scratch/users/bchen45/code/python_general/python_general/images/'

#path /scratch/users/bchen45/HLA_prediction/RNN_data/training_files/
# HLADRB10101fix_val_withIEDB_1to1_tr_1_val.csv
# HLADRB10101fix_val_withIEDB_tr_1_val.csv
# HLADRB10101simplev1_tr_1_val.csv
# HLADRB10101simplev2_fix_HLA_tr_1_val.csv
# HLADRB10101val_check_fix_HLA_decluster_tr_1_val.csv
# HLADRB10101val_check_fix_HLA_tr_1_val.csv
#This version will calculate AUC

#initialize the lists
#n_train is actually not used


def plot_scatter(x,y,x_name,y_name,title,filename1,path0=''):
    plt.figure()
    plt.plot(x,y,'b.')
    plt.ylabel(y_name)
    plt.xlabel(x_name)
    plt.ylim([0,1])
    plt.xlim([0,1])
    plt.title(title)
    save(path0+filename1)

for filename0 in open(path_train+'file_names_IEDB_MCL_1to1.csv','r'):
    filename0 = filename0.rstrip()
    p_train = []
    p_val = []
    n_train = []
    n_val = []
    #input format
    #0 - negative training 1 - positive training 2 - negative validation 3 - positive validation        
    for line0 in open(path_train+filename0,'r'):
        line0 = line0.rstrip()
        str0 = line0.split('\t')[0]
        if line0.split('\t')[1] == '1':
            p_train.append(str0)
        elif line0.split('\t')[1] == '0': 
            n_train.append(str0)
        elif line0.split('\t')[1] == '2':
            n_val.append(str0)
        elif line0.split('\t')[1] == '3':
            p_val.append(str0)
    
    print('Positive training peptide number: '+str(len(p_train)))
    p_score = []
    
    
    for x in p_val:
        score0 = 0
        for y in p_train:
            if Levenshtein.ratio(x,y)>score0:
                score0 = Levenshtein.ratio(x,y)
        p_score.append(score0)
    n_score = []
    for x in n_val:
        score0 = 0
        for y in p_train:
            if Levenshtein.ratio(x,y)>score0:
                score0 = Levenshtein.ratio(x,y)
        n_score.append(score0)
    # t_score = []
    # for x in p_train:
    #     score0 = 0
    #     for y in p_train:
    #         if Levenshtein.ratio(x,y)>score0:
    #             score0 = Levenshtein.ratio(x,y)
    #     t_score.append(score0)
        
    #print(t_score)
    best_i = 0
    best_accuracy = 0
    test_at_best_i = 0
    i = 0.1
    n_iterations = 100
    tpr = []
    fpr = []
    report0 = False
    for _ in range(0,n_iterations):
        i += 0.9/n_iterations
        #print('cut_off='+str(i))
        tp = 0
        tn = 0
        fp = 0
        fn = 0
        for x in p_score:
            if x>i:
                tp +=1
            else:
                fn +=1
        for x in n_score:
            if x>i:
                fp +=1
            else:
                tn +=1
        
        if tp+fn>0 and tp+fp>0 and fp+tn>0:
            recall = float(tp)/(tp+fn)
            precision = float(tp)/(tp+fp)
            tpr.append(float(tp)/(tp+fn))
            fpr.append(float(fp)/(fp+tn))
    
        else:
            recall = 0
            precision = 0
        #accuracy = F1 score
        if precision > 0 and recall >0:
            accuracy = 2*recall*precision/(recall+precision)
        else:
            accuracy = 0
    #     tp_t = 0
    #     for x in t_score:
    #         if x>i:
    #             tp_t +=1
    #     recall_t = float(tp_t)/len(t_score)
        if accuracy>best_accuracy:
            best_accuracy = accuracy
            best_i = i
            #test_at_best_i = recall_t
        if abs(i-0.65)<0.1 and not report0:
            print('F1 socre='+str(accuarcy)+' at cut-off=0.65')
    
        #print('precision='+str(precision)+'\trecall='+str(recall)+'\tF1 score='+str(accuracy))
    print('best F1 score='+str(best_accuracy)+' achieved at cut_off='+str(best_i))
    name1 = filename0.split('_')[0]
    #AUC
    from sklearn import metrics
    #print fpr
    #print tpr
    auc0 = metrics.auc(fpr,tpr,reorder=True)
    print('AUC= '+str(auc0)) 
    plot_scatter(fpr, tpr, '1-specificity', 'Sensitivity', name1+' AUC='+str(auc0), name1+'_AUC',path_AUC)

    #print('best test recall at best i ='+str(test_at_best_i))      
