from sklearn import metrics
import fileinput


def(p_list,n_list,n_iteration):

value_dict = dict()
p_list = []
n_list = []

for line0 in fileinput.input():
    if not 'allele' in line0 and len(line0) > 3 :
        #print line0
        line0 = line0.split(',')
        num0 = int(line0[1])
        val0 = float(line0[-2])
        if not num0 in value_dict or value_dict[num0] > val0:
            value_dict[num0] = val0

for key0, value0 in value_dict.iteritems():
    if key0 % 2 == 0:
        n_list.append(value0)
    else:
        p_list.append(value0)

n_iterations = 1000
max_cutoff = 10000
best_i = 0
best_accuracy = 0
test_at_best_i = 0
i = 0
tpr = []
fpr = []
report0 = False    

for _ in range(0,n_iterations):
    i += max_cutoff/n_iterations
    #print('cut_off='+str(i))
    tp = 0
    tn = 0
    fp = 0
    fn = 0

    for x in p_list:
        if x<i:
            tp +=1
        else:
            fn +=1
    for x in n_list:
        if x<i:
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
    if abs(i-500)<10 and not report0:
        print('F1 socre='+str(accuracy)+' at cut-off='+str(i))
        report0 = True
    
        #print('precision='+str(precision)+'\trecall='+str(recall)+'\tF1 score='+str(accuracy))
print('best F1 score='+str(best_accuracy)+' achieved at cut_off='+str(best_i))
    #name1 = filename0.split('_')[0]
    #AUC
from sklearn import metrics
    #print fpr
    #print tpr
auc0 = metrics.auc(fpr,tpr,reorder=True)
print('AUC= '+str(auc0)) 
    
    
    
    
    

