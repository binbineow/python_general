import fileinput
import Levenshtein

#path /scratch/users/bchen45/HLA_prediction/RNN_data/training_files/
# HLADRB10101fix_val_withIEDB_1to1_tr_1_val.csv
# HLADRB10101fix_val_withIEDB_tr_1_val.csv
# HLADRB10101simplev1_tr_1_val.csv
# HLADRB10101simplev2_fix_HLA_tr_1_val.csv
# HLADRB10101val_check_fix_HLA_decluster_tr_1_val.csv
# HLADRB10101val_check_fix_HLA_tr_1_val.csv

#initialize the lists
#n_train is actually not used
p_train = []
p_val = []
n_train = []
n_val = []



#input format
#0 - negative training 1 - positive training 2 - negative validation 3 - positive validation        
for line0 in fileinput.input():
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
for _ in range(0,n_iterations):
    i += 0.9/n_iterations
    print('cut_off='+str(i))
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
    if tp+fn>0 and tp+fp>0:
        recall = float(tp)/(tp+fn)
        precision = float(tp)/(tp+fp)

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

    print('precision='+str(precision)+'\trecall='+str(recall)+'\tF1 score='+str(accuracy))
print('best F1 score='+str(best_accuracy)+' achieved at cut_off='+str(best_i))  
#print('best test recall at best i ='+str(test_at_best_i))      
