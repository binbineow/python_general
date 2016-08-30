def cal_auc(p_list,n_list,n_iterations):
    tpr = []
    fpr = []
    from sklearn import metrics
    max_cutoff = max (p_list+n_list)
    i = min(p_list+n_list)
    for _ in range(0,n_iterations):
        i += float(max_cutoff)/n_iterations
        #print('cut_off='+str(i))
        tp = 0.1
        tn = 0.1
        fp = 0.1
        fn = 0.1
    
        for x in p_list:
            if x>i:
                tp +=1
            else:
                fn +=1
        for x in n_list:
            if x>i:
                fp +=1
            else:
                tn +=1
            
        if tp+fn>0 and tp+fp>0 and fp+tn>0:
            #recall = float(tp)/(tp+fn)
            #precision = float(tp)/(tp+fp)
            #print(str(tpr)+''+str(fpr))
            tpr.append(float(tp)/(tp+fn))
            fpr.append(float(fp)/(fp+tn))
    return metrics.auc(fpr,tpr,reorder=True)



