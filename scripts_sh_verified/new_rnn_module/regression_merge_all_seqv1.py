   
from runn_support.py import *

def regression_merge_all_seq():
    rnn_master = RNN_master()    
    [X_train_fixed,train_seq,y_train] = pickle.load(open(path_data+train_file0,'r'))
    [X_val_fixed,val_seq,y_val] = pickle.load(open(path_data+val_file0,'r'))
    ########encoding
    len_feature = len(X_train_fixed[0])
    model = make_model(len_feature)
    X_train_variable = encoding_data(train_seq,MAXLEN)
    X_val_variable = encoding_data(val_seq,MAXLEN)
    print(X_train_fixed[0])
    print(X_train_variable[0])
    r_best = 0
    output_perf2(['Iteration','Training PCC','Training p-val','Val PCC','Val p-val'])
    print('start training')
    for n0 in range(0,n_iteration+1):
        #fit    
        #print(y_train)
        model.fit([X_train_fixed,X_train_variable], y_train, batch_size=BATCH_SIZE, verbose=vb0, nb_epoch=nb0)      
        #calculate the performance
        #calculate Pearson Correltion Coeficient 
        y_train_pred = model.predict([X_train_fixed,X_train_variable],batch_size=BATCH_SIZE)
        y_train_pred = y_train_pred.reshape(y_train_pred.shape[0])
        [r0_train,pval0_train] = pearsonr(y_train_pred,y_train)
        
        y_predicted = model.predict([X_val_fixed,X_val_variable],batch_size=BATCH_SIZE)
        y_predicted = y_predicted.reshape(y_predicted.shape[0])
        print(y_predicted)
        [r0, pval0] = pearsonr(y_predicted,y_val)
        #print('PCC in validation'+str(r0))
        #save performance
        output_perf2([n0,r0_train,pval0_train,r0,pval0])
        #print performance
        #print([n0,r0_train,pval0_train,r0,pval0])
        #print('Predicted binding aff')
        #print(y_predicted[0:10])
        #print('Measured binding aff')
        #print(y_val[0:10])
        #save the model
        if r0 > r_best+0.005:
            model.save_weights(path_save+file_name0+out_name+'_weight.h5',overwrite=True)
            r_best = r0

#print out parameteres
output_perf2([input_info]) 
regression_merge_all_seq()
