from rnn_master import *

#testing
path_para = '/home/stanford/rbaltman/users/bchen45/code/slurm/'
para_file = 'run_now_merge_class.txt'
name0 = 'test'
rnn_master = RNN_master()
rnn_master.get_input(path_para+para_file)
rnn_master.create_out_file_names()
#rnn_master.mask0

#test run
#get data
[mhc_train,seq_train,label_train] = pickle.load(open(rnn_master.path_data+rnn_master.train_file0,'r'))
[mhc_val,seq_val,label_val] = pickle.load(open(rnn_master.path_data+rnn_master.val_file0,'r'))
#get max length
rnn_master.MAXLEN = 35
#max([len(seq0) for seq0 in seq_train]+[len(seq0) for seq0 in seq_val])
print(rnn_master.MAXLEN)
#encode x 
len_hla = len(mhc_train[0])
x_mhc_train = encoding_fixed(mhc_train,len_hla,rnn_master.dict_aa,len(rnn_master.chars),add_placeh=3)
x_seq_train = encoding_data(seq_train,rnn_master.MAXLEN,rnn_master.dict_aa,len(rnn_master.chars))
x_mhc_val = encoding_fixed(mhc_val,len_hla,rnn_master.dict_aa,len(rnn_master.chars),add_placeh=3)
x_seq_val = encoding_data(seq_val,rnn_master.MAXLEN,rnn_master.dict_aa,len(rnn_master.chars))
#process y
y_train = encoding_y(label_train)
y_val = encoding_y(label_val)
label_train = np.array([int(y0) for y0 in label_train])
label_val = np.array([int(y0) for y0 in label_val])
#generate model
fix_len = len(x_mhc_train[0])
rnn_master.fix_len = fix_len

#create model and test
#rnn_master.drop_out_c = 0.4
#rnn_master.BATCH_SIZE = 512
model_merge = rnn_master.make_model('class2')
print(model_merge.summary())
auc_val = cal_performance(model_merge,[x_mhc_val,x_seq_val],label_val)

#load weight from 
path_save = rnn_master.path_save
file0 = rnn_master.file_name0+rnn_master.out_name+'lstm2_weight.h5'
model_merge.load_weights(path_save+file0)

#train model just to see
n_iteration = 1
rnn_master.nb0 = 1
for _ in range(0,n_iteration):
    model_merge.fit([x_mhc_train,x_seq_train],y_train,
                    batch_size=rnn_master.BATCH_SIZE,
                    verbose=1,epochs=rnn_master.nb0,
                    validation_data=([x_mhc_val,x_seq_val],y_val))
#train model
n_iteration = 10
rnn_master.nb0 = 5
vb0 = 0
auc_best = 0.83
for i in range(0,n_iteration):
    print(i)
    model_merge.fit([x_mhc_train,x_seq_train],y_train,
                    batch_size=rnn_master.BATCH_SIZE,
                    verbose=vb0,epochs=rnn_master.nb0,
                    validation_data=([x_mhc_val,x_seq_val],y_val))
    auc_val = cal_performance(model_merge,[x_mhc_val,x_seq_val],label_val)
    #save if better
    if auc_val > auc_best+0.01:
        auc_best = auc_val
        model_merge.save_weights(rnn_master.path_save+ 
                         rnn_master.file_name0+rnn_master.out_name+name0+'_weight.h5',
                         overwrite=True)
        
