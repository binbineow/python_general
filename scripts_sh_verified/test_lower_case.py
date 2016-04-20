import fileinput
path_save_train = '/scratch/users/bchen45/HLA_prediction/RNN_data/training_files/'

def test_lower(str0):
    for a in str0:
        if not a in 'ARNDCQEGHILKMFPSTWYVBZX':
            print str0
            break

for file_name0 in open(path_save_train+'file_names_IEDBv1.csv'):
    #model = model1
    file_name0 = file_name0.rstrip()
    for line0 in open(path_save_train+file_name0):
        line0 = line0.split('\t')[0]
        test_lower(line0)
        

