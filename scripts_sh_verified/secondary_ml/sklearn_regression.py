#import
import sklearn
import cPickle as pickle
import numpy as np

###scaling
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
data_x_new_scaled = scaler.fit_transform(data_x_new)

#linear Regression
from sklearn.linear_model import LinearRegression, BayesianRidge
lr_model = LinearRegression()


def cv_score(data_x,data_y,model0,scoring0 = 'mean_squared_error',cv0=10):
    from sklearn.model_selection import cross_val_score
    import numpy as np
    scores0 = cross_val_score(model0,data_x,data_y,cv=cv0,scoring=scoring0)
    print('average '+scoring0+' = '+str(np.average(scores0)))
    return scores0

scores_lr = cv_score(data_x_new2,data_y,lr_model)
scores_br = cv_score(data_x_new2, data_y, BayesianRidge())


from sklearn.naive_bayes import GaussianNB

def cv_score(data_x,data_y,model0,scoring0 = 'f1',cv0=10):
    from sklearn.model_selection import cross_val_score
    import numpy as np
    scores0 = cross_val_score(model0,data_x,data_y,cv=cv0,scoring=scoring0)
    print('average '+scoring0+' = '+str(np.average(scores0)))
    return scores0

from sklearn.model_selection import train_test_split
x_train, x_val, y_train, y_val = train_test_split(data_x_new2,data_y)

def make_bin(array0,cut_off=0):
    list_out = np.zeros(len(array0))
    mask1 = array0 > cut_off
    list_out[mask1] = 1
    return list_out