def plot_2_lines(y_pred,y_reco,name0,path0='/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/variable_region_plots/'):
    import numpy as np
    import matplotlib.pyplot as plt
    #example name0 = MCL001_H_chain
    len0 = lne(y_pred) #len y1 = y2
    x = np.linspace(1, len0+1, len0+1)
    plt.fill(x,y_pred,'b',x,y_reco,'r',alpha=0.3)
    plt.savefig(path0+name0+'_pred_vs_recovered.png', bbox_inches='tight')
