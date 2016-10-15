def plot_2_lines(y_pred,y_reco,name0,path0='/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/variable_region_plots/'):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')
    #example name0 = MCL001_H_chain

    len0 = len(y_pred) #len y1 = y2
    x = np.linspace(1, len0, len0)
    #print(x)
    y_pred = np.array(y_pred)
    #print(y_pred)
    y_reco = np.array(y_reco)
    #print(y_reco)
    ax = plt.gca()
    plt.fill(x,y_pred,'b*',x,y_reco,'r*',alpha=0.3)
    ax.set_xlabel(name0+' Amino Acid Sequence')
    ax.set_ylabel('Peptides recovered/Predicted scores')
    ax.set_title(name0+' Peptide Recovered vs. MARIA Prediction Scores')
    plt.savefig(path0+name0+'_pred_vs_recovered.png', bbox_inches='tight')
    plt.close()