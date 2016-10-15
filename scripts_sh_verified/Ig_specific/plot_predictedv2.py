def plot_2_lines(y_pred,y_reco,name0,chain0,path0='/home/stanford/rbaltman/users/bchen45/data/MCL_data/ig_specific/variable_region_plots/'):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    import subprocess
    matplotlib.use('Agg')
    #example name0 = MCL001_H_chain
    cmd_line = 'mkdir '+path0+chain0
    cmd0 = subprocess.Popen(cmd_line,shell=True)      
    cmd0.wait()
    len0 = len(y_pred) #len y1 = y2
    x = np.linspace(0, len0+1, len0+2)
    #print(x)
    y_pred = np.array([0]+y_pred+[0])
    #print(y_pred)
    y_reco = np.array([0]+y_reco+[0])
    #print(y_reco)
    ax = plt.gca()
    alpha0 = 0.3
    pred0, = plt.fill(x,y_pred,'b',alpha=alpha0,label='MARIA Predicted')
    reco0, = plt.fill(x,y_pred,'r',alpha=alpha0,label='MHCII Peptide Recovered')
    plt.legend(handles=[pred0,reco0,])
    #plt.fill(x,y_pred,'b*',alphax,y_reco,'r*',alpha=0.2)
    ax.set_xlabel(name0+' Amino Acid Sequence')
    ax.set_ylabel('Peptides recovered/Predicted scores')
    ax.set_title(name0+' Peptide Recovered vs. MARIA Prediction Scores')
    plt.savefig(path0+chain0+'/'+name0+'_pred_vs_recovered.png', bbox_inches='tight')
    plt.close()