import cPickle as pickle
import plotly.plotly as py
import plotly.graph_objs as go
import numpy as np
path0 = '/Users/binbineow2/Documents/Machine_Learning/HLA_prediction/Figure_MCL/heatmap/'
pid_list = pickle.load(open(path0+'mcl_patient2.list','r'))
z = pickle.load(open(path0+'normal_j_matrixv2.np','r'))
#z = pickle.load(open(path0+'heatmap_drb.mat','r'))
#for n0 in range(0,len(z)):
#    #z[n0,n0] = 0
#    z[n0] = list(reversed(z[n0]))
z = np.flipud(z)
x = pid_list
#x = list(reversed(pid_list))
y = list(reversed(pid_list))

annotations = []
for n, row in enumerate(z):
    for m, val in enumerate(row):
        var = z[m][n]
        annotations.append(
            dict(
                 
                text=str(val),
                x=x[m], y=y[n],
                xref='x1', yref='y1',
                font=dict(color='white' if val > 0.5 else 'black'),
                showarrow=False)
            )

colorscale = [[0, '#FFFFFF'], [0.5, '#003366'],[1, '#003366']]  # custom colorscale
#trace = go.Heatmap(x=x, y=y, z=z, colorscale=colorscale, showscale=False)
trace = go.Heatmap(x=x, y=y, z=z, colorscale=colorscale)
fig = go.Figure(data=[trace])
fig['layout'].update(
    title="Overlapping of peptides between patients (Standard Jaccard Index)",
    #title="HLA-DRB1 Sequence similarity between patients (Levenshtein Ratio)",
    #annotations=annotations,
    xaxis=dict(ticks='', side='top'),
    # ticksuffix is a workaround to add a bit of padding
    yaxis=dict(ticks='', ticksuffix='  '),
    width=900,
    height=900,
    autosize=True
)

url = py.plot(fig, filename='Boostrapped Peptide Overlap Heatmap', height=1000)
