import cPickle as pickle
import plotly.plotly as py
import plotly.graph_objs as go
from scipy.stats.stats import pearsonr
path0 = '/Users/binbineow2/Documents/Machine_Learning/HLA_prediction/Figure_MCL/heatmap/'
pid_list = pickle.load(open(path0+'heatmap_pid.list','r'))
z_pep = pickle.load(open(path0+'heatmap_mhc2_J.mat','r'))
z_DRB1 = pickle.load(open(path0+'heatmap_drb.mat','r'))
pep_plot = []
DRB1_plot = []
for n0 in range(0,len(z_pep)):
    for m0 in range(n0,len(z_pep[n0])):
        if not (z_pep[n0][m0] == 1 and z_DRB1[n0][m0] == 1):
            pep_plot.append(z_pep[n0][m0])
            DRB1_plot.append(z_DRB1[n0][m0])

[r,p] = pearsonr(pep_plot, DRB1_plot)
print ([r,p])
#######plotting
import plotly.plotly as py
import plotly.graph_objs as go


# Create a trace
trace = go.Scatter(
    x = DRB1_plot,
    y = pep_plot,
    mode = 'markers'
)
data=[trace]
fig = go.Figure(data=[trace])
fig['layout'].update(
    title="Regression between Class II peptide similarity and HLA-DRB1 similarity\nR=0.50  p=9.4e-11",
    #title="HLA-DRB1 Sequence similarity between patients (Levenshtein Ratio)",
    #annotations=annotations,
    xaxis=dict(title = 'HLA-DRB1 similarity score', ticks='-', side='bottom'),
    # ticksuffix is a workaround to add a bit of padding
    yaxis=dict(title = 'MHCII Peptide similarity score',ticks='-'),
    width=900,
    height=900,
    autosize=True
)
url = py.plot(fig, filename='Regressionv2',height=500)

# Plot and embed in ipython notebook!
#py.iplot(data, filename='basic-scatter')

#plot_url = py.plot(data, filename='basic-line')
    


