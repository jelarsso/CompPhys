import numpy as np
import plotly.graph_objects as go
import pandas as pd
from IPython import embed


def read_file(filename):
    return np.loadtxt(filename)
    



def animate(filename):
    datas = read_file(filename)
    x = np.linspace(0,1,101)
    y = datas[0,:]
    
    fig = go.Figure(
    data=[go.Scatter(x=x, y=y)],
    layout=go.Layout(
        xaxis=dict(range=[np.min(x), np.max(x)], autorange=False),
        yaxis=dict(range=[np.min(datas), np.max(datas)], autorange=False),
        title="Start Title",
        updatemenus=[dict(
            type="buttons",
            buttons=[dict(label="Play",
                          method="animate",
                          args=[None,{"frame": {"duration": 1, "redraw": True},"fromcurrent": False, "transition": {"duration": 0}}])])]
    ),
    frames=[go.Frame(data=[go.Scatter(x=x, y=datas[i,:])],layout=go.Layout(title_text=f"Frame {i}/{datas[:,0].size}")) for i in range(0,datas.shape[0],100)])
    fig.show()

animate("dump.data")

