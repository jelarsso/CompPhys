import numpy as np
import plotly.graph_objects as go
import subprocess as sb
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


def compare_p5c():
    sb.run(["./p5c",str(0.01),str(10000)])
    data_f = read_file("forward_euler.data")
    data_b = read_file("backward_euler.data")
    data_cn = read_file("cnicholson.data")
    x = np.linspace(0,1,data_f.shape[1])
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x,y=data_f[-1,:],mode="lines",name="Forward"))
    fig.add_trace(go.Scatter(x=x,y=data_b[-1,:],mode="lines",name="Backward"))
    fig.add_trace(go.Scatter(x=x,y=data_cn[-1,:],mode="lines",name="Cranky"))
    fig.show()




if __name__=="__main__":
    #animate("dump.data")
    compare_p5c()
