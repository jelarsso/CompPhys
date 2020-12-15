import numpy as np
import plotly.graph_objects as go
import plotly.subplots as ps
import subprocess as sb
import pandas as pd
from IPython import embed


def read_file(filename):
    return np.loadtxt(filename)
    



def animate(filename):
    sb.run(["./p5c",str(0.01),str(10000)])
    datas = read_file(filename)
    x = np.linspace(0,1,101)
    y = datas[0,:]
    
    fig = go.Figure(
    data=[go.Scatter(x=x, y=y+x)],
    layout=go.Layout(
        xaxis=dict(range=[np.min(x), np.max(x)], autorange=False),
        yaxis=dict(range=[np.min(datas+x), np.max(datas+x)], autorange=False),
        title="Start Title",
        updatemenus=[dict(
            type="buttons",
            buttons=[dict(label="Play",
                          method="animate",
                          args=[None,{"frame": {"duration": 1, "redraw": True},"fromcurrent": False, "transition": {"duration": 0}}])])]
    ),
    frames=[go.Frame(data=[go.Scatter(x=x, y=datas[i,:]+x)],layout=go.Layout(title_text=f"Frame {i}/{datas[:,0].size}")) for i in range(0,datas.shape[0],100)])
    fig.show()


def compare_p5c():
    for dx in [0.1,0.01]:
        T = 0.1
        dt = 0.5*dx**2
        nt = int(round(T/dt))
        sb.run(["./p5c",str(dx),str(nt)])
        data_f = read_file("forward_euler.data")
        data_b = read_file("backward_euler.data")
        data_cn = read_file("cnicholson.data")
        data_an = read_file("analytical.data")
        x = np.linspace(0,1,data_f.shape[1])
        nx = data_f.shape[1]
        print(nx)
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x,y=x+data_f[-1,:],mode="lines",name="Forward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=x+data_b[-1,:],mode="lines",name="Backward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=x+data_cn[-1,:],mode="lines",name="Crank-Nicholson"))
        fig.add_trace(go.Scatter(x=x,y=x+data_an[-1,:],mode="lines",name="Analytical"))
        
        
        fig.update_xaxes(title="x")
        fig.update_yaxes(title="u(T,x)")
        fig.update_layout(font_family="lmodern",title_text=f"Solutions after T={nt} timesteps, dx = {dx}, alpha = 0.5",font_size=12)
        fig.write_image(f"p5c_comparisons_long_dx{dx}.pdf",width=600*1.41,height=600,scale=2)
        fig.show()
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x,y=x+data_f[nt//10,:],mode="lines",name="Forward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=x+data_b[nt//10,:],mode="lines",name="Backward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=x+data_cn[nt//10,:],mode="lines",name="Crank-Nicholson"))
        fig.add_trace(go.Scatter(x=x,y=x+data_an[nt//10,:],mode="lines",name="Analytical"))
        fig.update_xaxes(title="x")
        fig.update_yaxes(title="u(T,x)")
        fig.update_layout(font_family="lmodern",title_text=f"Solutions after T={nt//8} timesteps, dx = {dx}, alpha = 0.5",font_size=12)
        fig.write_image(f"p5c_comparisons_short_dx{dx}.pdf",width=600*1.41,height=600,scale=2)
        fig.show()
        
    
def animate_2d():
    data = np.loadtxt("fe2d.data")
    nt = 1001
    nx = 101
    print(data.shape)
    data = data.reshape((nt,nx,nx))
    

    x = np.linspace(0,1,nx)
    y = np.linspace(0,1,nx)
    z = data[0,:,:]
    
    fig = go.Figure(
    data=[go.Surface(z=z.T,x=x, y=y)],
    layout=go.Layout(
        xaxis=dict(range=[np.min(x), np.max(x)], autorange=False),
        yaxis=dict(range=[np.min(y), np.max(y)], autorange=False),
        title="Start Title",
        updatemenus=[dict(
            type="buttons",
            buttons=[dict(label="Play",
                          method="animate",
                          args=[None,{"frame": {"duration": 1, "redraw": True},"fromcurrent": False, "transition": {"duration": 0}}])])]
    ),
    frames=[go.Frame(data=[go.Surface(z=data[i,:,:].T,x=x, y=y)],layout=go.Layout(title_text=f"Frame {i+1}/{data[:,0,0].size}")) for i in range(0,data.shape[0],10)])
    fig.show()


def show_differences_litho():
    nt = 1601
    nx = 160
    islice = 80
    data = np.loadtxt("litho_enriched.data")
    data_enriched = data.reshape((nt,nx,nx))
    data = np.loadtxt("litho_no_Q.data")
    data_steady = data.reshape((nt,nx,nx))
    data = np.loadtxt("litho_Q_pb.data")
    data_steay2 = data.reshape((nt,nx,nx))

    x = np.linspace(0,1,nx)

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x,y=(data_steay2[-1,islice,:]-data_steady[-1,islice,:])/data_steady[-1,islice,:],mode="lines",name="Relative unenriched"))
    fig.add_trace(go.Scatter(x=x,y=(data_enriched[-1,islice,:]-data_steady[-1,islice,:])/data_steady[-1,islice,:],mode="lines",name="Relative enriched"))
    fig.update_xaxes(title="Depth / L")
    fig.update_yaxes(title="Relative Difference / ")
    fig.update_layout(font_family="lmodern",title_text=f"Difference from the steady state after {nt} timesteps",font_size=12)
    #fig.write_image(f"p5c_comparisons_long_dx{dx}.pdf",width=600*1.41,height=600,scale=2)
    fig.show()


    ds = data_steady[-1,islice,:]
    ds2 = data_steay2[-1,islice,:]
    de = data_enriched[:,islice,:]

    fig = go.Figure(
    data=[go.Scatter(x=x, y=(de[0,:]-ds)/ds)],
    layout=go.Layout(
        xaxis=dict(range=[np.min(x), np.max(x)], autorange=False),
        yaxis=dict(range=[np.min((de[:,:]-ds)/ds), np.max((de[:,:]-ds)/ds)], autorange=False),
        title="Start Title",
        updatemenus=[dict(
            type="buttons",
            buttons=[dict(label="Play",
                          method="animate",
                          args=[None,{"frame": {"duration": 1, "redraw": True},"fromcurrent": False, "transition": {"duration": 0}}])])]
    ),
    frames=[go.Frame(data=[go.Scatter(x=x, y=(de[i,:]-ds)/ds)],layout=go.Layout(title_text=f"Frame {i}/{nt}")) for i in range(0,nt,10)])
    fig.show()
    embed()





    

if __name__=="__main__":
    #animate("cnicholson.data")
    compare_p5c()
    #animate_2d()
    #show_differences_litho()
