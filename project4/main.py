import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np
import subprocess as sb
from IPython import embed

def read_datafile(filename):
    with open(filename,"r") as infile:
        l = infile.readline().split(" ")
        n_spins = int(l[1]) 
        n_mc = int(l[3])
    return n_spins,n_mc,np.loadtxt(filename)



def p4a_analytical():
    filename = "p4c.data"
    n_mc = 10000
    start_temp = 0.5
    stop_temp = 5
    step_temp = 0.1

    beta = 1/np.arange(start_temp,stop_temp,step_temp)
    Z = (12 + 4*np.cosh(8*beta))
    mean_E = -32*np.sinh(8*beta)/Z
    mean_E2 = 4*64*np.cosh(8*beta)/Z
    mean_M = (16+8*np.exp(8*beta))/Z
    mean_M2 = 32*(1+np.exp(8*beta))/Z

    fig = make_subplots(rows=2, cols=2, subplot_titles=("Average Energy","Mean Magnetization", "Specific heat", "Susceptibility"))
    
    fig.add_trace(go.Scatter(x=1/beta,y=mean_E),row=1,col=1)
    fig.add_trace(go.Scatter(x=1/beta,y=mean_M),row=1,col=2)
    fig.add_trace(go.Scatter(x=1/beta,y=beta*beta*(mean_E2-mean_E**2)),row=2,col=1)
    fig.add_trace(go.Scatter(x=1/beta,y=beta*(mean_M2-mean_M**2)),row=2,col=2)
    
    sb.run(["./spins2", filename, str(n_mc), str(start_temp),str(stop_temp),str(step_temp)])
    ns,nmc,data = read_datafile(filename)

    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,1]),row=1,col=1)
    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,5]),row=1,col=2)
    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,2]),row=2,col=1)
    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,4]),row=2,col=2)

    fig.show()

def p4c_comapre_nmc():
    filename = "p4c.data"
    n_mc = np.logspace(1,7,50)
    temp = 1
    d = []

    beta = 1
    Z = (12 + 4*np.cosh(8*beta))
    mean_E = -32*np.sinh(8*beta)/Z
    mean_E2 = 4*64*np.cosh(8*beta)/Z
    mean_M = (16+8*np.exp(8*beta))/Z
    mean_M2 = 32*(1+np.exp(8*beta))/Z

    fig = make_subplots(rows=2, cols=2, subplot_titles=("Average Energy","Mean Magnetization", "Specific heat", "Susceptibility"))
    fig.add_trace(go.Scatter(x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[mean_E,mean_E]),row=1,col=1)
    fig.add_trace(go.Scatter(x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[mean_M,mean_M]),row=1,col=2)
    fig.add_trace(go.Scatter(x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[beta*beta*(mean_E2-mean_E**2),beta*beta*(mean_E2-mean_E**2)]),row=2,col=1)
    fig.add_trace(go.Scatter(x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[beta*(mean_M2-mean_M**2), beta*(mean_M2-mean_M**2)]),row=2,col=2)
    
    for n in n_mc:
        sb.run(["./spins2", filename, str(n), str(temp),str(temp),str(1)])
        ns,nmc,data = read_datafile(filename)
        d.append(data[:])
    d = np.asarray(d)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,1]),row=1,col=1)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,5]),row=1,col=2)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,2]),row=2,col=1)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,4]),row=2,col=2)

    fig.show()







if __name__ == "__main__":
    p4c_comapre_nmc()