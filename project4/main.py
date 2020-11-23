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
    filename = "p4a.data"
    nspins = 2
    n_mc = 100000
    start_temp = 1.0
    stop_temp = 2.6
    step_temp = 0.05

    beta = 1/np.arange(start_temp,stop_temp+step_temp,step_temp)
    Z = (12 + 4*np.cosh(8*beta))
    mean_E = -32*np.sinh(8*beta)/Z
    mean_E2 = 4*64*np.cosh(8*beta)/Z
    mean_M = (16+8*np.exp(8*beta))/Z
    mean_M2 = 32*(1+np.exp(8*beta))/Z

    fig = make_subplots(rows=2, cols=2,  subplot_titles=("Average Energy","Mean Magnetization", "Specific heat", "Susceptibility"),vertical_spacing=0.2,horizontal_spacing=0.10)
     
    fig.add_trace(go.Scatter(name="Analytical",x=1/beta,y=mean_E/nspins/nspins,line=dict(color="Crimson")),row=1,col=1)
    fig.add_trace(go.Scatter(x=1/beta,y=mean_M/nspins/nspins,showlegend=False,line=dict(color="Crimson")),row=1,col=2)
    fig.add_trace(go.Scatter(x=1/beta,y=beta*beta*(mean_E2-mean_E**2)/nspins/nspins,showlegend=False,line=dict(color="Crimson")),row=2,col=1)
    fig.add_trace(go.Scatter(x=1/beta,y=beta*(mean_M2-mean_M**2)/nspins/nspins,showlegend=False,line=dict(color="Crimson")),row=2,col=2)

    #sb.run(["./spins2", filename, str(n_mc), str(start_temp),str(stop_temp),str(step_temp)])
    ns,nmc,data = read_datafile(filename)

    fig.add_trace(go.Scatter(name="Simulated",x=data[:,0],y=data[:,1],mode="markers",marker=dict(color="MediumPurple")),row=1,col=1)
    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,5],mode="markers",showlegend=False,marker=dict(color="MediumPurple")),row=1,col=2)
    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,2],mode="markers",showlegend=False,marker=dict(color="MediumPurple")),row=2,col=1)
    fig.add_trace(go.Scatter(x=data[:,0],y=data[:,4],mode="markers",showlegend=False,marker=dict(color="MediumPurple")),row=2,col=2)

    fig.update_xaxes(title="Temperature / J/k")
    fig.update_yaxes(title="Energy / J",row=1,col=1)
    fig.update_yaxes(title="Magnetization / #",row=1,col=2)
    fig.update_yaxes(title="Specific Heat / J",row=2,col=1)
    fig.update_yaxes(title="Susceptibility / ",row=2,col=2)
    fig.update_yaxes(title_standoff=1)
    fig.update_layout(title_text="2*2 grid. All values per spin.",font_family="lmodern",font_size=12)
    fig.write_image("exval4c.pdf",width=600*1.41,height=600,scale=2)
    fig.show()

def p4c_comapre_nmc():
    filename = "p4c.data"
    n_mc = np.logspace(1,7,50)
    temp = 1
    d = []


    nspins = 2

    beta = 1
    Z = (12 + 4*np.cosh(8*beta))
    mean_E = -32*np.sinh(8*beta)/Z
    mean_E2 = 4*64*np.cosh(8*beta)/Z
    mean_M = (16+8*np.exp(8*beta))/Z
    mean_M2 = 32*(1+np.exp(8*beta))/Z

    fig = make_subplots(rows=2, cols=2, subplot_titles=("Average Energy","Mean Magnetization", "Specific heat", "Susceptibility"))
    fig.add_trace(go.Scatter(mode="lines",name="Analytical",x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[mean_E/nspins/nspins,mean_E/nspins/nspins],line=dict(color="Crimson")),row=1,col=1)
    fig.add_trace(go.Scatter(mode="lines",showlegend=False,x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[mean_M/nspins/nspins,mean_M/nspins/nspins],line=dict(color="Crimson")),row=1,col=2)
    fig.add_trace(go.Scatter(mode="lines",showlegend=False,x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[beta*beta*(mean_E2-mean_E**2)/nspins/nspins,beta*beta/nspins/nspins*(mean_E2-mean_E**2)],line=dict(color="Crimson")),row=2,col=1)
    fig.add_trace(go.Scatter(mode="lines",showlegend=False,x=[np.log10(n_mc[0]),np.log10(n_mc[-1])], y=[beta/nspins/nspins*(mean_M2-mean_M**2), beta/nspins/nspins*(mean_M2-mean_M**2)],line=dict(color="Crimson")),row=2,col=2)
    
    for n in n_mc:
        sb.run(["./spins2", filename, str(n), str(temp),str(temp),str(1)])
        ns,nmc,data = read_datafile(filename)
        d.append(data[:])
    d = np.asarray(d)
    fig.add_trace(go.Scatter(name="Simulated",x=np.log10(n_mc),y=d[:,1],mode="lines+markers",marker=dict(color="MediumPurple")),row=1,col=1)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,5],mode="lines+markers",showlegend=False,marker=dict(color="MediumPurple")),row=1,col=2)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,2],mode="lines+markers",showlegend=False,marker=dict(color="MediumPurple")),row=2,col=1)
    fig.add_trace(go.Scatter(x=np.log10(n_mc),y=d[:,4],mode="lines+markers",showlegend=False,marker=dict(color="MediumPurple")),row=2,col=2)

    fig.update_xaxes(title="log(Number of Monte Carlo cycles)")
    fig.update_yaxes(title="Energy / J",row=1,col=1)
    fig.update_yaxes(title="Magnetization / #",row=1,col=2)
    fig.update_yaxes(title="Specific Heat / J",row=2,col=1)
    fig.update_yaxes(title="Susceptibility / ",row=2,col=2)
    fig.update_yaxes(title_standoff=1)
    fig.update_layout(font_family="lmodern",title_text="2*2 grid, T=1 J/k, All Values per Spin.",font_size=12)
    fig.write_image("nmc4c.pdf",width=600*1.41,height=600,scale=2)
    fig.show()

def p4e_pde():
    nmc = 100_000
    equiltime = 10_000
    a=sb.run(["./pde", "pde.data", str(nmc), str(equiltime), "1", "2.4", "1.4"],capture_output=True)
    print(a)
    with open("pde.data","r") as infile:
        infile.readline()
        infile.readline()
        infile.readline()
        pde_1 = np.asarray([float(x) for x in infile.readline().split()])
        infile.readline()
        infile.readline()
        pde_24 = np.asarray([float(x) for x in infile.readline().split()])

    fig = make_subplots(rows=1, cols=2, subplot_titles=("T=1 J/k","T=2.4 J/k"))
    fig.add_trace(go.Histogram(showlegend=False,x=pde_1, histnorm='probability density'),row=1,col=1)
    fig.add_trace(go.Histogram(showlegend=False,x=pde_24, histnorm='probability density'),row=1,col=2)
    fig.update_yaxes(title="Probability density / P(E)")
    fig.update_xaxes(title="Energy / J")
    fig.update_layout(title_text="Probability Density - P(E)")
    fig.write_image("pde4e.pdf",width=600*1.41,height=600,scale=2)
    fig.show()

def p4f_many_spin():
    n_mc = 10_000_000
    equiltime = 100_000
    Ls = [40,60,80,100]
    start_temp = [2.25,2.25,2.25,2.25]
    stop_temp = [2.35,2.35,2.35,2.35]
    step_temp = 0.005

    for i,k in enumerate(Ls):    
        filename = f"p4f_l{Ls[i]}_dT{step_temp}_narrower10mill.data"
        fig = make_subplots(rows=2, cols=2, subplot_titles=("Average Energy","Mean Magnetization", "Specific heat", "Susceptibility"))
        sb.run(["./para", filename, str(n_mc),str(equiltime), str(Ls[i]), str(start_temp[i]),str(stop_temp[i]),str(step_temp)])
        ns,nmc,data = read_datafile(filename)
        sort = np.argsort(data[:,0]) # due to parallelization
        fig.add_trace(go.Scatter(showlegend=False,mode="markers",x=data[sort,0],y=data[sort,1]),row=1,col=1)
        fig.add_trace(go.Scatter(showlegend=False,mode="markers",x=data[sort,0],y=data[sort,5]),row=1,col=2)
        fig.add_trace(go.Scatter(showlegend=False,mode="markers",x=data[sort,0],y=data[sort,2]),row=2,col=1)
        fig.add_trace(go.Scatter(showlegend=False,mode="markers",x=data[sort,0],y=data[sort,4]),row=2,col=2)
        fig.update_xaxes(title="Temperature / J/k")
        fig.update_yaxes(title="Energy / J",row=1,col=1)
        fig.update_yaxes(title="Magnetization / #",row=1,col=2)
        fig.update_yaxes(title="Specific Heat / J",row=2,col=1)
        fig.update_yaxes(title="Susceptibility / ",row=2,col=2)
        fig.update_yaxes(title_standoff=1)
        fig.update_layout(title_text=f"{Ls[i]}*{Ls[i]}"+" grid. All values per spin.",font_family="lmodern",font_size=12)
        fig.write_image(f"p4f_l{Ls[i]}_dT{step_temp}_narrower10mill.pdf",width=600*1.41,height=600,scale=2)
        fig.show()
        print(f"L = {Ls[i]}")
        print(f"Maximum value for C_V at {data[sort,0][np.argmax(data[sort,2])]} for suscep {data[sort,0][np.argmax(data[sort,4])]}")

def p4g():
    T_cs = []
    Ls = []

    """0.001, nmc = 1e5
    40 Maximum value for C_V at 2.274 for suscep 2.324
    60 Maximum value for C_V at 2.279 for suscep 2.301
    80 Maximum value for C_V at 2.277 for suscep 2.293
	100 Maximum value for C_V at 2.252 for suscep 2.293
    """

    """ 0.001 nmc = 1e6
    L = 40
    Maximum value for C_V at 2.304 for suscep 2.323
    L = 60
    Maximum value for C_V at 2.283 for suscep 2.295
    L = 80
    Maximum value for C_V at 2.281 for suscep 2.300
    L = 100
    Maximum value for C_V at 2.28 for suscep 2.280
    """

    """
    L = 40
    Maximum value for C_V at 2.285 for suscep 2.32
    L = 60
    Maximum value for C_V at 2.285 for suscep 2.305
    L = 80
    Maximum value for C_V at 2.28 for suscep 2.295
    L = 100
    Maximum value for C_V at 2.275 for suscep 2.29
    """

    L = np.array((40,60,80,100))
    Tc = np.array((2.32,2.305,2.295,2.29))

    print(np.polyfit(1/L,Tc,1,full=False,cov=True))
    print(np.sqrt(2.31801817e-06))

def equil():  
    N = np.linspace(1, 1e3, int(1e3))

    f = open("output1.data", "r")
    f.readline(); f.readline()

    E = np.zeros_like(N)
    M = np.zeros_like(N)
    for i in range(int(1e3)):
        sent = f.readline()
        words = sent.split()
        E[i] = words[1]
        M[i] = words[5]

    fig = make_subplots(rows=2, cols=1, subplot_titles=("Average energy per spin","Average absolute magnetization per spin"), shared_xaxes=True, vertical_spacing=0.1)

    fig.add_trace(go.Scatter(x=N, y=E, name='Average energy per spin'),row=1,col=1)
    fig.add_trace(go.Scatter(x=N,y=M, name='Average absolute magnetization per spin'),row=2,col=1)

    # prøver noe

    f = open("output2.data", "r")
    f.readline(); f.readline()

    E = np.zeros_like(N)
    M = np.zeros_like(N)
    for i in range(int(1e3)):
        sent = f.readline()
        words = sent.split()
        E[i] = words[1]
        M[i] = words[5]

    fig.add_trace(go.Scatter(x=N, y=E, name='Average energy per spin'),row=1,col=1)
    fig.add_trace(go.Scatter(x=N,y=M, name='Average absolute magnetization per spin'),row=2,col=1)

    # check it

    fig.update_xaxes(title_text="Monte Carlo cycles", row=2, col=1)
    fig.update_yaxes(title_text="Energy / J", row=1, col=1)
    fig.update_yaxes(title_text="Magnetization / #", row=2, col=1)

    fig.update_layout(font_family="lmodern",font_size=12)
    fig.update_layout(showlegend=False)
    fig.write_image("testbror.pdf",width=600*1.41,height=600,scale=2)

    fig.show()

    #fig.write_image("fig1.pdf")

def accept():
    N = np.linspace(1, 1e3, int(1e3))

    f1 = open("output1.data", "r")
    f1.readline(); f1.readline()
    f2 = open("output2.data", "r")
    f2.readline(); f2.readline()

    acpt1 = np.zeros_like(N)
    acpt2 = np.zeros_like(N) 
    for i in range(int(1e3)):
        sent1 = f1.readline()
        words1 = sent1.split()
        acpt1[i] = words1[6]

        sent2 = f2.readline()
        words2 = sent2.split()
        acpt2[i] = words2[6]

    fig = make_subplots(rows=2, cols=1, subplot_titles=("Accepted configurations at temperature T=1kT/J","Accepted configurations at temperature T=2.4kT/J"), shared_xaxes=True, vertical_spacing=0.1)

    fig.add_trace(go.Scatter(x=N, y=acpt1, name='Accepted configurations at temperatur T=1kT/J'),row=1,col=1)
    fig.add_trace(go.Scatter(x=N,y=acpt2, name='Accepted configurations at temperatur T=2.4kT/J'),row=2,col=1)

    fig.update_xaxes(title_text="Monte Carlo cycles", row=2, col=1)
    fig.update_yaxes(title_text="Accepted configurations ", row=1, col=1)
    fig.update_yaxes(title_text="Accepted configurations ", row=2, col=1)

    fig.update_layout(font_family="lmodern",font_size=12)
    fig.update_layout(showlegend=False)
    #fig.write_image("testbror.pdf",width=600*1.41,height=600,scale=2)

    fig.show()

if __name__ == "__main__":
    p4a_analytical()
    p4c_comapre_nmc()
    p4e_pde()
    #p4f_many_spin()
    #accept()
    #p4g()
