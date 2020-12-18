import numpy as np
import plotly.graph_objects as go
import plotly.subplots as ps
import plotly.express as px
import subprocess as sb
import pandas as pd
from IPython import embed


def read_file(filename):
    # Reads the dumpfile
    return np.loadtxt(filename)


def animate(filename):
    # Animates the one-dimensional solution
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
    # plots all solutions in one dimension with different dx
    for dx in [0.1,0.01]:
        T = 0.1
        dt = 0.5*dx**2
        nt = int(round(T/dt))
        sb.run(["./p5c",str(dx),str(nt),str(.5)])
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
        #fig.show()
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x,y=x+data_f[nt//10,:],mode="lines",name="Forward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=x+data_b[nt//10,:],mode="lines",name="Backward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=x+data_cn[nt//10,:],mode="lines",name="Crank-Nicholson"))
        fig.add_trace(go.Scatter(x=x,y=x+data_an[nt//10,:],mode="lines",name="Analytical"))
        fig.update_xaxes(title="x")
        fig.update_yaxes(title="u(T,x)")
        fig.update_layout(font_family="lmodern",title_text=f"Solutions after T={nt//10} timesteps, dx = {dx}, alpha = 0.5",font_size=12)
        fig.write_image(f"p5c_comparisons_short_dx{dx}.pdf",width=600*1.41,height=600,scale=2)
        #fig.show()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x,y=np.abs(data_an[nt//10,:] - data_f[nt//10,:]),mode="lines",name="Forward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=np.abs(data_an[nt//10,:] - data_b[nt//10,:]),mode="lines",name="Backward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=np.abs(data_an[nt//10,:] - data_cn[nt//10,:]),mode="lines",name="Crank-Nicholson"))
        fig.update_xaxes(title="x")
        fig.update_yaxes(title="abs(u(T,x) - u'(T,x))")
        fig.update_layout(font_family="lmodern",title_text=f"Difference from analytical after T={nt//10} timesteps, dx = {dx}, alpha = 0.5",font_size=12)
        fig.write_image(f"p5c_compare_anal_direct_short{dx}.pdf",width=600*1.41,height=600,scale=2)
        #fig.show()

        fig = go.Figure()
        fig.add_trace(go.Scatter(x=x,y=np.abs(data_an[-1,:] - data_f[-1,:]),mode="lines",name="Forward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=np.abs(data_an[-1,:] - data_b[-1,:]),mode="lines",name="Backward-Euler"))
        fig.add_trace(go.Scatter(x=x,y=np.abs(data_an[-1,:] - data_cn[-1,:]),mode="lines",name="Crank-Nicholson"))
        fig.update_xaxes(title="x")
        fig.update_yaxes(title="abs(u(T,x))")
        fig.update_layout(font_family="lmodern",title_text=f"Difference from analytical after T={nt} timesteps, dx = {dx}, alpha = 0.5",font_size=12)
        fig.write_image(f"p5c_compare_anal_direct_long{dx}.pdf",width=600*1.41,height=600,scale=2)
        #fig.show()


        
    
def animate_2d():
    # Animates two dimensional solutions
    data = np.loadtxt("fe2d.data")
    nt = 1001
    nx = 101
    data = data.reshape((nt,nx,nx))
    

    x = np.linspace(0,1,nx)
    y = np.linspace(0,1,nx)
    z = data[0,:,:]
    
    fig = go.Figure(
    data=[go.Surface(z=z,x=x, y=y)],
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
    frames=[go.Frame(data=[go.Surface(z=data[i,:,:],x=x, y=y)],layout=go.Layout(title_text=f"Frame {i+1}/{data[:,0,0].size}")) for i in range(0,data.shape[0],10)])
    fig.show()



def p5f():
    #make plots of the two-dimensional forward euler solutions and errors.
    sb.run(["./p5f",str(0.01),str(1000),str(0.25)])
    nt = 1000+1
    nx = 100+1

    data_f = np.loadtxt("fe2d.data")
    data_f = data_f.reshape((nt,nx,nx))

    data_a = np.loadtxt("analytical2d.data")
    data_a = data_a.reshape((nt,nx,nx))

    x = np.linspace(0,1,nx)
    xx,yy = np.meshgrid(x,x)

    fig = ps.make_subplots(rows=1,cols=2,subplot_titles=("Forward-Euler", "Analytical"))
    fig.append_trace(go.Contour(z=yy+data_f[-1,:,:],x=x,y=x,coloraxis="coloraxis"),row=1,col=1)
    fig.append_trace(go.Contour(z=yy+data_a[-1,:,:],x=x,y=x,coloraxis="coloraxis"),row=1,col=2)
    fig.update_yaxes(title="x")
    fig.update_xaxes(title="y")
    
    fig.update_layout(coloraxis=dict(colorscale='Bluered_r'),coloraxis_colorbar=dict(title="u(x,y)"), showlegend=False)
    fig.update_layout(font_family="lmodern",title_text=f"Solutions after T={nt-1} timesteps, dx = {0.01}, alpha = 0.25",font_size=12)
    fig.write_image(f"p5f_2d_sol_long.pdf",width=600*1.41,height=600,scale=2)
    fig.show()
    
    stds = []
    nxs = [ 4**i - 1 for i in range(1,5)]#np.linspace(10,50,20) #[0.1,0.01,0.001]
    nt = 1000
    dt = 0.001/nt
    for nx in nxs:
        dx = 1/(nx+1)
        alpha = dt/dx/dx
        out = sb.run(["./p5f",str(dx),str(nt),str(alpha)],capture_output=True)
        nx = int(str(out.stdout).split(" ")[6])
        print("nx = ", nx+1, " nt = ", nt, " dx = ", dx, " dt = ", dt, " alpha = ", alpha)

        data_f = np.loadtxt("fe2d.data")
        data_f = data_f.reshape((nt+1,nx+1,nx+1))

        data_a = np.loadtxt("analytical2d.data")
        data_a = data_a.reshape((nt+1,nx+1,nx+1))

        #x = np.linspace(0,1,nx-1,endpoint=False)
        x = np.arange(0,1+1e-10,dx)
        xx,yy = np.meshgrid(x,x)

        print("x ", x.shape)
        print("x = ", x)
        print("x slice ",x[::(nx)//4][1:-1])
        print("x slice ",x[::(nx)//4][1:-1].shape)
        print("data : ", data_f[-1,::(nx)//4,::(nx)//4][1:-1,1:-1])
        #print(data_a[0,1:-1:nx//4,1:-1:nx//4].shape)

        stds.append(np.sum((data_a[-1,::(nx)//4,::(nx)//4][1:-1,1:-1]-data_f[-1,::(nx)//4,::(nx)//4][1:-1,1:-1])**2))
    
    #these values are from a previous run of the above for loop:
    stds = [0.0019045588662245931, 1.1071996881704011e-07, 1.845193450238234e-12, 3.98819614687517e-15]
    dxs = [0.25,       0.0625,     0.015625,   0.00390625]

    print(stds)
    print(dxs)
    import scipy.stats as ss
    print(ss.linregress(np.log10(dxs),np.log10(stds)))
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=np.log10(dxs),y=np.log10(np.asarray(stds))))
    fig.update_xaxes(title="log(dx)")
    fig.update_yaxes(title="log(standard deviation)")
    fig.update_layout(font_family="lmodern",title_text="The standard deviation of the difference between the analytical and numerical solution.",font_size=12)
    fig.write_image("p5f_std_diff_dxs_1.pdf",width=600*1.41,height=600,scale=2)
    fig.show()
    
    
    stds = []
    dts = np.linspace(0.0025,0.001,10) #[0.1,0.01,0.001]
    nx = 10
    #alpha = 0.25
    for dt in dts:
        dx = 1/nx #np.sqrt(dt)/alpha
        nt = int(1/dt)
        alpha = dt/dx/dx
        out = sb.run(["./p5f",str(dx),str(nt),str(alpha)],capture_output=True)
        nx = int(str(out.stdout).split(" ")[6])
        print("nx = ", nx, " nt = ", nt, " dx = ", dx, " dt = ", dt, " alpha = ", alpha)


        data_f = np.loadtxt("fe2d.data")
        data_f = data_f.reshape((nt+1,nx+1,nx+1))

        data_a = np.loadtxt("analytical2d.data")
        data_a = data_a.reshape((nt+1,nx+1,nx+1))

        x = np.linspace(0,1,nx+1)
        stds.append(np.std(data_a-data_f))
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=dts,y=stds))
    fig.update_xaxes(title="dt")
    fig.update_yaxes(title="standard deviation")
    fig.update_layout(font_family="lmodern",title_text="The standard deviation of the difference between the analytical and numerical solution.",font_size=12)
    fig.write_image("p5f_std_diff_dts.pdf",width=600*1.41,height=600,scale=2)
    fig.show()
    
    

    

def show_differences_litho():
    #plot the quantities of interest for the lithosphere simulation.    
    nt = 401
    nx = 80
    ny = 159
    islice = 80
    print(nt*nx*ny)
    sb.run(["make","lithos"])
    sb.run(["./litho", "0.01", "1"])

    data = np.loadtxt("litho_enriched.data")
    data_enriched = data.reshape((nt,nx,ny))
    data = np.loadtxt("litho_no_Q.data")
    data_steady = data.reshape((nt,nx,ny))
    data = np.loadtxt("litho_Q_pb.data")
    data_steay2 = data.reshape((nt,nx,ny))

    x = np.linspace(0,1,nx)
    y = np.linspace(0,2,ny)

    fig = ps.make_subplots(rows=2,cols=1,shared_xaxes=True,vertical_spacing=0.02)
    fig.append_trace(go.Scatter(x=x,y=data_steady[-1,:,islice],mode="lines",name="Steady state"),row=1,col=1)
    fig.append_trace(go.Scatter(x=x,y=data_steay2[-1,:,islice],mode="lines",name="Pre-enriched",line=dict(color="#EF553B")),row=1,col=1)
    fig.append_trace(go.Scatter(x=x,y=data_enriched[-1,:,islice],mode="lines",name="Enriched",line=dict(color="#00CC96")),row=1,col=1)

    fig.update_yaxes(title="Temperature / C",row=1,col=1)
    
    
    fig.append_trace(go.Scatter(x=x,y=(data_steay2[-1,:,islice]-data_steady[-1,:,islice]),mode="lines",line=dict(color="#EF553B"),showlegend=False),row=2,col=1)
    fig.append_trace(go.Scatter(x=x,y=(data_enriched[-1,:,islice]-data_steady[-1,:,islice]),mode="lines",line=dict(color="#00CC96"),showlegend=False),row=2,col=1)
    fig.update_yaxes(title="Absolute Difference Temperature / C", row=2,col=1)
    fig.update_xaxes(title="Depth / L",row=2,col=1)
    
    fig.update_layout(font_family="lmodern",title_text=f"Lithosphere slice through center after 1 Gy",font_size=12)
    fig.write_image(f"litho_temps.pdf",width=600*1.41,height=600,scale=2)
    fig.show()
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x,y=(100*data_steay2[-1,:,islice]-data_steady[-1,:,islice])/data_steady[-1,:,islice],mode="lines",name="Relative pre-enriched"))
    fig.add_trace(go.Scatter(x=x,y=(100*data_enriched[-1,:,islice]-data_steady[-1,:,islice])/data_steady[-1,:,islice],mode="lines",name="Relative enriched"))
    fig.update_xaxes(title="Depth / L")
    fig.update_yaxes(title="Relative Difference / %")
    fig.update_layout(font_family="lmodern",title_text=f"Relative Difference from the steady state after 1 Gy",font_size=12)
    fig.write_image(f"litho_relative_differences.pdf",width=600*1.41,height=600,scale=2)
    fig.show()
    
    fig = ps.make_subplots(rows=1,cols=2,subplot_titles=("After 0.25 Gy","After 1 Gy"))
    fig.append_trace(go.Contour(x=x,y=y,z=data_enriched[100,:,:],coloraxis="coloraxis"),row=1,col=1)
    fig.append_trace(go.Contour(x=x,y=y,z=data_enriched[-1,:,:],coloraxis="coloraxis"),row=1,col=2)
    fig.update_layout(font_family="lmodern",title_text=f"Contour plot of the temperature in the Lithosphere",font_size=12)
    fig.update_xaxes(title_text="y / L", row=1, col=1)
    fig.update_yaxes(title_text="x / L", row=1, col=1)
    fig.update_xaxes(title_text="y / L", row=1, col=2)
    fig.update_yaxes(title_text="x / L", row=1, col=2)
    fig.update_layout(coloraxis_colorbar=dict(title="Temperature / C"))
    fig.write_image(f"litho_temp_contour.pdf",width=600*1.41,height=600,scale=2)
    fig.show()
    
    
    ds = data_steady[-1,:,islice]
    ds2 = data_steay2[-1,:,islice]
    de = data_enriched[:,:,islice]

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


def animate_contour():
    #Animate a contour plot
    nt = 401
    nx = 80
    ny = 159
    print(nt*nx*ny)
    sb.run(["make","lithos"])
    sb.run(["./litho", "0.01", "1"])
    data = np.loadtxt("litho_enriched.data")
    data_enriched = data.reshape((nt,nx,ny))

    x = np.linspace(0,1,nx)
    y = np.linspace(0,2,ny)


    fig = go.Figure(
    data=[go.Contour(x=y, y=x,z=data_enriched[0,:,:])],
    layout=go.Layout(
        xaxis=dict(range=[np.min(y), np.max(y)], autorange=False),
        yaxis=dict(range=[np.min(x), np.max(x)], autorange=False),
        #zaxis=dict(range=[np.min(data_enriched), np.max(data_enriched)], autorange=False),
        title="Start Title",
        updatemenus=[dict(
            type="buttons",
            buttons=[dict(label="Play",
                          method="animate",
                          args=[None,{"frame": {"duration": 1, "redraw": True},"fromcurrent": False, "transition": {"duration": 0}}])])]
    ),
    frames=[go.Frame(data=[go.Contour(x=y, y=x, z=data_enriched[i,:,:])],layout=go.Layout(title_text=f"Time-step {i}/{nt}")) for i in range(0,nt)])
    fig.show()




def p5d_dt_analysis():
    # plots numerical deviations / error as a function of timesteps
    sum_f = []
    sum_b = []
    sum_cn = []
    xx = []    
    for alpha in [.498,  .45, .4, .35, .3, .25, .2, .15, .1]:
        dx = .01
        T = .01
        dt = alpha*dx**2
        nt = int(round(T/dt))
        sb.run(["./p5c",str(dx),str(nt),str(alpha)])
        
        data_f = read_file("forward_euler.data")
        x = np.linspace(0,1,data_f.shape[1])
        data_f = x+read_file("forward_euler.data")
        data_b = x+read_file("backward_euler.data")
        data_cn = x+read_file("cnicholson.data")
        data_an = x+read_file("analytical.data")

        pf = []; pb = []; pcn = []
        for i in range(len(data_f[-1,:])):
            pf.append(abs(data_f[-1,i] - data_an[-1,i])) 
            pb.append(abs(data_b[-1,i] - data_an[-1,i])) 
            pcn.append(abs(data_cn[-1,i] - data_an[-1,i])) 

        sum_f.append(np.std(pf))
        sum_b.append(np.std(pb))
        sum_cn.append(np.std(pcn))
        xx.append(dt)
    
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=xx,y=sum_f,mode="lines+markers",name="Forward-Euler"))
    fig.add_trace(go.Scatter(x=xx,y=sum_b,mode="lines+markers",name="Backward-Euler"))
    fig.add_trace(go.Scatter(x=xx,y=sum_cn,mode="lines+markers",name="Crank-Nicolson"))
    fig.update_xaxes(title=f"Timestep / dt")
    fig.update_yaxes(title=f"Stanard deviation / s")
    fig.update_layout(font_family="lmodern",title_text=f"Standard deviation from analytical solution after T={nt} timesteps, dx = {dx}",font_size=12)
    fig.write_image(f"p5d_best_sol.pdf",width=600*1.41,height=600,scale=2)   
    fig.show()


def p5d_log_log_dt():
    # Log-log plot of deviations as function of timestep
    sum_b = []
    sum_cn = []
    xx = []    
    sb.run(["make","5c"])

    for alpha in [ 30, 25, 20, 15, 10, 1]:
        dx = .01
        T = .01
        dt = alpha*dx**2
        nt = int(round(T/dt))
        sb.run(["./p5c",str(dx),str(nt),str(alpha)])
        
        data_f = read_file("forward_euler.data")
        x = np.linspace(0,1,data_f.shape[1])
        data_f = x+read_file("forward_euler.data")
        data_b = x+read_file("backward_euler.data")
        data_cn = x+read_file("cnicholson.data")
        data_an = x+read_file("analytical.data")

        pf = []; pb = []; pcn = []
        for i in range(len(data_f[-1,:])):
            pb.append(abs(data_b[-1,i] - data_an[-1,i])) 
            pcn.append(abs(data_cn[-1,i] - data_an[-1,i])) 

        sum_b.append(np.std(pb))
        sum_cn.append(np.std(pcn))
        xx.append(dt)
    
    lgx = np.log10(xx)
    lgb = np.log10(sum_b)
    lgcn = np.log10(sum_cn)

    fit1 = np.polyfit(lgx, lgb, 1, full=True)
    line1 = fit1[0][0]*lgx+fit1[0][1]

    fit2 = np.polyfit(lgx, lgcn, 1, full=True)
    line2 = fit2[0][0]*lgx+fit2[0][1]

    print('Stigning for BE:', fit1[0][0], '+-', fit1[1][0])
    print('Stigning for CN:', fit2[0][0], '+-', fit2[1][0])
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lgx,y=lgb,mode="lines+markers",name="Backward-Euler",line = dict(color='royalblue', width=2)))
    fig.add_trace(go.Scatter(x=lgx,y=lgcn,mode="lines+markers",name="Crank-Nicolson",line = dict(color='firebrick', width=2)))
    fig.add_trace(go.Scatter(x=lgx,y=line1,name="Backward-Euler reg",line = dict(color='royalblue', width=2, dash='dash')))
    fig.add_trace(go.Scatter(x=lgx,y=line2,name="Crank-Nicolson reg",line = dict(color='firebrick', width=2, dash='dash')))    
    fig.update_xaxes(title=f"Timestep / log(dt)")
    fig.update_yaxes(title=f"Standard deviation / log(s)")
    fig.update_layout(font_family="lmodern",title_text=f"Standard deviation from analytical solution, dx = {dx} (LOG-LOG)",font_size=12)
    fig.write_image(f"p5d_error_rate.pdf",width=600*1.41,height=600,scale=2)   
    fig.show()


def p5d_log_log_dx():
    # Log-log plot of deviations as function of steplength
    sum_b = []
    sum_cn = []
    xx = []    
    lens = []
    sb.run(["make","5c"])


    for alpha in [80, 70, 60, 50, 40, 30, 20, 10, 8, 7, 6, 5, 3, 2, 1]:
        dt = .01
        T = 1
        dx = np.sqrt(dt/alpha)
        nt = int(round(T/dt))
        sb.run(["./p5c",str(dx),str(nt),str(alpha)])
        
        data_f = read_file("forward_euler.data")
        x = np.linspace(0,1,data_f.shape[1])
        data_f = x+read_file("forward_euler.data")
        data_b = x+read_file("backward_euler.data")
        data_cn = x+read_file("cnicholson.data")
        data_an = x+read_file("analytical.data")


        print('nx = ',np.sqrt(alpha/dt))
        pf = []; pb = []; pcn = []
        for i in range(len(data_f[-1,::alpha])):
            pb.append(abs(data_b[-1,i] - data_an[-1,i])) 
            pcn.append(abs(data_cn[-1,i] - data_an[-1,i])) 
        
        sum_b.append(np.std(pb))
        sum_cn.append(np.std(pcn))
        xx.append(dx)
        lens.append(len(data_f[-1]))

    sum_b = np.array(sum_b)
    sum_cn = np.array(sum_cn)
    xx = np.array(xx)    
    lens = np.array(lens)
    sum_b = sum_b/lens
    sum_cn = sum_cn/lens

    lgx = np.log10(xx)
    lgb = np.log10(sum_b)
    lgcn = np.log10(sum_cn)

    fit1 = np.polyfit(lgx, lgb, 1, full=True)
    line1 = fit1[0][0]*lgx+fit1[0][1]
    print('residuals 1: ', fit1[1][0])

    fit2 = np.polyfit(lgx, lgcn, 1, full=True)
    line2 = fit2[0][0]*lgx+fit2[0][1]
    print('residuals 2: ', fit2[1][0])

    print('Stigning for BE:', fit1[0][0])
    print('Stigning for CN:', fit2[0][0])
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lgx,y=lgb,mode="lines+markers",name="Backward-Euler",line = dict(color='royalblue', width=2),))
    fig.add_trace(go.Scatter(x=lgx,y=lgcn,mode="lines+markers",name="Crank-Nicolson",line = dict(color='firebrick', width=2)))
    fig.add_trace(go.Scatter(x=lgx,y=line1,
        name="Backward-Euler reg",
        line = dict(color='royalblue', width=2, dash='dash')))
    fig.add_trace(go.Scatter(x=lgx,y=line2,
        name="Crank-Nicolson reg",
        line = dict(color='firebrick', width=2, dash='dash'))) 
  
    fig.update_xaxes(title=f"Steplength / log(dx)")
    fig.update_yaxes(title=f"Standard deviation / log(s)")
    fig.update_layout(font_family="lmodern",title_text=f"Standard deviation from analytical solution after T={nt} timesteps, dt = {dt} (LOG-LOG)",font_size=12)
    fig.write_image(f"p5d_error_rate_dx.pdf",width=600*1.41,height=600,scale=2)   
    fig.show()


if __name__=="__main__":
    #Choose one!:

    #animate("cnicholson.data")
    compare_p5c()
    #animate_2d()
    #show_differences_litho()
    #p5d_dt_analysis()
    #p5d_log_log_dt
    #p5d_log_log_dx
