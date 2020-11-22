import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

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
