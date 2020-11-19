import plotly.graph_objects as go
from plotly.subplots import make_subplots
import numpy as np

N = np.linspace(1, 1e3, int(1e3))
print(N)

f = open("output1.data", "r")
f.readline(); f.readline()

E = np.zeros_like(N)
M = np.zeros_like(N)
for i in range(int(1e3)):
    sent = f.readline()
    words = sent.split()
    E[i] = words[1]
    M[i] = words[3]

fig = make_subplots(rows=2, cols=1, subplot_titles=("Simulation of 20x20 spins with temperature T=2.4kT/J"," "), shared_xaxes=True, vertical_spacing=0.06)
    
fig.add_trace(go.Scatter(x=N, y=E, name='Average Energy per spin'),row=1,col=1)
fig.add_trace(go.Scatter(x=N,y=M, name='Average Magnetization per spin'),row=2,col=1)

fig.update_xaxes(title_text="Monte Carlo cycles", row=2, col=1)
fig.update_yaxes(title_text="Energy / J", row=1, col=1)
fig.update_yaxes(title_text="Magnetization / #", row=2, col=1)

fig.show()

#fig.write_image("fig1.pdf")
