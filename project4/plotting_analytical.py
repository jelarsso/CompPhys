import plotly.graph_objects as go
import numpy as np

N = np.linspace(1, 1e2, int(1e2))
print(N)

f = open("output1.data", "r")
f.readline(); f.readline()

E = np.zeros_like(N)
M = np.zeros_like(N)
for i in range(int(1e2)):
    sent = f.readline()
    words = sent.split()
    E[i] = words[1]
    M[i] = words[3]

fig = go.Figure(data = go.Scatter(x=N, y=(E)))
fig.show()
fig = go.Figure(data = go.Scatter(x=N, y=(M)))
fig.show()

#fig.write_image("fig1.pdf")
