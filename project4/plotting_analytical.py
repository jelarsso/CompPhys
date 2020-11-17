import plotly.graph_objects as go
import numpy as np


beta = 1/np.linspace(1,100,100)
Z = (12 + 4*np.cosh(8*beta))
mean_M = (16+8*np.exp(8*beta))/Z
mean_M2 = 32*(1+np.exp(8*beta))/Z


fig = go.Figure(data = go.Scatter(x=1/beta, y=(mean_M2 - mean_M**2)*beta))
fig.show()