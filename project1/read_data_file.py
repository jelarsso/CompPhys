import numpy as np
import matplotlib.pyplot as plt



data = []
with open("output.data") as fl:
    for l in fl:
        data.append(float(l))

data = np.asarray(data)

x = np.linspace(0,1,len(data))

analytisk = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)



plt.plot(data)
plt.plot(analytisk)
plt.show()