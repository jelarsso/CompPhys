import numpy as np
import matplotlib.pyplot as plt
from IPython import embed


data = np.loadtxt("output.data")
v = data.reshape((10000,2,3))

earthpos = np.zeros((3,))

#earthpos = v[:,1,:] - v[:,0,:]

#for i in range(1):
plt.plot(v[:,1,0],v[:,1,1])
plt.plot(v[:,0,0],v[:,0,1],"rx")
#print(v)
#plt.plot(np.linalg.norm(earthpos,axis=1))
plt.show()