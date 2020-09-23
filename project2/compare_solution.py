import numpy as np
import matplotlib.pyplot as plt
from IPython import embed

data  = np.loadtxt("data.out")

solved_eigvals = data[0,:]
sort_keys = np.argsort(solved_eigvals)
solved_eigvals = solved_eigvals[sort_keys]
solved_eigvectors = data[1:,:]

print(solved_eigvals)

N = np.size(solved_eigvals)

h = 1/N
d = 2/h**2
a = -1/h**2

j = np.linspace(1,N-1,N,dtype=np.int64)

anal_eigvals = d + 2*a*np.cos(j*np.pi/N)
anal_eigvector = np.zeros_like(solved_eigvectors)

for i in range(N):
    anal_eigvector[i] = np.sin(i*j*np.pi/N)


plt.plot(anal_eigvals)
plt.plot(solved_eigvals)
plt.show()



plt.plot(-(solved_eigvectors[sort_keys][2,:])/np.linalg.norm((solved_eigvectors[sort_keys][2,:])))
plt.plot(anal_eigvector[2,:]/np.linalg.norm(anal_eigvector[2,:]))
plt.show()
