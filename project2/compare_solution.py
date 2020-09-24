import numpy as np
import matplotlib.pyplot as plt
from IPython import embed

def analytical_solution(N):
    h = 1/N
    d = 2/h**2
    a = -1/h**2

    j = np.linspace(1,N,N,dtype=np.int64)
    print(j)


    anal_eigvals = d + 2*a*np.cos(j*np.pi/(N+1))
    anal_eigvector = np.zeros_like(solved_eigvectors)

    for i in range(1,N+1):
        anal_eigvector[i-1] = np.sin(i*j*np.pi/(N+1))
    return anal_eigvals,anal_eigvector

data  = np.loadtxt("data.out",skiprows=1)

solved_eigvals = data[0,:]
sort_keys = np.argsort(solved_eigvals)
solved_eigvals = solved_eigvals[sort_keys]
solved_eigvectors = data[1:,:][sort_keys]

print(solved_eigvals)

N = np.size(solved_eigvals)

anal_eigvals,analytical_solution = analytical_solution(N)

plt.plot(anal_eigvals)
plt.plot(solved_eigvals)
plt.show()



plt.plot(-(solved_eigvectors[6,:])/np.linalg.norm((solved_eigvectors[6,:])))
plt.plot(anal_eigvector[6,:]/np.linalg.norm(anal_eigvector[6,:]))
plt.legend(["solv","anal"])
plt.show()

