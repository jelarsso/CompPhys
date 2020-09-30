import numpy as np
import matplotlib.pyplot as plt
import subprocess as sb
from IPython import embed


def analytical_solution(N):
    """
        N: int: positive integer. Finds the N first eigenvalues and eigenvectors
    Return the analytical solutions to the buckling beam problem specified in the article.
        analytical_eigvals: array of dimension (1,N): the analytically calculated eigenvalues
        analytical_eigvector: array of dimension (N,N): the matrix of the basis eigenvectors. The eigenvectors are analytical_eigvector[i,:]
    """
    h = 1/N
    d = 2/h**2
    a = -1/h**2

    j = np.linspace(1,N,N,dtype=np.int64)
    print(j)


    analytical_eigvals = d + 2*a*np.cos(j*np.pi/(N+1))
    analytical_eigvector = np.zeros((N,N))

    for i in range(1,N+1):
        analytical_eigvector[i-1] = np.sin(i*j*np.pi/(N+1))
    return analytical_eigvals,analytical_eigvector


sb.run(["./main.o", "100","1e-8"])
sb.run(["./arma.o", "100"])

data  = np.loadtxt("data.out",skiprows=1)

jacobi_rot_eigvals = data[0,:]
sort_keys = np.argsort(jacobi_rot_eigvals)
jacobi_rot_eigvals = jacobi_rot_eigvals[sort_keys]
jacobi_rot_eigvectors = data[1:,:][sort_keys]



arma_data  = np.loadtxt("arma_data.out")

arma_eigvals = arma_data[0,:]
arma_sort_keys = np.argsort(arma_eigvals)
arma_eigvals = arma_eigvals[arma_sort_keys]
arma_eigvectors = data[1:,:][arma_sort_keys]

N = np.size(jacobi_rot_eigvals)

analytical_eigvals,analytical_eigvector = analytical_solution(N)
embed()
fig,ax = plt.subplots(2,1)
fig.suptitle("Comparison between the analytical, jacobi rotation and armadillo eigenvalues")
ax[0].set_title("Eigenvalues")
ax[0].set_xlabel("Eigenvalue number")
ax[0].set_ylabel("Eigenvalue")
ax[0].plot(analytical_eigvals)
ax[0].plot(jacobi_rot_eigvals)
ax[0].plot(arma_eigvals)
ax[0].legend(["Analytical", "Jacobi rotation", "Armadillo"])

ax[1].set_title("Relative difference compared to analytically obatined eigenvalues")
ax[1].plot((analytical_eigvals-jacobi_rot_eigvals)/analytical_eigvals)
ax[1].plot((analytical_eigvals-arma_eigvals)/analytical_eigvals)
ax[1].set_xlabel("Eigenvalue number")
ax[1].set_ylabel("Relative error")
ax[1].legend(["Jacobi rotation", "Armadillo"])
plt.show()



jr = np.zeros((N+2))
al = np.zeros((N+2))
jr[1:-1] = (jacobi_rot_eigvectors[0,:])/np.linalg.norm(jacobi_rot_eigvectors[0,:])
al[1:-1] = analytical_eigvector[0,:]/np.linalg.norm(analytical_eigvector[0,:])
plt.plot(jr)
plt.plot(al)
plt.legend(["Jacobi rotation","Analytical"])
plt.title("Eigenvector corresponding to the lowest eigenvalue")
plt.xlabel("Component $x_i$")
plt.ylabel("Value $u_i$")
plt.show()

