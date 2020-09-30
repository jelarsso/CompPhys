import numpy as np
import matplotlib.pyplot as plt
import subprocess as sb
from time import time_ns
from IPython import embed


def e2b_number_of_iterations():
    """
    Finds the number of iterations needed in the jacobi rotation algorithm, and plots them.
    """
    n_vals = []
    iterations_needed = []
    legend = []
    fig,ax = plt.subplots(1,1)
    for i,tol in enumerate([1e-6,1e-8,1e-16,1e-32]):
        for n in range(200):
            n_vals.append(n)
            result = sb.run(["./main.o", str(n), str(tol)],capture_output=True)
            iterations_needed.append(int(result.stdout))
        ax.plot(n_vals,iterations_needed)
        legend.append(f"log(tol)={np.log10(tol):.1f}")
        iterations_needed = []
        n_vals = []
    
    n_vals = range(200)
    ax.fill_between(n_vals,3*np.asarray(n_vals)**2,5*np.asarray(n_vals)**2,where=True,facecolor="red",alpha=0.5)
    fig.suptitle("Convergence rate for Jacobi rotation method")
    ax.set_xlabel("size of matrix")
    ax.set_ylabel("number of iterations needed")
    ax.legend(legend+["theoretical"],loc="upper left")
    plt.show()



def e2b_armadillo_vs_own():
    """
    Times the jacobi_rotation_algorithm vs the armadillo implementation, and plots the results.
    """
    arma_times = []
    own_times = []

    for i in range(1,200):
        print(i)
        start_time = time_ns()
        sb.run(["./arma.o",str(i)])
        arma_times.append((time_ns()-start_time)/1e9)

        start_time = time_ns()
        sb.run(["./main.o",str(i),str(1e-8)])
        own_times.append((time_ns()-start_time)/1e9)
    
    plt.plot(arma_times)
    plt.plot(own_times)
    plt.legend(["Armadillo implementation", "Jacobi rotation algorithm"])
    plt.xlabel("size of matrix")
    plt.ylabel("elapsed time [s]")
    plt.show()

def e2c_quantum():
    """
    Find the 10 first eigenvalues from the data file eigen.data. To use this, the program corresponding to quantum_dots_two_electron.cpp must be executed first.
    """
    data = np.loadtxt("eigen.data",skiprows=1)
    sortkeys = np.argsort(data[0,:])
    evals = data[0,:][sortkeys]
    evecs = data[1:,:][sortkeys]

    for i in range(10):
        print(f"Eigenvalues {i+1} :  {evals[i]}")


def e2e_two_electron():
    """
    Find and plot the wave function for the lowest eigenvalue for the two electron problem with varying omega_r
    """
    omega_r = [0.01,0.5,1,5]
    first_eig = []
    for omega in omega_r:
        sb.run(["./quant_two.o",str(400),str(1e-8),str(omega)])
        data = np.loadtxt("two_electron.data",skiprows=1)
        sorts = np.argsort(data[0,:])
        first_eig.append(data[1:,:][sorts][0])


    rho = np.linspace(0,100,first_eig[0].size+2)
    for i,omega in enumerate(omega_r):
        print(data[0,:][sorts][i])
        dim = first_eig[0].size + 2
        vals = np.zeros(dim)
        vals[1:-1] = first_eig[i]
        plt.plot(rho,vals/np.linalg.norm(vals))
    
    plt.legend(["$\\omega_r = " + str(o) + "$" for o in omega_r])
    plt.xlabel("$\\rho_i$")
    plt.ylabel("wave function $\\psi(\\rho_i)$")
    plt.title("Wave function for two electrons")
    plt.savefig("figs/wave_two_lowomega.pdf",format="pdf")
    plt.show()

#Choose one (or all)
#e2b_armadillo_vs_own()
#e2b_number_of_iterations()
#e2c_quantum()
# #e2e_two_electron()