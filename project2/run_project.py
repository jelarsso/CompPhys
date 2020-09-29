import numpy as np
import matplotlib.pyplot as plt
import subprocess as sb
from time import time_ns



def e2b_number_of_iterations():
    n_vals = []
    iterations_needed = []
    legend = []
    fig,ax = plt.subplots(1,1)
    for i,tol in enumerate([1e-6]):
        for n in range(200):
            n_vals.append(n)
            result = sb.run(["./main.out", str(n), str(tol)],capture_output=True)
            iterations_needed.append(int(result.stdout))
        ax.plot(n_vals,iterations_needed)
        legend.append(f"log(tol)={np.log10(tol):.1f}")
        iterations_needed = []
        n_vals = []
    
    n_vals = range(100)
    ax.fill_between(n_vals,3*np.asarray(n_vals)**2,5*np.asarray(n_vals)**2,where=True,facecolor="red",alpha=0.5)
    fig.suptitle("Convergence rate for Jacobi rotation method")
    ax.set_xlabel("size of matrix")
    ax.set_ylabel("number of iterations needed")
    ax.legend(legend+["theoretical"],loc="upper left")
    plt.show()


def e2b_armadillo_vs_own():
    arma_times = []
    own_times = []

    for i in range(1,200):
        print(i)
        start_time = time_ns()
        sb.run(["./arma.out",str(i)])
        arma_times.append((time_ns()-start_time)/1e9)

        start_time = time_ns()
        sb.run(["./main.out",str(i),str(1e-8)])
        own_times.append((time_ns()-start_time)/1e9)
    
    plt.plot(arma_times)
    plt.plot(own_times)
    plt.xlabel("size of matrix")
    plt.ylabel("elapsed time [s]")
    plt.show()

def e2c_quantum():
    data = np.loadtxt("eigen_data.out",skiprows=1)
    sortkeys = np.argsort(data[0,:])
    evals = data[0,:][sortkeys]
    evecs = data[1:,:][sortkeys]

    for i in range(10):
        print(f"Eigenvalues {i+1} :  {evals[i]}")


def e2e_two_electron():
    omega_r = [0.01]#,0.5,1,5]
    first_eig = []
    for omega in omega_r:
        sb.run(["./quantum_two.out",str(400),str(1e-8),str(omega)])
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


#e2b_armadillo_vs_own()
#e2b_number_of_iterations()
e2c_quantum()
#e2e_two_electron()