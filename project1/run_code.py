import numpy as np
import matplotlib.pyplot as plt
import subprocess as sb
import time


def read_file(filename):
    """
    Read the file containing the result of the c++ programs.
    """
    data = []
    with open(filename) as fl:
        for l in fl:
            data.append(float(l))
    data = np.asarray(data)
    return data

def analytical_result(n):
    x = np.linspace(0,1,n,endpoint=True)
    return x,1 - (1-np.exp(-10))*x-np.exp(-10*x)

"""
Each of the functions corresponds to the tasks in project 1.
"""


def func_1b():
    for n in [10,100,1000]:
        sb.run(["./general.out",str(n)])
        calc_result = read_file("general_output.data")
        xaxis,analy_result = analytical_result(calc_result.size)
        plt.plot(xaxis,calc_result)
    plt.plot(xaxis,analy_result)
    
    plt.title("Solution of the linear system using the general algorithm.")
    plt.xlabel("x")
    plt.ylabel("$f_{sol}(x)$")
    plt.legend(["$n=10$","$n=100$","$n=1000$","Analytical"])
    plt.show()


def func_1c(reps=10):
    N=range(1,8)
    general_cpu_time = []
    special_cpu_time = []
    for n in N:
        cpu_time=0
        for i in range(reps):
            start_time = time.time_ns()
            sb.run(["./general.out",str(10**n)])
            end_time = time.time_ns()
            cpu_time += end_time-start_time
        cpu_time/=reps*10**9
        general_cpu_time.append(cpu_time)
        
    for n in N:
        cpu_time=0
        for i in range(reps):
            start_time = time.time_ns()
            sb.run(["./special.out",str(10**n)])
            end_time = time.time_ns()
            cpu_time += end_time-start_time
        cpu_time/=reps*10**9
        special_cpu_time.append(cpu_time)
    
    special_cpu_time = np.asarray(special_cpu_time)
    general_cpu_time = np.asarray(general_cpu_time)
    plt.plot(N,np.log10(general_cpu_time),"rx")
    plt.plot(N,np.log10(special_cpu_time),"bx")
    #plt.plot(N,general_cpu_time-special_cpu_time)
    plt.title("CPU time for general algorithm vs specialized algorithm")
    plt.xlabel("log(n)")
    plt.ylabel("log(CPU time [s])")
    plt.legend(["General algorithm", "Special algorithm"])
    plt.show()
    

def func_1d():
    print("relative errors for special algorithm")
    for n in [10**i for i in range(1,8)]:
        sb.run(["./special.out",str(n)])
        calc_result = read_file("special_output.data")
        xaxis,analy_result = analytical_result(calc_result.size)
        error = np.log10(np.abs((calc_result[1:-1]-analy_result[1:-1])/analy_result[1:-1])) #exclude endpoints since they are zeros
        print(f"log(n) = {np.log10(n)} : {np.max(error)}")
        plt.plot(np.log10(n),np.log10(np.max(error)),"rx")
    plt.title("Max relative error for the specialized algorithm")
    plt.xlabel("log(n)")
    plt.ylabel("log(rel. error)")
    plt.show()



def func_1e(reps=10):
    N=range(1,5)
    print("n | special | arma")
    for n in N:
        arma_cpu_time=0
        for i in range(reps):
            start_time = time.time_ns()
            sb.run(["./arma.out",str(10**n)])
            end_time = time.time_ns()
            arma_cpu_time += end_time-start_time
        arma_cpu_time/=reps*10**9

    
        special_cpu_time=0
        for i in range(reps):
            start_time = time.time_ns()
            sb.run(["./special.out",str(10**n)])
            end_time = time.time_ns()
            special_cpu_time += end_time-start_time
        special_cpu_time/=reps*10**9

        print(f"log(n) = {n} | {special_cpu_time} | {arma_cpu_time} ")
    


print("1b")
func_1b()
print("1c")
func_1c(reps=10)
print("1d")
func_1d()
print("1e")
func_1e()
