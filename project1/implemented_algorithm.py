import numpy as np

#hvordan sette andre randbetingelser?

n=11
h = 1/(n+1)

a = np.zeros((n-1,))-1
b = np.zeros((n,))+2
c = np.zeros((n-1,))-1

d = np.zeros((n,))
l = np.zeros((n-1,))


x = np.linspace(0,1,n)#np.zeros((n,))
f = 100*np.exp(-10*x)*h**2
print(h,f)
sol = 1 - (1-np.exp(-10))*x-np.exp(-10*x)


#set f


#find LU-decomp
d[0] = b[0]

for i in range(1,n):
    l[i-1] = a[i-1]/d[i-1]
    d[i] = b[i-1] - l[i-1]*c[i-1]

print(l,d)

"""
L = np.zeros((n,n))
U = np.zeros((n,n))

for i in range(n):
    L[i,i] = 1
    U[i,i] = d[i]

for i in range(n-1):
    L[i+1,i] = l[i]
    U[i,i+1] = c[i]

print(L,"\n\n",U)
print("\n\n", np.matmul(L,U))
"""
#find u:
u = np.zeros((n,))
v = np.zeros((n,))

u[0] = f[0]

for i in range(1,n):
    u[i] = f[i] - l[i-1]*u[i-1]


u[-1] = u[-1]/d[-1]
v[-1] = u[-1]

for i in range(n-2,1,-1):
    v[i] = (u[i] - v[i+1]*c[i])/d[i]


print(v)

#import matplotlib.pyplot as plt
#plt.plot(v)
#plt.plot(sol)
##plt.plot(f)
#plt.legend(["algo","anal"])
#plt.show()

print(np.sum(np.abs(sol-v))/n)