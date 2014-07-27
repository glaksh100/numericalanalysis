from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
import scipy.special as sps 

def gauss_quad(f,n,a,b):
    nodes=sps.legendre(n).weights[:,0]
    V=np.zeros([n,n])
    for i in range(0,n):
        for j in range(0,n):
            V[i][j]=nodes[i]**j
    A=V.T
    q=np.zeros([n,1])
    for i in range(0,n):
        q[i]=(b**(i+1) - a**(i+1))/(i+1)
    sol=la.solve(A,q)
    
    sum=0
    for i in range(0,n):
        sum=sum+sol[i]*f(nodes[i])
    return sum

def f(x):
	return np.sin(2*np.pi*x)
def g(x):
        return np.abs(x)
def gauss_quad_large(f,n,a,b):
    nodes=sps.legendre(n).weights[:,0]
    weights=sps.legendre(n).weights[:,1]
    return np.dot(weights,f(nodes))

ans_f=[]
e_f=[]
e_g=[]
ans_g=[]
n= np.linspace(1,6,6)
for i in range(0,6):
    ans1=gauss_quad(f,int(n[i]),-1,1)
    ans_f.append(ans1)
    e_f.append(np.abs(ans1-0))
    ans2=gauss_quad(g,int(n[i]),-1,1)
    ans_g.append(ans2)
    e_g.append(np.abs(ans2-1))
pt.plot(np.log(n),e_f,label='Error for f')
pt.plot(np.log(n),np.log(e_g),label='Error for g')
pt.xlabel("log n")
pt.ylabel("log E")
pt.title("Solving for weights using n=6")
pt.legend()
#print ans_g
ans_f1=[]
e_f1=[]
e_g1=[]
ans_g1=[]
n1= np.linspace(1,100,100)
for i in range(0,100):
    ans11=gauss_quad_large(f,int(n1[i]),-1,1)
    ans_f1.append(ans11)
    e_f1.append(np.abs(ans11-0))
    ans21=gauss_quad_large(g,int(n1[i]),-1,1)
    ans_g1.append(ans21)
    e_g1.append(np.abs(ans21-1))
pt.figure(2)
pt.plot(np.log(n1),e_f1,label='Error for f')
pt.plot(np.log(n1),np.log(e_g1),label='Error for g')
pt.xlabel("log n")
pt.ylabel("log E")
pt.title("Using scipy weights for n=100")
p_g=np.log(e_g[4]/e_g[2])/np.log(n[2]/n[4])
print p_g
pt.legend()
pt.show()
