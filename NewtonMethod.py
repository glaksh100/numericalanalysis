from __future__ import division
import numpy as np
import numpy.linalg as la
def newton_method(f,df,init,tol):
    x=[]
    x.append(init)
    e=[]
    k=0
    a=0
    tol1=1
    while(la.norm(tol1)>=tol):
        xnew=x[k]-f(x[k])/df(x[k])
        x.append(xnew)
        
        tol1=x[k+1]-x[k]
        #print tol1
        a=x[k+1]
        k+=1
    n=len(x)
    
    e.append(la.norm(a-x[0]))
    i=1
    while(i<n-2):
        enew=la.norm(a-x[i])
        e.append(enew)
        c1=e[i]/(e[i-1]**1)
        c2=e[i]/(e[i-1]**2)
        #print i,x[i],c1,c2
        i+=1
    return a
def f(x):
    return x**2-1
def df(x):
    return 2*x
x1=newton_method(f,df,10e6,1e-10)

def g(x):
    return (x-1)**4
def dg(x):
    return 4*(x-1)**3
x2=newton_method(g,dg,10,1e-10)


def h(x):
    return x-np.cos(x)
def dh(x):
    return 1+np.sin(x)
x3=newton_method(h,dh,1,1e-10)