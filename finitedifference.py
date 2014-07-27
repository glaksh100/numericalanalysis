from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
import scipy.special as sps 

def finite_difference(f,a,b,h):
    n=(b-a)/h+1
    x=np.linspace(a,b,n)
    print len(x)
    n=len(x)
    df=[]
    e=[]
    for i in range(0,n):
        d=(4*f(x[i]+h)-3*f(x[i])-f(x[i]+2*h))/(2*h)
        er=(d-np.cos(x[i]))
        df.append(d)
        e.append(er)
        e1=np.array(e)
        maxerr=e1.max()
    return maxerr
def f(x):
    return np.sin(x)
k=3
h=[]
for i in range(3,15):
    h.append(2**-i)
m=len(h)

maxe=[]
for i in range(0,m):
    val=finite_difference(f,-1,1,h[i])
    maxe.append(val)
pt.plot(h,np.log(maxe))
pt.xlabel("h")
pt.ylabel("E")
pt.title("Plot of E vs h for one sided difference")
pt.show()