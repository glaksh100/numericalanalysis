from __future__ import division
import numpy as np
import numpy.linalg as la
import scipy.optimize as spo
import matplotlib.pyplot as pt
import math as m
y0=np.array([20,51.58,68.73,75.46,74.36,67.09,54.73,37.98,17.28])
t0=np.arange(0,2.25,0.25)
def f(x,t):
    return x[0]+x[1]*t+x[2]*(t**2)+x[3]*m.exp(x[4]*t)
def ff(x,t):
    F=[]
    for i in range(0,len(t)):
        F.append(f(x,t[i]))
    return F
def jacobian(x,t=t0):
    J=[]
    for i in range(0,len(t)):
            val1=-1
            val2=-t[i]
            val3=-t[i]**2
            val4=-m.exp(x[4]*t[i])
            val5=-x[3]*t[i]*m.exp(x[4]*t[i])
            Jtemp=np.array([val1,val2,val3,val4,val5])
            J.append(Jtemp)
    return np.array(J)
def residual(x,y=y0,t=t0):
    r=[]
    for i in range(0,len(t)):
        rtemp=y[i]-f(x,t[i])
        r.append(rtemp)
    return np.array(r)
            

def gauss_newton(x0,y=y0,t=t0):
    tol=1e-14
    tol1=1
    k=0
    x=[]
    x.append(x0)
    while(tol1>=tol and k<500):
            
        Jx=jacobian(x[k],t)
        rx=residual(x[k])
        #print Jx,rx
        s=la.lstsq(Jx,-rx)[0]
        xnew=x[k]+s
        x.append(xnew)
        tol1=la.norm(x[k+1]-x[k])
        k+=1
    print k
    return x[k]        
x0=np.array([0,0,0,0,1])  
ans=gauss_newton(x0)         
g=ff(ans,t0)
    
p1=pt.plot(t0,y0,label='Data')
pt.legend()
pt.plot(t0,g,label='First Starting vector')
pt.legend()
pt.show()
pt.figure(2)
x1=np.array([1,0,0,0,0])  
ans1=gauss_newton(x1)         
g1=ff(ans1,t0)
    
pt.plot(t0,y0,label='Data')
pt.plot(t0,g1,label='Second Starting Vector')
pt.legend()
pt.figure(3)
x2=np.array([1,0,0,1,0])  
ans2=gauss_newton(x2)         
g2=ff(ans2,t0)
    
pt.plot(t0,y0,label='Data')
pt.plot(t0,g2,label='Third Starting Vector')
pt.legend()
pt.show()