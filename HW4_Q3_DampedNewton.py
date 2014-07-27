from __future__ import division
import numpy as np
import numpy.linalg as la
import scipy.optimize as spo
import matplotlib.pyplot as pt
def f(x):
    return (x[0]**2+x[1]-11)**2 + (x[0]+x[1]**2 -7)**2
def df(x):
    a=2*(x[0]**2+x[1]-11)*2*x[0]+2*(x[0]+x[1]**2 -7)
    b=2*(x[0]**2+x[1]-11)+2*(x[0]+x[1]**2 -7)*2*x[1]
    return np.array([a,b])
def hessian(x):
    val1=12*(x[0]**2)+4*x[1]-42
    val2=4*x[0]+4*x[1]
    val3=4*x[0]+12*(x[1]**2)-26
    return np.array([[val1,val2],[val2,val3]])
def damped_newton(x0):
    x=[]
    x.append(x0)
    k=0
    
    tol=1e-12
    tol1=1
    while(tol1>tol):
        s=la.solve(hessian(x[k]),-df(x[k]))
        alpha=spo.line_search(f,df,x[k],-df(x[k]))
        c=alpha[0]
        if(alpha[0]==None):
            c=1
        
        xnew=x[k]+s*c
        x.append(xnew)
        tol1=la.norm(x[k+1]-x[k])
        a=x[k+1]
        k+=1
    return a,x
p=np.array([2,2])
q,r=damped_newton(p)
it_array=np.array(r)
xmesh, ymesh = np.mgrid[-5:5:50j,-5:5:50j]
fmesh = f(np.array([xmesh, ymesh]))
pt.axis("equal")
pt.contour(xmesh, ymesh, fmesh, 200)
pt.plot(it_array.T[0], it_array.T[1], "x-")

p2=np.array([2,-1])
q2,r2=damped_newton(p2)
it_array2=np.array(r2)
xmesh, ymesh = np.mgrid[-5:5:50j,-5:5:50j]
fmesh = f(np.array([xmesh, ymesh]))
pt.axis("equal")
pt.contour(xmesh, ymesh, fmesh)
pt.plot(it_array2.T[0], it_array2.T[1], "x-")

p3=np.array([-2,2])
q3,r3=damped_newton(p3)
it_array3=np.array(r3)
xmesh, ymesh = np.mgrid[-5:5:50j,-5:5:50j]
fmesh = f(np.array([xmesh, ymesh]))
pt.axis("equal")
pt.contour(xmesh, ymesh, fmesh)
pt.plot(it_array3.T[0], it_array3.T[1], "x-")

p4=np.array([-2,-2])
q4,r4=damped_newton(p4)
it_array4=np.array(r4)
xmesh, ymesh = np.mgrid[-5:5:50j,-5:5:50j]
fmesh = f(np.array([xmesh, ymesh]))
pt.axis("equal")
pt.contour(xmesh, ymesh, fmesh)
pt.plot(it_array4.T[0], it_array4.T[1], "x-")

pt.show()   