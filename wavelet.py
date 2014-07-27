from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
import scipy.optimize as spo
def forward_euler(f,h,t0=0,tf=1):
    y=[]
    n=(tf-t0)/h 
    t=np.linspace(0,1,n)
    #print t
    n=len(t)
    y.append(1)
    for k in range(0,n):
        val=y[k]+h*f(t[k],y[k])
        y.append(val)
    return y
def backward_euler(f,h,t0=0,tf=1):
    y=[]
    n=(tf-t0)/h + 1
    t=np.linspace(0,1,n)
    #print t
    n=len(t)
    y.append(1)
    for k in range(0,n-1):
        val=spo.newton(g,y[k],args=(y[k],t[k+1],h))
        y.append(val)
    return y
def rk_4(f,h,t0=0,tf=1):
    y=[]
    n=(tf-t0)/h 
    t=np.linspace(0,1,n)
    n=len(t)
    y.append(1)
    for k in range(0,n):
        k1=f(t[k],y[k])
        k2=f((t[k]+h/2),(y[k]+h*k1/2))
        k3=f((t[k]+h/2),(y[k]+h*k2/2))
        k4=f((t[k]+h),(y[k]+h*k3))
        val=y[k]+h*(k1+2*k2+2*k3+k4)/6
        y.append(val)
    return y
def f(x,y):
    return (-200*x*(y**2))
def g(x,x0,t,h):
    return (200*t*h*(x**2) + x - x0)
def true_y(t):
    return 1/(1+100*(t**2))

h=np.array([0.125,0.25,0.5,1])
ans_fe=[]
ans_be=[]
ans_rk=[]
e_fe=[]
e_be=[]
e_rk=[]
t_1=np.linspace(0,1,9)
t_2=np.linspace(0,1,5)
t_3=np.linspace(0,1,3)
t_4=np.linspace(0,1,2)
sum=0
sum=0
n=len(h)
for i in range(0,n):
    val1=forward_euler(f,h[i])
    val2=backward_euler(f,h[i])
    val3=rk_4(f,h[i])
    ans_fe.append(val1)
    ans_be.append(val2)
    ans_rk.append(val3)
    l_1=len(ans_fe[i])
    l_2=len(ans_be[i])
    l_3=len(ans_rk[i])
    e_1=np.abs(ans_fe[i][l_1-1] - true_y(1))
    e_2=np.abs(ans_be[i][l_2-1] - true_y(1))
    e_3=np.abs(ans_rk[i][l_3-1] - true_y(1))
    e_fe.append(e_1)
    e_be.append(e_2)
    e_rk.append(e_3)
    #print ans_be
pt.figure(1)
pt.plot(ans_be[0],t_1,label='h=0.125') 
pt.plot(ans_be[1],t_2,label='h=0.25')
pt.plot(ans_be[2],t_3,label='h=0.5')
pt.plot(ans_be[3],t_4,label='h=1') 
pt.xlabel("t")
pt.ylabel("y")
pt.title("backward_euler") 
print ans_rk
pt.legend()
pt.show() 

pt.figure(2)
pt.plot(ans_rk[0],t_1,label='h=0.125') 
pt.plot(ans_rk[1],t_2,label='h=0.25')
pt.plot(ans_rk[2],t_3,label='h=0.5')
pt.plot(ans_rk[3],t_4,label='h=1') 
pt.xlabel("t")
pt.ylabel("y")
pt.title("rk_4") 
pt.legend()
pt.show() 

pt.figure(3)
pt.plot(ans_fe[0],t_1,label='h=0.125') 
pt.plot(ans_fe[1],t_2,label='h=0.25')
pt.plot(ans_fe[2],t_3,label='h=0.5')
pt.plot(ans_fe[3],t_4,label='h=1') 
pt.xlabel("t")
pt.ylabel("y")
pt.title("forward_euler") 
pt.legend()
pt.show() 

pt.figure(4)
pt.plot(np.log(h),np.log(e_fe), label ='forward_euler')
pt.plot(np.log(h),np.log(e_be), label ='backward_euler')
pt.plot(np.log(h),np.log(e_rk), label ='rk_4')
pt.xlabel("log h")
pt.ylabel("log E")
pt.legend()
pt.show()

#print ans_be[0][len(ans_be[0])-1)]

