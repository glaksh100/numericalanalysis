from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
import scipy.sparse as sps
import scipy.sparse.linalg as sla
import cmath
from mpl_toolkits.mplot3d import Axes3D

n=400
delx=1/n
u=np.eye(n)
C=(np.roll(u, -1, axis=-1) - np.roll(u, 1, axis=-1))/(2*delx)
dt1=10/4000
eig_val_C,eig_vec_C=la.eig(C)
scaled_eig_val_C=dt1*eig_val_C
max_eig_val_C=scaled_eig_val_C.max()
min_eig_val=scaled_eig_val_C.min()
print max_eig_val_C, min_eig_val
#pt.plot(scaled_eig_val_C.real,scaled_eig_val_C.imag,'bo')
#pt.title("The point lies almost on the boundary of the stability region for Forward Euler")


#pt.show()
du = (np.roll(u, -1, axis=-1) - np.roll(u, 0, axis=-1))/delx
dt= 1/400
eig_val_u,eig_vec_u=la.eig(du)
scaled_eig_val_u=dt*eig_val_u
#pt.plot(scaled_eig_val_u.real,scaled_eig_val_u.imag,'bo')
#pt.title("The point lies almost on the boundary of the stability region for Forward Euler")
#pt.plot(dt,0,'go')
#print scaled_eig_val_u.min() , scaled_eig_val_u.max()
#pt.show()

def u(t,y):
	return np.sin(2*np.pi*t)
def w(t,y):
	return np.abs((4*t)%2 -1)
def z(t,y):
	if(t>0 and t<0.5):
		return 1
	else:
		return 0
def fw_euler_step(y, t, h, f):
    return y + h * f(t, y)
def rk4_step(y, t, h, f):
    k1 = f(t, y)
    k2 = f(t+h/2, y + h/2*k1)
    k3 = f(t+h/2, y + h/2*k2)
    k4 = f(t+h, y + h*k3)
    return y + h/6*(k1 + 2*k2 + 2*k3 + k4)

def ce_rk(delt,tf,g,u0):
	nt=(tf-0)/delt
	tspan=np.linspace(0,tf,nt)
	u_sol=[]
	u_sol.append(u0)
	nt=len(tspan)
	for i in range(0,nt-1):
		val=rk4_step(u_sol[i],tspan[i],delt,g)
		u_sol.append(val)
	return u_sol
def ce_euler(delt,tf,g,u0):
	nt=(tf-0)/delt
	tspan=np.linspace(0,tf,nt)
	u_sol=[]
	u_sol.append(u0)
	nt=len(tspan)
	for i in range(0,nt-1):
		val=fw_euler_step(u_sol[i],tspan[i],delt,g)
		u_sol.append(val)
	return u_sol
def upwind_rk(delt,tf,g,u0):
	nt=(tf-0)/delt
	tspan=np.linspace(0,tf,nt)
	u_sol=[]
	u_sol.append(u0)
	nt=len(tspan)
	for i in range(0,nt-1):
		val=rk4_step(u_sol[i],tspan[i],delt,g)
		u_sol.append(val)
	return u_sol

def upwind_euler(delt,tf,g,u0):
	nt=(tf-0)/delt
	tspan=np.linspace(0,tf,nt)
	u_sol=[]
	u_sol.append(u0)
	nt=len(tspan)
	for i in range(0,nt-1):
		val=fw_euler_step(u_sol[i],tspan[i],delt,g)
		u_sol.append(val)
	return u_sol
a=rk4_step(0,1,0.1,u)
nts=4000
ts=np.linspace(0,10,nts)
b=upwind_rk(1/400,10,z,1)
mesh=np.linspace(0,1,nts)
print len(b)
pt.plot(ts,b)
#pt.title("Upwind with RK - i")
pt.subplot(111, projection="3d")
pt.gca().plot_wireframe(mesh, ts, b)
pt.show()






