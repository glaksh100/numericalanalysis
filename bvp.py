from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
import scipy.sparse as sps
import scipy.sparse.linalg as sla
import math
def bvp_solve(p,q,f,x,g,h1):
	n = len(x)
	h = x[1] - x[0]
	second_deriv=sps.diags([1,-2,1], offsets=np.array([-1,0,1])+1, shape=(n-2, n))/(h**2)
	factor_u= sps.diags([q(x[1:])], offsets=[1], shape=(n-2, n))
	first_deriv_term1=sps.diags([-p(x[1:n-1])], offsets=np.array([-1])+1, shape=(n-2, n))/(2*h)
	first_deriv_term2=sps.diags([p(x[1:n-1])], offsets=np.array([1])+1, shape=(n-2, n))/(2*h)
	first_deriv=first_deriv_term2+first_deriv_term1
	#print first_deriv.todense()
	A_int=factor_u+second_deriv+first_deriv
	A = sps.vstack([sps.coo_matrix(([1], ([0],[0])), shape=(1, n)),A_int, sps.coo_matrix(([1], ([0],[n-1])), shape=(1, n)),])
	A = sps.csr_matrix(A)
	
	B=A.todense()
	cond=la.cond(B)
	rhs=f(x)
	rhs[0]=g
	rhs[-1]=h1
	sol=sla.spsolve(A,rhs)

	return sol,cond

def u0(x):
	return 1/3*(np.exp(-4*x)) +2/3*(np.exp(2*x))
def p0(x):
	return 2
def q0(x):
	return -8
def f0(x):
	return x*0
def u(x):
	return np.sin(np.log(x))
def p(x):
	return 2/x
def q(x):
	return -2/(x**2)
def f(x):
	return (1/(x**2))*(np.cos(np.log(x))-np.sin(np.log(x)))

def f2(x):
	return 2*(np.sin(x)) + 4*np.cos(x)
def u2(x):
	return -np.sin(x) + 3*np.cos(x)
def p2(x):
	return -1
def q2(x):
	return -2
gk=u0(0)
h0=u0(1)

g=u(1)
h1=u(2)
g2=u2(0)
h2=u2(np.pi/2)
e0=[]
e=[]
e2=[]
emax=[]
n=[]
c=[]
i=0
for k in range(3,11):
	
	n.append(2**k)
	mesh0=np.linspace(0,1,n[i])
	mesh=np.linspace(1, 2, n[i])
	mesh2 = np.linspace(0, np.pi/2, n[i])
	sol0,cond0=bvp_solve(p0,q0,f0,mesh0,gk,h0)
	c.append(cond0)
	#sol,cond=bvp_solve(p,q,f,mesh,g,h1)
	sol2,cond2=bvp_solve(p2,q2,f2,mesh2,g2,h2)
	#e0.append(la.norm(u0(mesh)-sol0,ord=np.inf)/la.norm(u0(mesh)))
	#e.append(la.norm(u(mesh)-sol,ord=np.inf)/la.norm(u(mesh)))
	e2.append(la.norm(u2(mesh2)-sol2,ord=np.inf)/la.norm(u2(mesh)))
	i+=1
print c


#pt.plot(np.log(n),np.log(c))
#pt.figure(1)
pt.plot(np.log(n),np.log(e2),label='error')
#pt.plot(np.log(n), np.log(n),label='slope 2')
pt.xlabel("Log n")
pt.ylabel("Log e")
#pt.title("C plot ")
#pt.plot(mesh2,sol2,label='Approx sol')
#pt.plot(mesh2,u2(mesh2),label='True sol')
pt.title("iii- Error plot")
pt.legend()

pt.show()

#print A.todense()



