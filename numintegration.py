from __future__ import division
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as pt
import scipy.integrate as spi 

def f(x):
	return 4/(1+x**2)
def composite_midpoint(f,a,b,h):
	n=(b-a)/h + 1
	x=np.linspace(a,b,n)
	sum=0
	n=len(x)
	for i in range(0,n-1):
		sum=sum+h*f((x[i]+x[i+1])/2)
	
	return sum
def composite_trapezoidal(f,a,b,h):
	n=(b-a)/h + 1
	x=np.linspace(a,b,n)
	sum=0
	n=len(x)
	for i in range(0,n-1):
		sum=sum+h*(f(x[i])+f(x[i+1]))/2
	return sum
def composite_simpson(f,a,b,h):
	n=(b-a)/h + 1
	x=np.linspace(a,b,n)
	sum=0
	n=len(x)
	for i in range(0,n-1):
		sum=sum+h*(f(x[i])+f(x[i+1]) + 4*f((x[i+1]+x[i])/2))/6
	return sum

def monte_carlo(f,a,b,h):
	n=1/h+1
	sum=0
	x=np.linspace(a,b,n)
	n=len(x)
	for i in range(0,n):
		sum=sum+f(x[i])
	avg=sum/n
	return avg 

h=np.array([0.001,0.01,0.1,1])
ans_mp=[]
e_mp=[]
ans_trap=[]
e_trap=[]
ans_simp=[]
e_simp=[]
ans_mc=[]
e_mc=[]
n=len(h)
#print n
for i in range(0,n):
		ans1=composite_midpoint(f,0,1,h[i])
		e_mp.append(np.abs(ans1-np.pi))
		ans2=composite_trapezoidal(f,0,1,h[i])
		e_trap.append(np.abs(ans2-np.pi))
		print e_trap
		ans3=composite_simpson(f,0,1,h[i])
		e_simp.append(np.abs(ans3-np.pi))
		ans4=monte_carlo(f,0,1,h[i])
		e_mc.append(np.abs(ans4-np.pi))
		ans_mp.append(ans1)
		ans_trap.append(ans2)
		ans_simp.append(ans3)
		ans_mc.append(ans4)
		#print ans4

e_log_mp=[]
e_log_mp=np.log(e_mp)
pt.plot(np.log(h),e_log_mp,label='Midpoint')
pt.plot(np.log(h), np.log(e_trap), label='Trapezoidal')
pt.plot(np.log(h), np.log(e_simp),label='Simpsons')
pt.plot(np.log(h),np.log(e_mc),label='Monte Carlo')
pt.xlabel("log h")
pt.ylabel("log E")
pt.legend()
p_mp=np.log(e_mp[1]/e_mp[0])/np.log(h[1]/h[0])
print p_mp
p_trap=np.log(e_trap[1]/e_trap[0])/np.log(h[1]/h[0])
print p_trap
p_simp=np.log(e_simp[2]/e_simp[1])/np.log(h[2]/h[1])
print p_simp
p_mc=np.log(e_mc[1]/e_mc[0])/np.log(h[1]/h[0]) 
print p_mc
pt.show()
