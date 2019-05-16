import math
import numpy as np
import matplotlib.pyplot as plt
#Numerical#
###################
import math
import numpy as np
import matplotlib.pyplot as plt
######Constants#####
nu = 0.005
t0 = 0.0
tend = 0.6
x0 = 0.0
xend = 1.0
dt = 0.01
dx = 0.01
#####Arrays######
#T Array
tsteps = int((tend-t0)/dt)+1
T = [0] *tsteps
for n in range(0,tsteps-1):
    T[n+1] = T[n] +dt
#X Array
xsteps = int((xend-x0)/dx)+1
X = [0] * xsteps
for i in range(0,xsteps-1):
    X[i+1] = X[i] +dx
########################################
###Numerical#####
########################################
NUM = np.zeros([xsteps,tsteps]) #T(x,t) : T(i.n)
for i in range(0,xsteps):
    NUM[i,0] = math.sin(math.pi*6*X[i])
#RHS
RHS = np.zeros([xsteps,1])
RHS[0] = 0
RHS[xsteps-1]= 0
A=np.identity(xsteps)
#Initial Conditions
for n in range(0,tsteps):
    NUM[0,n] = 0
    NUM[1,n] = 0
for i in range (0,xsteps):
    NUM[i,0] = math.sin(math.pi*6*X[i])                #NUM
############Matrix###############
A[0,0]=1.0
A[xsteps-1,xsteps-1]= 1.0
for i in range(1,xsteps - 1):
  
    #SUB
    if (i>0):
        A[i,i-1] = -dt/dx/dx*nu
    #Main
    A[i,i] = 2*dt*nu/dx/dx+1
    
    #Super
    if(i<xsteps-1):
        A[i,i+1] = -dt/dx/dx*nu
##############
#Numerical Solving
#Filling RHS      
for n in range(0,tsteps-1):
    for i in range(0,xsteps-1):
        RHS[i]=NUM[i,n]- dt*NUM[i,n]/dx/2.0*(NUM[i+1,n]-NUM[i-1,n])+\
        dt*(-math.sin(6.0*math.pi*X[i])*math.exp(-T[n])+math.sin(6.0*math.pi*X[i])*math.exp(-T[n])*6.0*math.pi*math.cos(6.0*math.pi*X[i])*math.exp(-T[n])\
            +nu*36.0*(math.pi**2)*math.sin(6.0*math.pi*X[i])*math.exp(-T[n]))
    NUM[:,range(n+1,n+2)]=np.linalg.solve(A,RHS)#SOLVE
#####Analytical Validation#####
ANA = np.zeros([xsteps,tsteps])
  
for i in range(0,xsteps):
    for n in range(0,tsteps-1):
        ANA[i,n] = math.exp(-T[n])*math.sin(6.0*math.pi*X[i])
###Plots
#Plots
plt.figure()
#levels = np.linspace(0.0,1.0,50)
x,y = np.meshgrid(X,T,indexing='ij')
#cp = plt.contourf(x,y,NUM,50)
cp = plt.contourf(x,y,ANA,50)
#cp = plt.contourf(x,y,abs(NUM-ANA),50)
plt.colorbar(cp)
plt.show()


