import math
import numpy as np
import matplotlib.pyplot as plt
######Constants#####
a = 0.05
t0 = 0.0
tend = 1.5
x0 = 0.0
xend= 1.0
dt = 0.001
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
    NUM[i,0] = 1.0
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
    NUM[i,0] = 1                 #NUM
    NUM[i,0] = 0                    #Validation 
###Starter Step####
for i in range(1,xsteps-1):
    NUM[i,1] =-a*dt*(NUM[i+1,0] + 2*NUM[i,0]+NUM[i-1,0])/2.0/dx + 0.5*math.sin(3*math.pi*X[i])*(NUM[i+1,0]-NUM[i-1,0])/2.0/dx +NUM[i,0]
############Matrix###############
A[0,0]=1.0
A[xsteps-1,xsteps-1]= 1.0
for i in range(1,xsteps - 1):

    #SUB
    if (i>0):
        A[i,i-1] = a*dt/dx/dx + 0.5*math.sin(3*math.pi*X[i])*dt/2.0/dx
    #Main
    A[i,i] = -2*a*dt/dx/dx-1.5
    
    #Super
    if(i<xsteps-1):
        A[i,i+1] = a*dt/dx/dx - 0.5*math.sin(3*math.pi*X[i])*dt/2.0/dx

##############
#Numerical Solving
#Filling RHS      
for n in range(1,tsteps-1):
    for i in range(1,xsteps-1):
        RHS[i]=0.5*NUM[i,n-1]-2.0*NUM[i,n] +  dt*(-9.0*a*math.pi**2.0*T[n]*math.sin(3*math.pi*X[i])-\
                                                   math.sin(3*math.pi*X[i])\
                                                   -0.5*math.sin(3*math.pi*X[i])*3*math.pi*T[n]*math.cos(3*math.pi*X[i]))
    NUM[:,range(n+1,n+2)]=np.linalg.solve(A,RHS)#SOLVE
#####Analytical Validation#####
ANA = np.zeros([xsteps,tsteps])
 
for i in range(0,xsteps):
    for n in range(0,tsteps):
        ANA[i,n] = T[n]*math.sin(3*math.pi*X[i])
###Plots
plt.figure()
#levels = np.linspace(0.0,1.0,50)
x,y = np.meshgrid(X,T,indexing='ij')
cp = plt.contourf(x,y,NUM,50)
plt.colorbar(cp)
plt.show()
cp = plt.contourf(x,y,ANA,50)
plt.colorbar(cp)
plt.show()
cp = plt.contourf(x,y,abs(NUM-ANA),50)
plt.colorbar(cp)
plt.show()
