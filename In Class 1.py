
import math
import numpy as np
import matplotlib.pyplot as plt


#Constants
dx =0.001
dt =0.001


#Arrays

#T
t0 = 0.0
tend = 1.0
tsteps = int((tend-t0)/dt) +1
T = [0] * tsteps
for n in range(0,tsteps-1):
    T[n+1] = T[n] +dt

#X
x0 = 0.0
xend =2* math.pi
xsteps = int((xend-x0)/dx)+1
X = [0] * xsteps

for i in range(0,xsteps-1):
    X[i+1] = X[i] + dx

###Numerical Solutions###
NUM = np.zeros([tsteps,xsteps])
for n in range(0,tsteps):
    NUM[n,0] = 0
    
for i in range (0,xsteps):
    NUM[0,i]= (math.sin(3*X[i]))**2






for i in range(1,xsteps-1):
    
    for n in range(0,tsteps-1):
    

        NUM[n+1,i]= -NUM[n,i] * dt*((NUM[n,i]-NUM[n,i-1])/dx)+NUM[n,i] + dt* (-math.sin(T[n])*((math.sin(3.0*X[i]))**2.0) + \
                    6.0*((math.cos(T[n]))**2.0)*((math.sin(3.0*X[i]))**3.0)*math.cos(3.0*X[i]))



###Validation###
        
ANA = np.zeros([1,xsteps])

for i in range(0,xsteps):
    ANA[0,i] = math.cos(tend)*(math.sin(3*X[i])**2.0)


#plt.plot(X,NUM[tsteps-1,:], linewidth = 3, color = 'blue')
#plt.plot(X,ANA[0,:], linewidth = 3, color = 'green')

plt.plot(X,abs(NUM[tsteps-1,:] - ANA[0,:]), linewidth = 3, color = 'pink')
plt.show()



    
