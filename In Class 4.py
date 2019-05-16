import math
import numpy as np
import matplotlib.pyplot as plt




#Numerical#
###################
import math
import numpy as np
import matplotlib.pyplot as plt



######Constants#####


a = 0.009
b = 25000.0
h = 10.0
p = 2750.0
cp = 800.0
kl = 0.2
r0 = 0.0
rend = a

t0 = 0.0
tend = 2400.0

dt = 0.1
dr = 0.0001
#####Arrays######
#T Array
tsteps = int((tend-t0)/dt)+1
T = [0] *tsteps
for n in range(0,tsteps-1):
    T[n+1] = T[n] +dt


#r Array

rsteps = int((rend-r0)/dr)+1
R = [0] * rsteps
for i in range(0,rsteps-1):
    R[i+1] = R[i] +dr


########################################
###Numerical#####
########################################

NUM = np.zeros([rsteps,tsteps]) #T(r,t) : T(i.n)

for i in range(0,rsteps):
    NUM[i,0] = 0




#RHS
RHS = np.zeros([rsteps,1])


A=np.identity(rsteps)


#Initial Conditions
#for n in range(0,tsteps):
#    NUM[0,n] = 0
#    NUM[1,n] = 0
for i in range (0,rsteps):
    NUM[i,0] = 0                #NUM
                        #Validation 





############Matrix###############


####K####

al = 0.2/p/cp



A[0,0]= 1.0
A[0,1] = -1.0


A[rsteps-1,rsteps-1] = kl+dr*h       
A[rsteps-1,rsteps-2]= -kl


for i in range(1,rsteps - 1):
   
    #SUB
    if (i>0):
        A[i,i-1] = dt*al/(R[i]*2.0*dr) - dt*al/(dr**2.0)
    #Main
    A[i,i] = 1.0+2.0*dt*al/(dr**2.0)
    
    #Super
    if(i<rsteps-1):
        A[i,i+1] = -dt*al/(R[i]*2.0*dr)-dt*al/(dr**2.0)

##############
#Numerical Solving
#Filling RHS      
for n in range(1,tsteps-1):
    RHS[0] = 0
    RHS[rsteps-1]=-0.1*math.exp(-NUM[i,n])*NUM[i,n] + 0.1*math.exp(-NUM[i,n])*NUM[i-1,n]\
                   
    for i in range(1,rsteps-1):
        anl = (1.0/(p*cp))*-0.1*math.exp(-NUM[i,n])
        RHS[i]= NUM[i,n] - dt*anl *(NUM[i+1,n]**2.0+2.0*NUM[i+1,n]*NUM[i-1,n]+NUM[i-1,n]**2.0)/(4.0*dr**2.0)\
                +dt*anl/R[i]*((NUM[i+1,n]+NUM[i-1,n])/(2.0*dr))+dt*anl*((NUM[i+1,n]-2.0*NUM[i,n]+NUM[i-1,n])/(dr**2.0))\
                +dt*b/p/cp+dt*(math.cos(200.0*R[i])*math.cos(T[n])+0.1/p/cp*math.exp(math.cos(200.0*R[i])*math.sin(T[n]))*(-200.0*math.sin(200.0*R[i])*math.cos(T[n]))**2.0-(anl+al)/R[i]*(-200.0*math.sin(200.0*R[i])*math.cos(T[n]))-(anl+al)*(-200.0**2.0*math.cos(200.0*R[i])*math.cos(T[n]))-b/p/cp)
            
    NUM[:,range(n+1,n+2)]=np.linalg.solve(A,RHS)#SOLVE


#####Analytical Validation#####
ANA = np.zeros([rsteps,tsteps])
    
for i in range(0,rsteps):
    for n in range(0,tsteps):
        ANA[i,n] = math.cos(200.0*R[i])*math.sin(T[n])


###Plots

plt.plot(R,NUM[:,tsteps-1], linewidth = 3, color='blue')

plt.show()

#Plots
plt.figure()
#levels = np.linspace(0.0,1.0,50)
x,y = np.meshgrid(R,T,indexing='ij')
#cp = plt.contourf(x,y,abs(B_EX-U),levels)
#cp = plt.contourf(x,y,B_EX,levels)
cp = plt.contourf(x,y,NUM,50)
plt.colorbar(cp)
plt.show()
cp = plt.contourf(x,y,ANA,50)
plt.colorbar(cp)
plt.show()
