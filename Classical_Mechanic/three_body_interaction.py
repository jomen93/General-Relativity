#three body interaction

import math
from math import *
from visual import *


# definition of the evolution of the system with the gravitaional interaccion

def fv1(t,x,v,m,n):             # this funtion find the force over each body 

    vec=vector(0,0,0)
    dxi=[vec]*n;
    G=-1
    d=10e-3                     # Auxiliar variable for avoid the divergence
    							# if the distance is minor of "d" between two bodies the function put zero in this case

    for i in range(0,n):
        for j in range(0,n):
            if j!=i:            
                rij=sqrt(dot(x[i]-x[j],x[i]-x[j]))
                if rij>d:       
                    dxi[i]=dxi[i]+G*m[j]*(x[i]-x[j])/pow(rij,3)                                            
        
    return dxi


def fx(v,dt,n):                    #Velocity of the bodies 
    return v
               
def operation(kix,kiv,xio,vio,n,var,h):         # Auxiliary function for RK4
    vec=vector(0,0,0)
    xip=[vec]*n;
    vip=[vec]*n;
    kipx=[vec]*n;
    kipv=[vec]*n;
    for i in range (0,n):
        kipx[i]=h*kix[i]
        kipv[i]=h*kiv[i]
        if var==1:
            xip[i]=xio[i]+kipx[i]/2
            vip[i]=vio[i]+kipv[i]/2
        else:
            xip[i]=xio[i]+kipx[i]
            vip[i]=vio[i]+kipv[i]

    return xip,vip,kipx,kipv
    
def operation2(k1,k2,k3,k4,xio2,n):				# Auxiliary function for RK4
    vec=vector(0,0,0)
    xip2=[vec]*n;
        
    for i in range(0,n):
        xip2[i]=xio2[i]+(k1[i]+2*(k2[i]+k3[i])+k4[i])/6.
    
    return xip2

#Metodo Ruge-kutta 
def rk4(fv,fx,xi,vi,tf,ti,m,n):          #xi is a vector with all initial positions 
        								
        var=1;
        h=tf-ti
           
        k1x=fx(vi,h,n)
        k1v=fv(ti,xi,vi,m,n)
         
        xiv,viv,k1x,k1v=operation(k1x,k1v,xi,vi,n,var,h)                       
        k2x=fx(viv,h,n)
        k2v=fv(ti+h/2.,xiv,viv,m,n)
        
         
        xiv,viv,k2x,k2v=operation(k2x,k2v,xi,vi,n,var,h)
        k3x=fx(viv,h,n)
        k3v=fv(ti+h/2.,xiv,viv,m,n)
         
        var=0
        xiv,viv,k3x,k3v=operation(k3x,k3v,xi,vi,n,var,h)
        k4x=fx(viv,h,n)
        k4v=fv(ti+h,xiv,viv,m,n)

        xiv,viv,k4x,k4v=operation(k4x,k4v,xi,vi,n,var,h)
        
        var=3
        xiv=operation2(k1x,k2x,k3x,k4x,xi,n)
        xf=xiv

        viv=operation2(k1v,k2v,k3v,k4v,vi,n)
        vf=viv
        #print "Metodo fin",xi,xf
        return xf,vf
        
## In this part we make the scene for the vpython visualisation

nc=3        		# Number of bodies
ni=1000.         	# Number of iterations
Ttot=1.				 
dt=Ttot/ni 			# step time
ti=0.				# initial time 
				
#make vectors for the initial conditions of the bodies 

mc=[0]*nc
xic=[0]*nc
xfc=[0]*nc
vic=[0]*nc
vfc=[0]*nc

#Definition of the mass for each body 

mc[0]=1
mc[1]=1
mc[2]=1


               
# Body 1
vic[0]=0.5*vector(0.93240737,0.86473146,0)
xic[0]=vector(0.97000436, -0.243087530,0)
xfc[0]=xic[0]
# Body 2
vic[1]=0.5*vector(0.93240737,0.86473146,0)
xic[1]=-xic[0]
xfc[1]=xic[1]
# Body 3
#vic[2]=vector(-0.93240737,-0.86473146,0)
vic[2]=vector(0,0,0)
xic[2]=vector(0,0,0)
xfc[2]=xic[1]


# set-up the characteristic of the elements in the visualisation 

ball1 = sphere(pos=(0,0,0),radius=2.5e-2, color=color.red)
ball1.pos = xic[0]
ball1.trail = curve(color=ball1.color)
ball2 = sphere(pos=(0,0,0),radius=2.5e-2, color=color.green)  
ball2.pos = xic[1]
ball2.trail = curve(color=ball2.color)
ball3 = sphere(pos=(0,0,0),radius=2.5e-2, color=color.blue)
ball3.pos = xic[2]
ball3.trail = curve(color=ball3.color)

n=0
n1=100


while True:

    #scene.autoscale = False
    n=n+1
    rate(500)
    
    tf=ti+dt

    #Body 1 data actualization
   
    xfc,vfc=rk4(fv1,fx,xic,vic,tf,ti,mc,nc)
    
    ball1.pos = xfc[0]
    ball1.trail.append(pos=ball1.pos,retain=2100)
   
    ball2.pos = xfc[1]
    ball2.trail.append(pos=ball2.pos,retain=10000) #,retain=10000

    ball3.pos = xfc[2]
    ball3.trail.append(pos=ball3.pos,retain=2100) #,retain=10000

    # "retain" is the amount of trajectory visualisation in the scene of vpython 
  
    
    
    vic=vfc
    xic=xfc
    t=tf
