# -*- coding: cp1252 -*-
# Solar system visualisation 
import math
from math import *
from visual import *

# Evolution rule, equations of motion

def fv1(t,x,v,m,n):             

    vec=vector(0,0,0)
    dxi=[vec]*n;
    G=-6.673884e-11
    d=5e7                       

    for i in range(0,n):
        for j in range(0,n):
            if j!=i:            
                rij=sqrt(dot(x[i]-x[j],x[i]-x[j]))
                if rij>d:
                    dxi[i]=dxi[i]+G*m[j]*(x[i]-x[j])/pow(rij,3)      
                else :
                    dxi[i]=dxi[i]
        
    return dxi

# Position of the system 

def fx(v):                   
    return v



def operation(kix,kiv,xio,vio,n,var,h):               

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

#Auxiliary function to find the new data in position and velocity with RK

def operation2(k1,k2,k3,k4,xio2,n):
    vec=vector(0,0,0)
    xip2=[vec]*n;
        
    for i in range(0,n):
        xip2[i]=xio2[i]+(k1[i]+2*(k2[i]+k3[i])+k4[i])/6.
    
    return xip2



# METODO RUGE-KUTTA

def rk4(fv,fx,xi,vi,tf,ti,m,n):         
        
        var=1;                          

        h=tf-ti
           
        k1x=fx(vi)
        k1v=fv(ti,xi,vi,m,n)
         
        xiv,viv,k1x,k1v=operation(k1x,k1v,xi,vi,n,var,h)                        
        k2x=fx(viv)
        k2v=fv(ti+h/2.,xiv,viv,m,n)
        
         
        xiv,viv,k2x,k2v=operation(k2x,k2v,xi,vi,n,var,h)
        k3x=fx(viv)
        k3v=fv(ti+h/2.,xiv,viv,m,n)
         
        var=0
        xiv,viv,k3x,k3v=operation(k3x,k3v,xi,vi,n,var,h)
        k4x=fx(viv,)
        k4v=fv(ti+h,xiv,viv,m,n)

        xiv,viv,k4x,k4v=operation(k4x,k4v,xi,vi,n,var,h)
        
        
        xiv=operation2(k1x,k2x,k3x,k4x,xi,n)
        xf=xiv

        viv=operation2(k1v,k2v,k3v,k4v,vi,n)
        vf=viv
        
        return xf,vf
        



#SECCION 2-1:
# Variables definitions for the system    
    
nc=11                                           #Number of bodies 
ni=1.                                           #Number of iteratiosn
Ttot=36000.
dt=Ttot/ni                                      #Step tiem
ti=0.
es=30.                                          #scale for can be see the moon in the scene of vpython 
ca=pi/180.                                      #degrees to radians

# vectors for mass, intial position, final position, intial velocity and final velocity, respectively 

mc=[0]*nc               
xic=[0]*nc              
xfc=[0]*nc              
vic=[0]*nc              
vfc=[0]*nc              


# Setp-up intital conditions of the system 

#Mass

mc[0]=1.989e30                  #Sun
mc[1]=3.302e23                  #Mercury
mc[2]=4.869e24                  #Venus    
mc[3]=5.9736e24                 #Earth
mc[4]=6.4185e24                 #Mars
mc[5]=7.349e22                  #Moon
mc[6]=1.899e27                  #Jupiter
mc[7]=5.688e26                  #Saturn
mc[8]=8.686e25                  #Uranus 
mc[9]=1.024e26                  #Neptun
mc[10]=1.25e22                  #Pluto


       
#Sun
vic[0]=vector(0.,0.,0.)
xic[0]=vector(0.,0.,0.)
xfc[0]=xic[0]

#Mercury
vp=59.25e3                 #perihelion velocity
pp=46.00e9                 #perihelion position
a=29.12478*ca              
i=7.00487*ca
vic[1]=vector(-vp*sin(a),vp*cos(a),0)
xic[1]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[1]=xic[1]

#Venus
vp=35.26e3                  #perihelion velocity
pp=107.48e9                 #perihelion postion 
a=54.8522*ca
i=3.39*ca
vic[2]=vector(-vp*sin(a),vp*cos(a),0)
xic[2]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[2]=xic[2]

#Earth

vp=30.288e3                  #perihelion velocity
pp=147.10e9                #perihelion postion 
a=282.94*ca
i=0.00*ca
vic[3]=vector(-vp*sin(a),vp*cos(a),0)
xic[3]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[3]=xic[3]

#Mars
vp=26.614e3                  #perihelion velocity
pp=206.67e9                  #perihelion postion 
a=286.46*ca
i=1.85*ca
vic[4]=vector(-vp*sin(a),vp*cos(a),0)
xic[4]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[4]=xic[4]

#Moon

vp=30.288e3-1000           #perihelion velocity
pp=147.10e9+3.84e8         #perihelion postion 
a=282.94*ca
i=0.00*ca
vic[5]=vector(-vp*sin(a),vp*cos(a),0)
xic[5]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[5]=xic[5]

#Jupiter
vp=13.73e3                  #perihelion velocity
pp=740.57e9                 #perihelion postion 
a=274.19*ca
i=1.3*ca
vic[6]=vector(-vp*sin(a),vp*cos(a),0)
xic[6]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[6]=xic[4]

#Saturn
vp=10.195e3                  #perihelion velocity
pp=1353.6e9                  #perihelion postion 
a=339.39*ca
i=2.48*ca
vic[7]=vector(-vp*sin(a),vp*cos(a),0)
xic[7]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[7]=xic[7]

#Uranus
vp=7.135e3                  #perihelion velocity
pp=2748.9e9                 #perihelion postion 
a=96.73*ca
i=0.769*ca
vic[8]=vector(-vp*sin(a),vp*cos(a),0)
xic[8]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[8]=xic[8]

#Neptun
vp=5.48e3                 #perihelion velocity
pp=4452.9e9                 #perihelion postion 
a=272.449*ca
i=1.769*ca
vic[9]=vector(-vp*sin(a),vp*cos(a),0)
xic[9]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[9]=xic[9]

#Pluto
vp=6.2618e3                  #perihelion velocity
pp=4438.6e9                #perihelion postion 
a=113.76*ca
i=17.15*ca
vic[10]=vector(-vp*sin(a),vp*cos(a),0)
xic[10]=vector(pp*cos(a)*cos(i),pp*sin(a)*cos(i),pp*sin(i))
xfc[10]=xic[10]

# Set-up scene of vpython visualisation, definitions of bodies
# For each object in the visualisation correspond one vector for the solution in the RK

ball1 = sphere(pos=xic[0],radius=5e9, color=color.yellow)
ball1.trail = curve(color=ball1.color)

ball2 = sphere(pos=xic[1],radius=3e9, color=(0.5,0.5,0.5))
ball2.trail = curve(color=ball2.color)

ball3 = sphere(pos=xic[2],radius=3e9, color=color.orange)
ball3.trail = curve(color=ball3.color)

ball4 = sphere(pos=xic[3],radius=3e9, color=color.blue)
ball4.trail = curve(color=ball4.color)

ball5 = sphere(pos=xic[4],radius=3e9, color=(0.5,0.2,0.))
ball5.trail = curve(color=ball5.color)

ball6 = sphere(pos=(xic[5]-xic[3])*es+xic[3],radius=1e9, color=(0.3,0.3,0.3))        #Scale for the moon visualisation
ball6.trail = curve(color=ball6.color)

ball7 = sphere(pos=xic[6],radius=4e9, color=(0.058,0.728,0.488))
ball7.trail = curve(color=ball7.color)

ball8 = sphere(pos=xic[7],radius=4e9, color=(0.8,0.4,0.4))
ball8.trail = curve(color=ball8.color)

ball9 = sphere(pos=xic[7],radius=4e9, color=(1,1,1))
ball9.trail = curve(color=ball9.color)

ball10 = sphere(pos=xic[9],radius=4e9, color=(0.,0.3,1.))
ball10.trail = curve(color=ball10.color)

ball11 = sphere(pos=xic[10],radius=4e9, color=(0.8,0.6,0.0))
ball11.trail = curve(color=ball11.color)


#Evolution

n=0

while True:                         

    
    #n=n+1                           
    rate(1e8)       # This important to visualisation !
    
    tf=ti+dt
    
    xfc,vfc=rk4(fv1,fx,xic,vic,tf,ti,mc,nc)         
    vic=vfc                                         
    xic=xfc
    ti=tf           
   
    # New data for position
    
    ball1.pos = xfc[0]
    ball1.trail.append(pos=ball1.pos)

    ball2.pos = xfc[1]
    ball2.trail.append(pos=ball2.pos,retain=1000)

    ball3.pos = xfc[2]
    ball3.trail.append(pos=ball3.pos,retain=1000)

    ball4.pos = xfc[3]
    ball4.trail.append(pos=ball4.pos,retain=1000)

    ball5.pos = xfc[4]
    ball5.trail.append(pos=ball5.pos,retain=5000)

    pl=(xfc[5]-xfc[3])*es+xfc[3]
    ball6.pos = pl
    ball6.trail.append(pos=ball6.pos,retain=400) #,retain=10000

    ball7.pos = xfc[6]
    ball7.trail.append(pos=ball7.pos)

    ball8.pos = xfc[7]
    ball8.trail.append(pos=ball8.pos)

    ball9.pos = xfc[8]
    ball9.trail.append(pos=ball9.pos)

    ball10.pos = xfc[9]
    ball10.trail.append(pos=ball10.pos)

    ball11.pos = xfc[10]
    ball11.trail.append(pos=ball11.pos)


    
