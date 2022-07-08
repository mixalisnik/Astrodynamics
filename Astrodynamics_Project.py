import numpy as np
from matplotlib import pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import solve_ivp
import math as m

##FUNCTIONS
#Rotation Matrix 3 
def R3(Om):
    R3=np.array([[np.cos(Om),np.sin(Om),0],[-np.sin(Om),np.cos(Om),0],[0,0,1]])
    return R3
#Rotation Matrix 1
def R1(i):
    R1=np.array([[1,0,0],[0,np.cos(i),np.sin(i)],[0,-np.sin(i),np.cos(i)]])
    return R1
#Keplerian Orbital Elements to Cartesian Coordinates and Velocities in ECI reference frame
def kep2car(kep,GM):
    #This function takes input keplerian orbital elements kep=(a,e,i,RAAN,ω,f) as well as GM of orbit
    #and returns cartesian coordinates (x,y,z) and velocities (vx,vy,vz) in ECI reference frame
    a,e,i,Om,w,f=kep
    p=a*(1-e**2)
    r=p/(1+e*np.cos(f))
    r_perifocal_vec=np.array([[r*np.cos(f)],[r*np.sin(f)],[0]])
    v_perifocal_vec=np.array([[-np.sqrt(GM/p)*np.sin(f)],[np.sqrt(GM/p)*(e+np.cos(f))],[0]])
    Tpqw_eci=np.transpose(np.dot(np.dot(R3(w),R1(i)),R3(Om)))
    r_eci=np.dot(Tpqw_eci,r_perifocal_vec)
    x=float(r_eci[0])
    y=float(r_eci[1])
    z=float(r_eci[2])
    v_eci=np.dot(Tpqw_eci,v_perifocal_vec)
    vx=float(v_eci[0])
    vy=float(v_eci[1])
    vz=float(v_eci[2])
    return [x,y,z,vx,vy,vz]

#Cartesian Coordinates and Velocities in ECI reference frame to Keplerian Orbital Elements 
def car2kep(car,GM):
    #This function takes input cartesian coordinates (x,y,z) and velocities (vx,vy,vz) in ECI as well as GM of orbit
    #and returns keplerian orbital elements kep=(a,e,i,RAAN,ω,f).
    x,y,z,vx,vy,vz=car
    r_vec=[x,y,z]
    r=np.linalg.norm(r_vec)
    v_vec=[vx,vy,vz]
    v=np.linalg.norm(v_vec)
    Energy=0.5*v**2 - GM/r
    a=-GM/(2*Energy)
    rv=np.dot(r_vec,v_vec)
    e1=(v**2/GM - 1/r)*x - rv/GM*vx
    e2=(v**2/GM - 1/r)*y - rv/GM*vy
    e3=(v**2/GM - 1/r)*z - rv/GM*vz
    e_vec=[e1,e2,e3]
    e=np.linalg.norm(e_vec)
    h_vec=np.cross(r_vec,v_vec)
    h=np.linalg.norm(h_vec)
    i=np.arccos(h_vec[2]/h)
    K=[0,0,1]
    N_vec=np.cross(K,h_vec)
    N=np.linalg.norm(N_vec)
    if N_vec[1]>=0:
        RAAN=np.arccos(N_vec[0]/N)
    else:
        RAAN=2*np.pi - np.arccos(N_vec[0]/N)
    if e_vec[2]>=0:
        w=np.arccos(np.dot(N_vec,e_vec)/(N*e))
    else:
        w=2*np.pi - np.arccos(np.dot(N_vec,e_vec)/(N*e))
    vr=rv/r
    if vr>=0:
        f=np.arccos(np.dot(e_vec,r_vec)/(e*r))
    else:
        f=2*np.pi - np.arccos(np.dot(e_vec,r_vec)/(e*r))
    return [a,e,i,RAAN,w,f]
#Define System of Differential Equations of orbit
def rhs_2bp(t,X):
    
    x,y,z,vx,vy,vz=X
    mu=4902.8
    
    r=np.sqrt(x**2+y**2+z**2)
    r3=r*r*r
    x_dot=vx
    y_dot=vy
    z_dot=vz
    vx_dot=-mu*x/r3
    vy_dot=-mu*y/r3
    vz_dot=-mu*z/r3
    return [x_dot,y_dot,z_dot,vx_dot,vy_dot,vz_dot]
#Cartesian Orbit Propagator
def orbit_propagator(kep1):
    #This function takes keplerian orbital elements (a,e,i,RAAN,w,f) as input, and returns an array of coordinates (x,y,z) in ECI
    #that will be used to plot the full elliptical orbit in three dimensions. 
    f=np.linspace(0,2*np.pi,1000)
    x11=[]
    y11=[]
    z11=[]   
    rr=[]
   
    
    for i1 in f:
        
        kep1=[kep1[0],kep1[1],kep1[2],kep1[3],kep1[4],i1]
        a3=kep2car(kep1,GM)
        x11.append(a3[0])
        y11.append(a3[1])
        z11.append(a3[2])
        rr.append(np.sqrt(a3[0]**2 +a3[1]**2 +a3[2]**2))
        coords_arr=[x11,y11,z11,rr]
    return coords_arr

def plot3dconfig():
    #Initializes a figure to plot the orbit in three dimensions
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(-3*kep1[0],3*kep1[0])
    ax.set_ylim(-3*kep1[0],3*kep1[0])
    ax.set_zlim(-3*kep1[0],3*kep1[0])
    ax.set_xlabel('X (km)')
    ax.set_ylabel('Y (km)')
    ax.set_zlabel('Z (km)')
    return ax
#Three Dimensional Elliptical Orbit
def draw_elipse_config(kep1,ax):
    #This function plot the Elliptical orbit with keplerian orbital elements (a,e,i,RAAN,w,f). 
    coords_arr=orbit_propagator(kep1)
 
    plt.title("Three Dimensional Elliptical Orbit EX2TASK2A")
    plt.grid()
    ax.plot(coords_arr[0],coords_arr[1],coords_arr[2])
    plt.show() 
    return

def dist_from_surf(kep1):
    #Calculation and graphical represantation of distance of satellite from the surface of the moon
    rmoon=1738.1
    
    f=np.linspace(0,2*np.pi,1000)
    coords_arr=orbit_propagator(kep1)
    r=np.array(coords_arr[3])-rmoon
    plt.figure()
    plt.plot(f,r)
    plt.title("Distance from the Moon Surface - True Anomaly EX2TASK2A")
    plt.xlabel("True Anomaly f (radians)")
    plt.ylabel("r (km)")
    plt.grid()
    plt.show()
    return
##----------------------------------MAIN ----------------------------------------------------------------------
 ## Exercise 1
a=32000 #Semi major axis of satellite orbit
e=0.2 #Eccentricity of satellite orbit
b=a*np.sqrt(1-e**2) #Semi minor axis of satellite orbit
R_e=6378 #Earth's Radius
mu_earth=398600.433 #mu of Earth in km^3/s^2
x_sq=(1-R_e**2/b**2)*a**2 # Sattelite Enters Earth's Shadow at point (x0,y0)=(-x_sq^(1/2),R_e) and exits Earth's Shadow at point (x0,y0)=(-x_sq^(1/2),-R_e)
x=-np.sqrt(x_sq)
x_perifocal=-x+a*e #Transfer x from ECI to Perifocal reference frame by adding a*e.
r_a=np.sqrt(x_perifocal**2 + R_e**2) #Find the distance from the perifocal (0,0) point
p=a*(1-e**2)
cosfa=(p/r_a - 1)/e #Calculate true anomaly at point a: (Earth's Shadow Entry Point )
fa=np.arccos(cosfa)
fb=2*np.pi-fa #Calculate true anomaly at point b: (Earth's Shadow Exit Point)
Ea=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(fa/2)) #Calculate eccentric anomaly at point a: (Earth's Shadow Entry Point)
Eb=2*(np.pi - np.arctan(-np.sqrt((1-e)/(1+e))*np.tan(fb/2))) #Calculate eccentric anomaly at point b: (Earth's Shadow Exit Point)
Ma=Ea-e*np.sin(Ea) #Calculate mean anomaly at point a: (Earth's Shadow Entry Point)
Mb=Eb-e*np.sin(Eb) #Calculate mean anomaly at point b: (Earth's Shadow Exit Point)
T=np.sqrt(4*np.pi**2*a**3/mu_earth) #Calculate Period of the satellite orbit
n=2*np.pi/T #Calculate n
dt=(Mb-Ma)/n #Find dt between points A and B (Time Spent in Earth's Shadow.)

print("\n\n------ EXERCISE 1 ------\n")
print("Time Spent in Earth's Shadow : ",np.round(dt,5),"sec\n\n------------------------") #Note: dt rounded at 5th decimal point
 ## Exercise 2 
 #Task 1
m_moon=7.3477*10**22 #Moon mass in kg
mu_moon=4902.8 #km^3/s^2
m_earth=5.972*10**24 #Earth mass in kg
w_moon=2.66186*10**-6 #Rotational speed of moon in rad/s
T_moon=2*np.pi/w_moon #Orbital period of moon in sec
al=(T_moon**2*mu_moon/(4*np.pi**2))**(1/3) #Semi major axis of a lunar-sychronous orbiter
a_moon=348400 # Moon's semi major axis of its orbit around the Earth in km
e_moon=0.0549 # Moon's eccentricity of its orbit around the Earth
r_hill=a_moon*(1-e_moon) * (m_moon/(3*m_earth))**(1/3) #Moon's Hill's Radius in km




 #Task 2a   
#Initial Keplerian Elements (Angles in Rad)
kep1=[5737.4,0.61,np.deg2rad(57.82),0,np.deg2rad(90),0] #(a,e,i,RAAN,w,f)
rmoon=1738.1 #radius of the moon
GM=4902.8 #mu of the moon in km^3/sec^2

#Draw Figures
ax=plot3dconfig() #Initialize figure
draw_elipse_config(kep1,ax) #Plot Initial Eliptical orbit
dist_from_surf(kep1) #Plot satellite's distance of the moon R(f)

 ## Task 2b
#Initial Conditions as Keplerian Elements
mu=4902.8 #mu of the moon in km^3/sec^2
IC=kep1 


#Initial Conditions as Cartesian Coordinates abd Velocities in ECI
ICCAR=kep2car(IC,mu)
t0=0
tmax=30*86400 #Max Integration Time: 30 days

#Propagate orbit in Cartesian
sol=solve_ivp(rhs_2bp,[t0,tmax],ICCAR,t_eval=np.linspace(0,tmax,100000),method='DOP853',atol=1e-13,rtol=1e-13)
celem=sol.y[:,-1] #Cartesian Elements at t=30days

FC=car2kep(celem,mu) #Transform cartesian elements at t=30days to keplerian orbital elements. 
aend=FC[0]
eend=FC[1]
iend=FC[2]
RAANend=FC[3]
wend=FC[4]
fend=FC[5]
kepend=[aend,eend,iend,RAANend,wend,fend]
#Check for changes in keplerian elements after t=30days
stringmat=['a','e','i','RAAN','w','f']
print("\n\n------ EXERCISE 2 ------\n")
print("-------- TASK 1 --------\n")
print("-We find that the Moon's Hill's Radius is r_hill= ",np.round(r_hill,5), 
      "km and that the semi major axis of a lunar-sychronous orbiter is a= ",np.round(al,5),
      "km. \n-We know that the Moon's Hill's Radius is the distance from the moon within which smaller bodies",
      " (i.e. satellites) would tend to orbit around it, unaffected from pertubative forces like Earth's",
      " Gravitational Pull. \n-We observe that a>r_hill. This means, that there is no lunar-sychronous ",
      "stable orbit that can be achieved, becuase such an orbit would be outside Moon's Hill's Sphere",
      " and thus would be rendered unstable due to significant pertubations (from the Earth, Sun etc).")
print("\n\n------ EXERCISE 2 ------\n")
print("-------- TASK 2 --------\n")
print("Initial and final values (t=30days) \nof keplerian orbital elements\n------------------------")
for i in range(6):
    print("Initial ",stringmat[i],": ",IC[i],"\tFinal ",stringmat[i],": ",kepend[i])

print("\n-The distance of the sattelite from the moon's surface as a function of true anomaly f, R(f)is a bell-curve",
      " \nAs we expect:\n-R(f) is maximum at around f=π, which is consistent with the apogee's true anomaly. ",
      " \n-R(f) is minimum at approximately f=0 (f=2π) which is consistent with the perigee's true anomaly.\n")
print("-As we see in the table above, all keplerian orbital elements remain the same after an orbital propagation of 30days.",
      " There is a small yet noticeable accuracy error that has to do with the computational precision of the machine.") 
print("-We expect that the orbital elements will not change",
      " with time since there are no pertubative forces in the",
      "system we are studying.\n\n")
print("\n\n------ EXERCISE 2 ------\n")
print("-------- TASK 3 --------\n")
print("-The groundtrack of the satellite at t=30days is almost all surface area between φmin,φmax, for Λ=[-π,π]")
print("-The satellite's orbit is limited between φmin=-1->φmax=1 rad approximately.")


 #Task 2c 
#Find the groundtrack of the satellite  
w_moon=2.66186*10**-6 #moon's rotational speed in rad/s
t0=0
tmax=3*86400 #Max integration time 3 days
t=0
lamda_mat=[] #Array that stores values of λ 
phi_mat=[] #Array that stores values of φ



#Transform cartesian coordinates in ECI to coordinates in ECEF reference frame. Find λ and corresponding φ and fill lamda 
#and phi arrays.
for i in range(np.size(sol.y,axis=1)):
    r_ecef=np.dot(R3(w_moon*sol.t[i]),sol.y[:3,i])
    x_ecef=r_ecef[0]
    y_ecef=r_ecef[1]
    z_ecef=r_ecef[2]
    lamda=m.atan2(r_ecef[1],r_ecef[0])
    phi=phi=m.atan2(r_ecef[2],(np.sqrt(r_ecef[0]**2+r_ecef[1]**2)))
    lamda_mat.append(lamda)
    phi_mat.append(phi)

#Plot the groundtrack for 3 days, by plotting the coords until t=3days (~10000th element of arrays)
plt.figure()
plt.grid()
plt.title("Satellite Groundtrack for 3 days EX2TASK2C")
plt.scatter(lamda_mat[:10000],phi_mat[:10000],marker='.')
plt.xlabel("Λ (rad)")
plt.ylabel("φ (rad)")
#Plot the groundtrack for 30 days
plt.figure()
plt.grid()
plt.title("Satellite Groundtrack for 30 days EX2TASK2C")
plt.scatter(lamda_mat[:],phi_mat[:],marker='.')
plt.xlabel("Λ (rad)")
plt.ylabel("φ (rad)")


#Exercise 3


def changeOrbitalPlane(a_i,e_i,i_i,RAAN_i,w_i,i_f,RAAN_f):
    #This function changes the orbital plane of the orbit. It takes initial keplerian elements as input,
    #as well as the target inclination and RAAN.
    #It returns the angle theta1 which is the true anomaly when the orbital plane change must happen.
    #Additional outputs Dv cost and w_f (target argument of perigee)
    mu=398600.433 #mu of earth in km^3/s^2
    DRAAN=RAAN_f-RAAN_i
    cosa=np.cos(i_i)*np.cos(i_f)+np.sin(i_i)*np.sin(i_f)*np.cos(DRAAN)
    a=m.acos(cosa)
    cosu2=(np.cos(i_i)-cosa*np.cos(i_f))/(-np.sin(a)*np.sin(i_f))
    cosu1=(np.cos(i_f)-cosa*np.cos(i_i))/(-np.sin(a)*np.sin(i_i)) 
    sinu2=np.sin(i_i)*np.sin(DRAAN)/np.sin(a)
    sinu1=np.sin(i_f)*np.sin(DRAAN)/np.sin(a)
    u1=m.atan2(sinu1,cosu1)
    u2=m.atan2(sinu2,cosu2)
    theta1=u1-w_i 
    w_f=u2-theta1
    p=a_i*(1-e_i**2)
    v_theta=np.sqrt(mu/p)*(1+e_i*np.cos(theta1))
    Dv=2*v_theta*np.sin((i_f-i_i)/2)
    kep=[a_i,e_i,i_f,RAAN_f,w_f,theta1]
    elem=[theta1,Dv,w_f]
    return kep,elem

def changePeriapsisArg(a_i,e_i,w_2,theta2,w_f):
    #This function changes the periapsis argument of the orbit. It takes initial keplerian elements as input,
    #as well as the target w and theta2 angle when the change happens.
    
    #Outputs: Dv cost and w_3 (target argument of perigee), theta2a and theta3a
    w_3=w_f
    Dw=w_f-w_2
    theta_2_a=Dw/2
    theta_3_a=2*np.pi-Dw/2
    p=a_i*(1-e_i**2)
    Dv2=2*np.sqrt(mu/p)*e_i*np.sin(Dw/2)
    return [Dv2,w_3,theta_2_a,theta_3_a]

def timeOfFlight(a,e,theta1,theta2):
    #This function takes keplerian elements a and e as input, as well as the theta1 and theta2 true anomaly values.
    #It calculates the time the orbit needs to get from fa=theta1 to true fb=theta2
    if np.sqrt((1-e)/(1+e))*np.tan(theta1/2)>0:
        Ea=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(theta1/2))
    else:
        Ea=2*(np.pi - np.arctan(-np.sqrt((1-e)/(1+e))*np.tan(theta1/2)))
    if np.sqrt((1-e)/(1+e))*np.tan(theta2/2)>0:
        Eb=2*np.arctan(np.sqrt((1-e)/(1+e))*np.tan(theta2/2))
    else:
        Eb=2*(np.pi - np.arctan(-np.sqrt((1-e)/(1+e))*np.tan(theta2/2)))
    Ma=Ea-e*np.sin(Ea)
    Mb=Eb-e*np.sin(Eb)
    
    T=np.sqrt(4*np.pi**2*a**3/mu_earth)
    n=2*np.pi/T
    if theta2>theta1:
        dt=(Mb-Ma)/n
    else:
        dt=(Mb-Ma)/n + T
    return dt

def hohmann_dv(a_i,e_i,a_t,e_t):
    mu=398600.433 #mu of earth in km^3/s^2
    uai=np.sqrt(mu/a_i*(1-e_i)/(1+e_i))  # Velocity of apogee in initial orbit 
    up1=np.sqrt(mu/a_i*(1+e_i)/(1-e_i))  # Velocity of perigee in initial orbit
    ua2=np.sqrt(mu/a_t*(1-e_t)/(1+e_t))  # Velocity of apogee in target orbit
    up2=np.sqrt(mu/a_t*(1+e_t)/(1-e_t))  # Velocity of perigee in target orbit
    a_2=(a_i*(1-e_i) + a_t*(1+e_t))/2    # Semi major axis of transfer orbit
    e_2=(a_t*(1+e_t) - a_i*(1-e_i))/(a_t*(1+e_t)+a_i*(1-e_i)) #Eccentricity of transfer orbit
    uat=np.sqrt(mu/a_2*(1-e_2)/(1+e_2))
    upt=np.sqrt(mu/a_2*(1+e_2)/(1-e_2))
    du1=upt-up1
    du2=ua2-uat
    dv=du1+du2
    return dv
    
mu_earth=398600.433
#Initial Conditions in Cartesian
car_man=[-9141.878,-1648.0758,4141.679,-1.153,-5.31,-2.898]
kep_0=car2kep(car_man,mu_earth) #Initial Keplerian Elements
kep_target=[12940,0.2173,0.8692,1.448,2.721,2.827]

#Manuever 1 Plane Change
kep11,elem=changeOrbitalPlane(kep_0[0],kep_0[1],kep_0[2],kep_0[3],kep_0[4],kep_target[2],kep_target[3])
ax=plot3dconfig()
draw_elipse_config(kep_0,ax) #Draw Initial Orbit

draw_elipse_config(kep11,ax) #Draw Second Orbit after orbital plane change.

kep_01=kep_0.copy()
kep_01[-1]=+kep11[-1] #Change true anomaly to theta1

#Manuever 2 Change Periapsis Argument

elems=changePeriapsisArg(kep11[0],kep11[1],kep11[4],kep11[5],kep_target[4])

kep11[-1]=+elems[-2] #Change the true anomaly to theta2a

kep11[-2]=elems[1] #Change argument of perigee to w_3
draw_elipse_config(kep11,ax) #Draw Third Orbit

draw_elipse_config(kep_target,ax) #Draw Target Elipse
plt.title("Orbits during Manuevers EX3")
#Calculate time of flight

#Time Of Flight
dt1=timeOfFlight(kep_0[0],kep_0[1],kep_0[-1],kep11[-1]) #dt1 : time of flight form t0 until we reach
                                                        # the plane change manuever point
dt2=timeOfFlight(kep_0[0],kep_0[1],kep11[-1],elems[-2]) #dt2 : time of flight from the change of orbital plane 
                                                        #until we reach the apsis rotation manuever point
dt3=timeOfFlight(kep_0[0],kep_0[1],elems[-1],2*np.pi)   #d3 : time of flight from the apsis rotation until we reach
                                                        #the Hohmann transfer manuever point
dt4=timeOfFlight(kep_target[0],kep_target[1],np.pi,kep_target[-1]) #dt4: Wait till we reach final point
DTall=dt1+dt2+dt3+dt4 
print("\n\n------ EXERCISE 3 ------\n\nThe time needed to reach the target orbit is :",DTall,"seconds\n")
print("Total ΔV cost is : ",elem[1]+elems[0]+hohmann_dv(kep_0[0],kep_0[1],kep_target[0],kep_target[1]),"km/s\n------------------------")

    

