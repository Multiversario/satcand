from scipy.integrate import solve_ivp,odeint
from scipy.optimize import minimize
import numpy as np
import sys

run_sim = 'T' #sys.argv[1]

#constants
G = 4.*np.pi**2
M_sun = 1.989e33 #Mass of Sun in g
AU = 149597870700. #AU in m
AU_km = AU/1000.
AU_cm = AU*100.

M_star = 1. #mass of host star (in M_sun)
M_E = 3.0035e-6 #M_sun
R_E = 4.2635e-5 #AU
rho_E = 5.515 # g/cm^3
rho_N = 1.64 # g/cm^3
rho_J = 1. # g/cm^3
rho_moon = 3.34 #g/cm^3

#Earth-like tidal parameters
alpha = 0.33  #moment of inertia
k2p = 0.299 #planet Love number
Q_p = 12 #tidal Quality factor

R_p = R_E  #Radius of planet (AU)
M_p = M_E  #Mass of planet 
a_p = 1# semimajor axis of planet (AU)
e_p = 0# eccentricity of planet
M_sat = 0.0123*M_E #Mass of Moon (M_sun)
R_H = a_p*((M_p+M_sat)/(3*M_star))**(1./3.) #Hill radius
mu = (M_sat/(M_p+M_sat)) #dynamical mass ratio for reduced mass factor
a_m = 60.*R_p #semimajor axis of moon (AU)
a_Roche = 2.44*R_E*(rho_E/rho_moon)**(1./3.)


#initial conditions
n_m0 = np.sqrt(G*(M_p+M_sat)/a_m**3) #mean motion of moon in rad/yr
n_p0 = np.sqrt(G*M_star/a_p**3) #mean motion of planet in rad/yr
O_p0 = 2.*np.pi/(1./365.25) #rotation rate of planet in rad/yr

n_crit = np.sqrt(3./(0.4061*(1.-1.1257*e_p))**3)*n_p0
n_Roche = np.sqrt(G*(M_p+M_sat)/a_Roche**3)
print(n_crit,n_Roche)

#tidal factors used for derivative functions SBO12 Eqn 10
#dnm_fact = -(9./2.)*(k2p*R_p**5/Q_p)*(G*M_sat)/(G*M_p)**(8./3.)  SBO12 Eqn 10
#dnp_fact = -(9./2.)*(k2p*R_p**5/Q_p)/((G*M_p)*(G*M_star)**(2./3.)) SBO12 Eqn 10
dnm_fact = -(9./2.)*(k2p*R_p**5/Q_p)*(M_sat/M_p)/(G*(M_p+M_sat))**(5./3.)*(1.-mu) #updated for larger mass ratios
dnp_fact = -(9./2.)*(k2p*R_p**5/Q_p)/((G*(M_p+M_sat))*(G*(M_star+M_p+M_sat))**(2./3.)) #updated for larger mass ratios
dOp_fact1 = -(3./2.)*(k2p*R_p**3/Q_p)*(G*M_sat)**2/(alpha*(G*M_p)**3)
dOp_fact2 = -(3./2.)*(k2p*R_p**3/Q_p)/(alpha*(G*M_p))

#tidal factors used for derivative functions SBO12 Eqn 14b
#dnm_lock1 = -M_p*(G*M_star)**(2./3.)
#dnm_lock2 = M_sat*(G*M_p)**(2./3.)
#dnm_lock3 = -3*alpha*M_p*R_p**2

dnm_lock1 = -M_p*(G*M_star)**(2./3.)/3. 
dnm_lock2 = (1.-mu)*M_sat*(G*(M_p+M_sat))**(2./3.)/3. #updated for larger mass ratios
dnm_lock3 = -alpha*R_p**2*M_p



def deriv(t,y): 
    #y = [n_m,n_p,O_p]
    sign_Op_nm = np.sign(y[2]-y[0])#(y[2]-y[0])/np.abs(y[2]-y[0])
    sign_Op_np = np.sign(y[2]-y[1])#(y[2]-y[1])/np.abs(y[2]-y[1])

    dnmdt = dnm_fact*y[0]**(16./3.)*sign_Op_nm
    dnpdt = dnp_fact*y[1]**(16./3.)*sign_Op_np
    dOpdt = dOp_fact1*y[0]**4*sign_Op_nm + dOp_fact2*y[1]**4*sign_Op_np
    if sign_Op_nm <= 0:
        sign_nm_np = np.sign(y[0]-y[1])
        dnpdt *= sign_nm_np
        dnmdt = dnm_lock1*y[1]**(-4./3.)*dnpdt/(dnm_lock2*y[0]**(-4./3.)+dnm_lock3)
        dOpdt = dnmdt
    if sign_Op_np <= 0:
        sign_np_nm = np.sign(y[1]-y[0])
        dnmdt *= sign_nm_np
        dnpdt = dnm_lock2*y[0]**(-4./3.)*dnmdt/(dnm_lock1*y[1]**(-4./3.)-dnm_lock3)
        dOpdt = dnpdt

    return [dnmdt,dnpdt,dOpdt]



t_fin = 1e11 #max lifetime considered
times = np.concatenate((np.arange(0,1000,100),np.arange(1000,1e6,1000),np.arange(1e6,t_fin+1e6,1e6))) #integration intervals
nsteps = len(times) #total number of steps
#set initial parameters
n_m,n_p,O_p = n_m0,n_p0,O_p0


T1 = 0
T = 0
if run_sim == 'T':
    fname = "output.txt"
    out = open(fname,'w')
    out.write("#time,n_sat,n_p,O_p,n_crit\n")
    out.close()
    moon_sync = False
    for i in range(0,nsteps-1):
        t_i = times[i]
        t_n = times[i+1]
        n_crit = np.sqrt(3./(0.4061*(1.-1.1257*e_p))**3)*n_p #calculate n_crit using RF2020
        a_sat = (G*(M_p+M_sat)/n_m**2)**(1./3.)/R_H #calculate satellite semimajor axis
        #output current state with time, n_sat, n_p, Omega_p, n_crit, a_sat
        out = open(fname,'a')
        out.write("%1.5e,%2.8e,%2.8e,%2.8e,%2.8e,%1.8f\n" % (t_i,n_m,n_p,O_p,n_crit,a_sat))
        out.close()
        
        delta_t = [t_i,t_n]
        temp = solve_ivp(deriv,delta_t,[n_m,n_p,O_p],method='RK45',atol=1e-12,rtol=1e-12)
        if (O_p-n_m)<0 and not moon_sync:
            T1 = t_i
            print("T1 =",t_i/1e9) #print moon synchronization time in Gyr
            moon_sync = True
        if np.abs(temp.y[0,-1]-n_m)<1e-8: #check to for collision
            break
        n_m,n_p,O_p = temp.y[0,-1],temp.y[1,-1],temp.y[2,-1]
        
        if n_m < n_crit: #check if lost due to outward migration
            break
        if n_m > n_Roche: #check if inside Roche limit
            break
    print("T =", t_i/1e9) #print moon escape/collision in Gyr
    T = t_i

####Plotting results
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker
import os

rcParams.update({'font.size': 22})

home = os.getcwd() + "/"

aspect = 16./9.
width = 12.
lw = 5
ms = 20
fs = 'x-large'

fig = plt.figure(1,figsize=(aspect*width,width),dpi=300)
ax = fig.add_subplot(111)
ax_top = ax.twiny()
ax_top.minorticks_on()
ax_top.tick_params(which='major',axis='both', direction='out',length = 12.0, width = 8.0,labelsize=fs)
ax_top.tick_params(which='minor',axis='both', direction='out',length = 8.0, width = 6.0)

ax_top.set_xlim(0,7)


#ax.axhline(n_crit,0,1,color='k',linestyle='-.',lw=lw)

fname = "output.txt"
data = np.genfromtxt(home+fname,delimiter=',',comments='#')
merge_idx = np.where(np.abs(data[:,3]-data[:,1])<=0)[0][0]

if T1 == 0:
    T1 = data[merge_idx,0]
if T == 0:
    T = data[-1,0]

ax.plot(data[:merge_idx,0]/1e10,data[:merge_idx,1],'-',color='k',lw=lw,label='$n_{sat}$')
ax.plot(data[:merge_idx,0]/1e10,data[:merge_idx,3],'--',color='k',lw=lw,label='$\Omega_p$')
ax.plot(data[:,0]/1e10,data[:,2],'-.',color='k',lw=lw,label='$n_p$')
ax.plot(data[merge_idx:,0]/1e10,data[merge_idx:,1],'-',color='k',lw=2*lw)
ax.plot(data[merge_idx:,0]/1e10,data[merge_idx:,3],':',color='w',lw=lw)

ax.legend(loc='best',fontsize=fs)

ax.text(0.25,0.14,r'$\Omega_{\rm p} = n_{\rm sat}$', color='k',fontsize='x-large',rotation=10,weight='bold',horizontalalignment='left',transform=ax.transAxes)
# ax.text(0.74,0.38,r'$\Omega_{\rm p} = n_{\rm sat}$', color='b',fontsize='x-large',rotation=75,weight='bold',horizontalalignment='left',transform=ax.transAxes)
# ax.text(0.74,0.85,r'$\Omega_{\rm p} = n_{\rm sat}$', color='m',fontsize='x-large',rotation=45,weight='bold',horizontalalignment='right',transform=ax.transAxes)
# ax.text(0.88,0.85,r'$\Omega_{\rm p} = n_{\rm sat}$', color='c',fontsize='x-large',rotation=45,weight='bold',horizontalalignment='right',transform=ax.transAxes)

ax.minorticks_on()
ax.tick_params(which='major',axis='both', direction='out',length = 12.0, width = 8.0,labelsize=fs)
ax.tick_params(which='minor',axis='both', direction='out',length = 8.0, width = 6.0)


ax.axvline(T1/1e10,color='r',linestyle='-',lw=lw)
ax.axvline(T/1e10,color='r',linestyle='-',lw=lw)

ax.set_xlim(0,7)
ax.set_ylim(0,700)


ax.set_xlabel("Time ( $\\times\;10^{10}$ yr)",fontsize=fs)
ax.set_ylabel("Angular frequency (rad/yr)",fontsize=fs)

fig.savefig("SBO12_tide_evol.png",bbox_inches='tight',dpi=300)
plt.close()