from scipy.integrate import solve_ivp
import numpy as np

def n_m(t):
    return (T_fact_inv*(G*M_sat)*t*(G*M_p)**(-8./3.) + n_m0**(-13./3.))**(-3./13.)
def n_p(t):
    return (T_fact_inv*t/((G*M_p)*(G*M_star)**(2./3.)) + n_p0**(-13./3.))**(-3./13.)
    
def deriv(t,y): 
    #y = [n_m,n_p,O_p]
    sign_Op_nm = (y[2]-y[0])/np.abs(y[2]-y[0])
    sign_Op_np = (y[2]-y[1])/np.abs(y[2]-y[1])

    dnmdt = dnm_fact*y[0]**(16./3.)
    dnpdt = dnp_fact*y[1]**(16./3.)
    dOpdt = dOp_fact1*y[0]**4*sign_Op_nm + dOp_fact2*y[1]**4*sign_Op_np
    if sign_Op_nm <= 0:
        sign_nm_np = (y[0]-y[1])/np.abs(y[0]-y[1])
        dnpdt *= sign_nm_np
        dnmdt = dnm_lock1*y[1]**(-4./3.)*dnpdt/(dnm_lock2*y[0]**(-4./3.)+dnm_lock3)  #Differential form of Eqn 14b SBO12
        dOpdt = dnmdt
    if sign_Op_np <= 0:
        sign_np_nm = (y[1]-y[0])/np.abs(y[1]-y[0])
        dnmdt *= sign_nm_np
        dnpdt = dnm_lock2*y[0]**(-4./3.)*dnmdt/(dnm_lock1*y[1]**(-4./3.)-dnm_lock3)  #Differential form of Eqn 15b SBO12
        dOpdt = dnpdt
    return [dnmdt,dnpdt,dOpdt]
    
def check_Eqn_28(T1):
    T_28_num  = ((G*M_p)*(G*M_sat)**(7./6.)*(G*M_star)**(2./3.))*(n_p0**(-13./3.)-(M_p/M_sat)**(13./6.)*n_m0**(-13./3.))
    T_28_denom = (G*M_p)**(1./2.)*(G*M_star)**(2./3.)-(G*M_sat)**(7./6.)
    T_28 = T_fact*T_28_num/T_28_denom
    if T1 <= T_28:
        return True
    else:
        return False
def check_Eqn_35(T1):
    T_35_num  = (G*M_p)**(8./3.)*(G*M_star)**(2./3.)*((f_crit**3/3.)**(13./6.)*n_p0**(-13./3.)-n_m0**(-13./3.))
    T_35_denom = (G*M_sat)*(G*M_star)**(2./3.)-(f_crit**3/3.)**(13./6.)*(G*M_p)**(5./3.)
    T_35 = T_fact*T_35_num/T_35_denom
    if T1 > T_35:
        return True
    else:
        return False
def check_Eqn_32(T1):
    X_32  = ((M_p/M_sat)**(13./6.)*n_m(T1)**(-13./3.)+ (G*M_p)**(1./2.)*(G*M_star)**(2./3.)*n_p(T1)**(-13./3.)/(G*M_sat)**(7./6.))**(1./13.)
    c_32 = 1./(alpha*R_p**2*(G*M_p))
    a_32_1 = (G*M_p)**(1/2.)*(G*M_sat)**(7./8.)
    a_32_2 = (G*M_sat)**(-7./26.)
    b_32 = ((G*M_sat)**(7./6.)+(G*M_p)**(1./2.)*(G*M_star)**(2./3.))**(1./13.)
    T_32  = a_32_1*b_32**12*X_32**4 - G*L0*X_32**3 + a_32_2*b_32**3/c_32
    if T_32 <=0:
        return True
    else:
        return False

def get_Moon_type(n_m0,n_p0,O_p0):
    #Find T1 by integrating Eqn 10 SBO12
    y_old = [n_m0,n_p0,O_p0]
    t_step = 1e7
    t_fin = 1e11
    nsteps = int(t_fin/t_step)
    T1 = 0
    for i in range(0,nsteps):
        t_i = i*t_step
        t_n = (i+1)*t_step
        delta_t = [t_i,t_n]

        temp = solve_ivp(deriv,delta_t,y_old,method='RK45',atol=1e-10,rtol=1e-10)
        y_new = [temp.y[0,-1],temp.y[1,-1],temp.y[2,-1]]
        if (y_old[2]-y_old[0])*(y_new[2]-y_new[0]) < 0:
            while (t_n - t_i)>1e6:
                m_o = 0.5*(t_i+t_n)
                delta_t = [t_i,m_o]
                temp = solve_ivp(deriv,delta_t,y_old,method='RK45',atol=1e-10,rtol=1e-10)
                y_s = [temp.y[0,-1],temp.y[1,-1],temp.y[2,-1]]
                if (y_old[2]-y_old[0])*(y_s[2]-y_s[0])<0:
                    t_n = m_o
                elif (y_new[2]-y_new[0])*(y_s[2]-y_s[0])<0:
                    t_i = m_o
            T1 = m_o
            break
        else:
            y_old = y_new

    print("T1 =",T1/1e9)

    T_5 = (3**(3./4.)*G*L0-4*((G*M_sat)**3*(G*M_p)**3*alpha*R_p**2)**(1./4.))**13*((G*M_sat)**(7./6.)+(G*M_p)**(1./2.)*(G*M_star)**(2./3.))**12
    T_6 = 3**(39./4.)*(G*M_p)**6*(G*M_star)**8*(G*L0)**13
    #SBO12 Decision Tree (Fig 17)
    if check_Eqn_35(T1):
        #Type IV Moon lifetime
        T_26_num = (G*M_p)**(8./3.)*((G*M_p/(f_crit*R_H)**3)**(-13./6.)-n_m0**(-13./3.))
        T_26 = T_fact*T_26_num/(G*M_sat)
        print("Type IV: lifetime = %1.3f x 10^10" % (T_26/1e10))
    elif check_Eqn_28(T1):
        #Type I Moon lifetime
        T_21_num = 3**(3./4.)*G*L0-4*((G*M_sat)**3*(G*M_p)**3*alpha*R_p**2)**(1./4.)
        T_21_denom = 3**(3./4.)*(G*M_p)*(G*M_star)**(2./3.)
        T_21 = T_fact*(G*M_p)*(G*M_star)**(2./3.)*((T_21_num/T_21_denom)**13-n_p0**(-13./3.))
        print("Type I Case 1: lifetime = %1.3f x 10^10" % (T_21/1e10))
    elif check_Eqn_32(T1):
        #Type I Moon lifetime
        T_21_num = 3**(3./4.)*G*L0-4*((G*M_sat)**3*(G*M_p)**3*alpha*R_p**2)**(1./4.)
        T_21_denom = 3**(3./4.)*(G*M_p)*(G*M_star)**(2./3.)
        T_21 = T_fact*(G*M_p)*(G*M_star)**(2./3.)*((T_21_num/T_21_denom)**13-n_p0**(-13./3.))
        print("Type I Case 2: lifetime = %1.3f x 10^10" % (T_21/1e10))
    elif T_5 > T_6:
        #Type II Moon lifetime
        T_22_A = (G*M_p)**(8./3.)/(G*M_sat)*n_m0**(-13./3.)
        T_22_B = (3**(3./4.)*G*L0-4*((G*M_sat)**3*(G*M_p)**3*alpha*R_p**2)**(1./4.))**13/(3**(39./4.)*(G*M_p)**12*(G*M_star)**8)
        T_22_C = (G*L0)**13/((G*M_p)**(1./2.)*(G*M_sat)**(7./6.)+(G*M_p)*(G*M_star)**(2./3.))**12
        T_22 = T_fact*(T_22_A+T_22_B-T_22_C)
        print("Type II: lifetime = %1.3f x 10^10" % ((2*T1 + T_22)/1e10))
    elif T_5 <= T_6:
        #Type III Moon lifetime
        T_24 = T_fact*(G*M_p)**(8./3.)/(G*M_sat)*n_m0**(-13./3.)
        print("Type III: lifetime = %1.3f x 10^10" % ((2*T1 + T_24)/1e10))

#constants
G = 4.*np.pi**2
M_E = 3.0035e-6 #M_sun
R_E = 4.2635e-5 #AU
rho_E = 5.515 # g/cm^3
rho_N = 1.64 # g/cm^3
rho_J = 1. # g/cm^3
rho_moon = 3.34 #g/cm^3

#assumed planet parameters
alpha = 0.3308  #moment of inertia
k2p = 0.299 #planet Love number
M_star = 1. #mass of host star (in M_sun)
Q_p = 12 #tidal Quality factor
R_p = 1.0*R_E  #Radius of planet (AU)
rho_p = rho_E
M_p = 1.*M_E#(rho_p/rho_E)*(R_p/R_E)**3*M_E  #Mass of planet assuming Neptune-like density
a_p = 1.0 # semimajor axis of planet (AU)

M_sat = 0.0123*M_E #Mass of Moon (M_sun)
R_H = a_p*((M_p+M_sat)/(3*M_star))**(1./3.)
mu = (M_sat/(M_p+M_sat)) #dynamical mass ratio for reduced mass factor
a_m = 60.*R_p #semimajor axis of moon (AU)
f_crit = 0.4061 #1.2**(-2./3.)*0.36 for SBO12 Fig 7

#initial conditions
n_m0 = np.sqrt(G*M_p/a_m**3) #mean motion of moon in rad/yr
n_p0 = 2.*np.pi#np.sqrt(M_star/a_p**3) #mean motion of planet in rad/yr
O_p0 = 2.*np.pi/(1./365.25) #rotation rate of planet in rad/yr
L0 = M_sat*(G*(M_p+M_sat))**(2./3.)*n_m0**(-1./3.) + alpha*R_p**2*M_p*O_p0 + M_p*(G*M_star)**(2./3.)*n_p0**(-1./3.)

#factors for differential equations
T_fact  = (2./39.)*(Q_p/(k2p*R_p**5))
T_fact_inv = (39./2.)*(k2p*R_p**5)/Q_p

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

get_Moon_type(n_m0,n_p0,O_p0)

