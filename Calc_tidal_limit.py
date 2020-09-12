import numpy as np
import sys
import os

def get_Crit_Q_QLRF(k2p,a_crit,a_m,R_p,M_m,M_p,M_s,T): #Quarles, Li, & Rosario-Franco (2020) Eqn 6
    n_m = np.sqrt(G*(M_p+M_m)/a_m**3)
    n_crit = np.sqrt(G*(M_p+M_m)/a_crit**3)
    top = (39./2.)*(k2p*R_p**5*T)*(G*M_m)
    bot = (G*M_p)*(G*(M_p+M_m))**(5./3.)*(n_crit**(-13./3.) - n_m**(-13./3.))
    return top/bot

def get_Crit_Q_CLP(k2p,k2m,a_crit,a_m,R_m,R_p,M_m,M_p,T): #Cheng, Lee, & Peale (2014) Eqns 16 & 19
    if (R_m/R_E) < 2:
        Q_m = 100
    else:
        Q_m = 1000
    n_m = np.sqrt(G*(M_p+M_m)/a_m**3)
    n_c = np.sqrt(G*(M_p+M_m)/a_crit**3)
    mr = M_p/(M_p+M_m)
    top = k2p*R_p**5*M_m**2
    bot1 = (2./39)*(n_c**(-13./3.)-n_m**(-13./3.))*mr*M_m*(G*(M_p+M_m))**(8./3.)/(G*T)
    bot2 = k2m*R_m**5*M_p**2/Q_m
    return top/bot1

G = 4.*np.pi**2
M_sun = 1.989e33 #Mass of Sun in g
AU = 149597870700. #AU in m
AU_km = AU/1000.
AU_cm = AU*100.

M_E = 3.0035e-6 #M_sun
R_E = 4.2635e-5 #AU
rho_E = 5.515 # g/cm^3
rho_N = 1.64 # g/cm^3
rho_J = 1. # g/c

f_crit = 0.4061 #C_1 from Rosario-Franco+ (2020)
#KOI planet radius
KOI_Rp = [3.32,2.78,4.54,1.10,3.224,5.559,10]
KOI_Rp_high = [0.85,0.39,0.34,0.05,0.213,0.252]
KOI_Rp_low = [0.64,0.38,0.31,0.04,0.159,0.889]
#KOI stellar mass
KOI_Ms = [1.175,0.871,1.406,0.890,1.45,1.34,10.4]
KOI_Ms_high = [0.058,0.142,0.086,0.009,0.601,0.054]
KOI_Ms_low = [0.064,0.142,0.086,0.011,0.271,0.051]
#KOI planet semimajor axis
KOI_ap = [0.4756,0.2897,0.5337,0.3183,0.2743,0.4039,0.87]
#KOI system age
KOI_Tp = [3.05,6.31,1.26,6.98,1.7,1.7,10]
KOI_Tp_high = [0.85,3.15,0.33,0.4,0.53,0.556]
KOI_Tp_low = [0.64,3.81,0.18,0.5,0.392,0.459]

#Probabalistic Masses from Chen & Kipping (2017)
KOI_Mp = [11.1,8.13,18.6,1.37,10.4,25.5,4*333]
KOI_Mp_high = [10.2,6.7,16.4,0.88,9.,24.2]
KOI_Mp_low = [5.5,3.67,8.4,0.44,4.71,12.6]

KOI_num = [268,303,1888,1925,2728,3220,1625]
M_frac = np.concatenate((np.arange(0.001,0.01,0.001),np.arange(0.01,0.1,0.005),np.arange(0.1,1.05,0.05)))
a_sat = np.concatenate((np.arange(2,10,0.5),np.arange(10,100,1),np.arange(100,300,50)))
Q_rng = [10,100,1000,10000,100000]

k2p_rng = [0.12,0.29,0.53]

home = os.getcwd() + "/"
model = sys.argv[1]
if not os.path.exists(home+model):
    os.makedirs(home+model)
fname = "Moon_lifetimes_"

for i in range(0,7):
    out = open(home+fname + "KOI%i.txt" % KOI_num[i],'w')
    out.write("#a_s,M_f,Q_p,T,stab\n")
    out.close()
    if KOI_Rp[i] < 2:
        k2p = 0.299  #Terrestrial
    elif KOI_Rp[i] > 8:
        k2p = 0.53   #Jovian
    else:
        k2p = 0.12   #Neptune
    for M_f in M_frac:
        for a_s in a_sat:
            R_p = KOI_Rp[i]*R_E
            M_p = KOI_Mp[i]*M_E
            R_H = KOI_ap[i]*((1.+M_f)*M_p/(3.*KOI_Ms[i]))**(1./3.)
            a_c = f_crit*R_H
            T_sys = (KOI_Tp[i])*1e9
            if M_f*M_p < (4*M_E):
                k2m = 0.299
                R_m = (M_f*M_p/M_E)**(1./3.)*R_E
            else:
                k2m = 0.12
                R_m = (M_f*M_p/M_E*(5.515/1.64))**(1./3.)*R_E

            if M_f >= 0.1:
                Q_c = get_Crit_Q_CLP(k2p,k2m,a_c,a_s*R_p,R_m,R_p,M_f*M_p,M_p,T_sys)
            else:
                Q_c = get_Crit_Q_QLRF(k2p,a_c,a_s*R_p,R_p,M_f*M_p,M_p,KOI_Ms[i],T_sys)
            if (a_s*R_p) > a_c:
                Q_c = -1 #satellite is unstable using RF2020
            out = open(home+fname + "KOI%i.txt" % KOI_num[i],'a')
            out.write("%1.1f,%1.3f,%1.3e\n" % (a_s,M_f,Q_c))
            out.close()
                