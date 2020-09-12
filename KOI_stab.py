import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.stats as st
from matplotlib import rcParams


rcParams.update({'font.size': 22})

def get_truncated_normal(mean=0,std=1, low=0, upp=10):
    a,b = -low/std,upp/std
    X = st.truncnorm.rvs(a,b, loc=mean, scale=1)
    return X

def RF_stab(ep,C_1,C_2):
    return C_1*(1.-C_2*ep)

def get_radius(R_p,M_star,a_p,KOI_mass_pl):
    KOI_RH = a_p*(KOI_mass_pl/(3.*M_star))**(1./3.)
    return KOI_RH

M_E = 3.0035e-6 #M_sun
R_E = 4.2635e-5 #AU
rho_E = 5.515 # g/cm^3
rho_N = 1.64 # g/cm^3
rho_J = 1. # g/cm^3
rho_moon = 3.34 #g/cm^3
rho_Mars = 3.93 #g/cm^3

#KOI planet radius
KOI_Rp = [3.32,2.78,4.54,1.10,3.224,5.559]
KOI_Rp_high = [0.85,0.39,0.34,0.05,0.213,0.252]
KOI_Rp_low = [0.64,0.38,0.31,0.04,0.159,0.889]
#KOI stellar mass
KOI_Ms = [1.175,0.871,1.406,0.890,1.45,1.34]
KOI_Ms_high = [0.058,0.142,0.086,0.009,0.601,0.054]
KOI_Ms_low = [0.064,0.142,0.086,0.011,0.271,0.051]
#KOI planet semimajor axis
KOI_ap = [0.4756,0.2897,0.5337,0.3183,0.2743,0.4039]

#Probabalistic Values from Chen & Kipping (2017)
KOI_Mp = [11.1,8.13,18.6,1.37,10.4,25.5]
KOI_Mp_high = [10.2,6.7,16.4,0.88,9.,24.2]
KOI_Mp_low = [5.5,3.67,8.4,0.44,4.71,12.6]

KOI_num = [268,303,1888,1925,2728,3220]

aspect = 2.
width = 12.
lw = 3
ms = 20
fs = 'x-large'

fig = plt.figure(1,figsize=(aspect*width,2*width),dpi=300)
ax11 = fig.add_subplot(321)
ax12 = fig.add_subplot(322)
ax21 = fig.add_subplot(323)
ax22 = fig.add_subplot(324)
ax31 = fig.add_subplot(325)
ax32 = fig.add_subplot(326)

ax_list = [ax11,ax12,ax21,ax22,ax31,ax32]
ep = np.arange(0,0.91,0.01)

yticks = [1,2,3,4,5,6,7,8,9,10,20,30,40,50]
yticklabels = ["","2","","4","","6","","8","","10","20","30","40","50"]
xticks = np.arange(0,1.,0.1)
KOI_stab_rvs = np.zeros(1000)

for i in range(0,6):
    ax = ax_list[i]
    KOI_Rpi = KOI_Rp[i]*R_E
    KOI_RH = get_radius(KOI_Rp[i],KOI_Ms[i],KOI_ap[i],KOI_Mp[i]*M_E)
    KOI_stab = RF_stab(ep,0.4061,1.1257)*KOI_RH/KOI_Rpi
    KOI_dens = KOI_Mp[i]/KOI_Rp[i]**3
    if i == 3:
        KOI_dens *= rho_E
    else:
        KOI_dens *= 3.829**3/17.15*rho_N
    Roche_rad = 2.44*(KOI_dens/rho_Mars)**(1./3.)
    print(Roche_rad)


    ax.fill_between(ep,KOI_stab,2*np.ones(len(ep)),KOI_stab>2*np.ones(len(ep)),color='k')
    for j in range(0,1000):
        KOI_Rp_rvs = get_truncated_normal(mean=KOI_Rp[i]*R_E,low=KOI_Rp_low[i]*R_E,upp=KOI_Rp_high[i]*R_E)
        KOI_Ms_rvs = get_truncated_normal(mean=KOI_Ms[i],low=KOI_Ms_low[i],upp=KOI_Ms_high[i])
        KOI_Mp_rvs = get_truncated_normal(mean=KOI_Mp[i]*M_E,low=KOI_Mp_low[i]*M_E,upp=KOI_Mp_high[i]*M_E)
        KOI_RH_rvs = get_radius(KOI_Rp_rvs,KOI_Ms_rvs,KOI_ap[i],KOI_Mp_rvs)
        if i == 3:
            KOI_stab_rvs[j] = RF_stab(0.62,0.4061,1.1257)*KOI_RH_rvs/KOI_Rp_rvs
        ax.plot(ep,RF_stab(ep,0.4061,1.1257)*KOI_RH_rvs/KOI_Rp_rvs,'-',color='gray',lw=lw,alpha=0.1)
    if i == 3:
        ax.text(0.25,0.01,'Roche Limit', color='k',fontsize='x-large',weight='bold',horizontalalignment='center',transform=ax.transAxes)
    ax.plot(ep,KOI_stab,'r-',lw=lw)
    ax.plot(ep,Roche_rad*np.ones(len(ep)),'w--',lw=lw)
    

    ax_list[i].text(0.99,0.9,'KOI %i.01' % KOI_num[i], color='k',fontsize='x-large',weight='bold',horizontalalignment='right',transform=ax.transAxes)
    ax.minorticks_on()
    ax.tick_params(which='major',axis='both', direction='out',length = 12.0, width = 8.0,labelsize=fs)
    ax.tick_params(which='minor',axis='both', direction='out',length = 6.0, width = 4.0)

    if i > 3:
        ax.set_xlabel("$e_{\\rm p}$",fontsize=fs)
    if i % 2 == 0:
        ax.set_ylabel("$a_{\\rm sat}$ (R$_{\\rm p}$)",fontsize=fs)
    ax.set_yscale('log')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)
    ax.set_xticks(xticks)
    ax.set_ylim(2,50)
    ax.set_xlim(0.,0.9)
    

fig.savefig("FW_KOIs_stab.png",bbox_inches='tight',dpi=300)
plt.close()
