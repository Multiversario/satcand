import numpy as np
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib import ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import griddata
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from matplotlib import rcParams
import os
import sys

def func(x,a,b,c):
    return a*x**2+b*x+c

def get_root(a,b,c,y0):
    return (-b - np.sqrt(b**2-4.*a*(c-y0)))/(2.*a)

def colorbar(mappable):
    last_axes = plt.gca()
    ax_m = mappable.axes
    fig1 = ax_m.figure
    divider = make_axes_locatable(ax_m)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig1.colorbar(mappable, cax=cax,orientation='vertical')
    cbar.ax.yaxis.set_ticks_position('right')
    plt.sca(last_axes)
    return cbar

home = os.getcwd() + "/"
fname = "/Moon_lifetimes_"

vmin = 10
vmax = 100
cmap = cm.gist_rainbow
my_cmap = cm.get_cmap(cmap)
norm = colors.Normalize(vmin,vmax)
cmmapable = cm.ScalarMappable(norm,my_cmap)
cmmapable.set_array(range(0,1))
my_cmap.set_under('w')
my_cmap.set_over('k')

#estimated values from Figure 3 3sigma curves; Kipping (2020)
KOI268_Mf = [1,0.6,0.4,0.3,0.2,0.1,0.06,0.03,0.01]
KOI268_ap = [3,4,5,6,9,16,30,50,120]

KOI303_ap = [4.2,6,9,22,50,100,250]
KOI303_Mf = [1,0.6,0.3,0.1,0.045,0.022,0.0075]

KOI1888_ap = [5,6,7,10,23,70,200]
KOI1888_Mf = [1,0.65,0.5,0.3,0.1,0.03,0.0082]

KOI1925_ap = [10,20,50,100,200,250]
KOI1925_Mf = [1,0.3,0.1,0.05,0.022,0.018]

KOI2728_ap = [6,8,10,20,30,50,100,250]
KOI2728_Mf = [1,0.6,0.45,0.2,0.1,0.065,0.03,0.01]

KOI3220_ap = [2,5,10,20,50,100,200]
KOI3220_Mf = [0.28,0.1,0.04,0.02,0.008,0.003,0.002]

KOI_ap = [KOI268_ap,KOI303_ap,KOI1888_ap,KOI1925_ap,KOI2728_ap,KOI3220_ap]
KOI_Mf = [KOI268_Mf,KOI303_Mf,KOI1888_Mf,KOI1925_Mf,KOI2728_Mf,KOI3220_Mf]

KOI_Ms = [1.175,0.871,1.406,0.890,1.45,1.34]
KOI_Mp = [11.1,8.13,18.6,1.37,10.4,25.5]
KOI_Rp = [3.32,2.78,4.54,1.10,3.224,5.559]
KOI_ap_med = [0.4756,0.2897,0.5337,0.3183,0.2743,0.4039]
f_crit = 0.4061

M_E = 3.0035e-6 #M_sun
R_E = 4.2635e-5 #AU
rho_E = 5.515 # g/cm^3
rho_moon = 3.3 # g/cm^3
rho_Mars = 3.93 #g/cm^3
rho_N = 1.64 # g/cm^3
rho_J = 1. # g/c

rcParams.update({'font.size': 22})

KOI_num = [268,303,1888,1925,2728,3220]

aspect = 2.
width = 12.
lw = 8
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


Mm_rng = np.concatenate((np.arange(0.001,0.01,0.001),np.arange(0.01,0.1,0.005),np.arange(0.1,1.05,0.05)))
yi = np.concatenate((np.arange(0.001,0.01,0.001),np.arange(0.01,0.105,0.005)))
xi = np.concatenate((np.arange(2,10,0.5),np.arange(10,100,1),np.arange(100,300,50)))
Q_rng = [10,100,1000,10000,100000]
ymax = [0.57,0.68,0.64,0.74,0.8,0.44]
color = ['k','r','g','b','c']


for i in range(0,6):
    ax = ax_list[i]
    dat_file = np.genfromtxt(home+fname+"KOI%i.txt" % KOI_num[i],delimiter=',',comments='#')
    popt, pcov = curve_fit(func, np.log10(KOI_ap[i]), np.log10(KOI_Mf[i])) #Fit curves to Kipping (2020) Figure 3

    y_obs = 10**func(np.log10(dat_file[:,0]),*popt)
    Obs_cut = np.where(dat_file[:,1]>y_obs)[0]
    dat_file[Obs_cut,2] = -1
    if i == 3:
        Q_cut = np.where(dat_file[:,2]>200)[0]
    else:
        Q_cut = np.where(dat_file[:,2]>2000)[0]
    dat_file[Q_cut,2] = -1
    dat_file[:,2] = np.ma.masked_where(dat_file[:,2] <= 0, dat_file[:,2])
    Q_cut = np.where(dat_file[:,2]>0)[0]
    dat_file[Q_cut,2] = 1000
    
    y_obs = 10**func(np.log10(xi),*popt)
    ax.plot(xi,y_obs,'k-',lw=lw,zorder=5)
    
    R_H = KOI_ap_med[i]*(KOI_Mp[i]*M_E/(3.*KOI_Ms[i]))**(1./3.)
    a_c = f_crit*R_H/(KOI_Rp[i]*R_E)
    ax.fill_betweenx(Mm_rng,2,a_c,color='b')
    zi = griddata((dat_file[:,0],dat_file[:,1]),(dat_file[:,2]),(xi[None,:],yi[:,None]),method = 'linear')
    CS = ax.pcolor(xi,yi,zi,cmap = cmap,vmin=vmin,vmax=vmax,norm=colors.LogNorm(vmin=10, vmax=100))
    
    xmax = 10**get_root(*popt,-1)

    #assume R_m = 0.5 R_E
    R_m = (Mm_rng*KOI_Mp[i]/rho_Mars)**(1./3.)
    Rm_cut = np.where(R_m>0.5)[0]

    if len(Rm_cut)>1:
        Rm_xmax = 10**get_root(*popt,np.log10(Mm_rng[Rm_cut[0]]))
        ax.fill_betweenx(np.arange(0.1,1.1,0.1),2,xmax,color='gray')
        ax.fill_betweenx(Mm_rng,a_c,250,color='r')
        ax.fill_betweenx(Mm_rng[Rm_cut],2,Rm_xmax,facecolor='none',edgecolor='w',hatch='X',lw=0)
        
        ax.axvline(a_c-0.25,0,ymax[i],color='c',linestyle='--',lw=8)

        
    else:
        xmax = 10**get_root(*popt,-1)
        ax.fill_betweenx(np.arange(0.1,1.1,0.1),2,xmax,color='gray')
        ax.fill_betweenx(Mm_rng,a_c,250,color='r')
    #y_stab = 10**func(np.log10(a_c),*popt)
    
    x_rng = np.arange(2,250,0.01)
    y_obs = 10**func(np.log10(x_rng),*popt)
    tide_limit = 0.1*len(x_rng)
    ax.fill_between(x_rng,y_obs,tide_limit,y_obs<tide_limit,color='w')

    ax_list[i].text(0.99,0.9,'KOI %i.01' % KOI_num[i], color='k',fontsize='x-large',weight='bold',horizontalalignment='right',transform=ax.transAxes)
    ax.minorticks_on()
    ax.tick_params(which='major',axis='both', direction='out',length = 12.0, width = 8.0,labelsize=fs)
    ax.tick_params(which='minor',axis='both', direction='out',length = 8.0, width = 6.0)
    if i > 3:
        ax.set_xlabel("$a_{\\rm sat}$  (R$_p$)",fontsize=fs)
    if i % 2 == 0:
        ax.set_ylabel("M$_{\\rm sat}$/M$_{\\rm p}$",fontsize=fs)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_ylim(0.001,1)
    ax.set_xlim(2.,250)
        

fig.savefig("FW_KOIs_limits.png",bbox_inches='tight',dpi=300)
plt.close()
