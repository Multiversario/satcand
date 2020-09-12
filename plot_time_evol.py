import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.ticker
import os

rcParams.update({'font.size': 22})

home = os.getcwd() + "/"
Q_rng = [10,100]
Mf_rng = [0.0123,0.3]
color = ['r','m','b','c']

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

ax_top.set_xscale('log')

ax_top.set_xlim(1e2,1e10)
locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
ax_top.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=np.arange(0,1.1,0.1),numticks=12)
ax_top.xaxis.set_minor_locator(locmin)
ax_top.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

ax.text(0.06,0.25,r'$e_p = 0.0$', color='k',fontsize='x-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.text(0.06,0.53,r'$e_p = 0.6$', color='k',fontsize='x-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.arrow(0.05, 0.23, 0., 0.1, width=0.01,head_width=0.03, head_length=0.03, fc='gray', ec='gray',transform=ax.transAxes)
ax.axhline(220,0,0.97,color='k',linestyle='-.',lw=lw)
ax.axhline(1194,0,0.97,color='k',linestyle=':',lw=lw)
cnt = 0
for Q_p in Q_rng:
    for Mf in Mf_rng:
        fname = "KOI1925/output_KOI1925_[%03i,%1.4f].txt" % (Q_p,Mf)
        data = np.genfromtxt(home+fname,delimiter=',',comments='#')
        
        merge_idx = np.where(np.abs(data[:,3]-data[:,1])<=0)[0][0]

        ax.plot(data[:merge_idx,0],data[:merge_idx,1],'-',color=color[cnt],lw=lw,label='$Q_p$ = %s, $M_{sat}/M_p$ = %s' % (Q_p,Mf))
        ax.plot(data[:merge_idx,0],data[:merge_idx,3],'--',color=color[cnt],lw=lw)
        
        ax.plot(data[merge_idx:,0],data[merge_idx:,1],'-',color=color[cnt],lw=2*lw)
        ax.plot(data[merge_idx:,0],data[merge_idx:,3],':',color='w',lw=lw)
        cnt += 1
        
ax.text(0.26,0.6,r'$n_{\rm sat}$', color='r',fontsize='x-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.text(0.38,0.6,r'$n_{\rm sat}$', color='b',fontsize='x-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.text(0.56,0.6,r'$\Omega_{\rm p}$', color='r',fontsize='x-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.text(0.685,0.6,r'$\Omega_{\rm p}$', color='b',fontsize='x-large',weight='bold',horizontalalignment='left',transform=ax.transAxes)

ax.text(0.62,0.38,r'$\Omega_{\rm p} = n_{\rm sat}$', color='r',fontsize='x-large',rotation=75,weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.text(0.74,0.38,r'$\Omega_{\rm p} = n_{\rm sat}$', color='b',fontsize='x-large',rotation=75,weight='bold',horizontalalignment='left',transform=ax.transAxes)
ax.text(0.74,0.85,r'$\Omega_{\rm p} = n_{\rm sat}$', color='m',fontsize='x-large',rotation=45,weight='bold',horizontalalignment='right',transform=ax.transAxes)
ax.text(0.88,0.85,r'$\Omega_{\rm p} = n_{\rm sat}$', color='c',fontsize='x-large',rotation=45,weight='bold',horizontalalignment='right',transform=ax.transAxes)

ax.legend(loc='best',fontsize='large',ncol=2)

ax.minorticks_on()
ax.tick_params(which='major',axis='both', direction='out',length = 12.0, width = 8.0,labelsize=fs)
ax.tick_params(which='minor',axis='both', direction='out',length = 8.0, width = 6.0)


ax.axvline(7e9,color='k',linestyle='-',lw=lw)
ax.set_yscale('log')
ax.set_xscale('log')

#ax.set_xticks([1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10])
ax.set_xlim(1e2,1e10)
ax.set_ylim(70,10000)

locmaj = matplotlib.ticker.LogLocator(base=10,numticks=12) 
ax.xaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=np.arange(0,1.1,0.1),numticks=12)
ax.xaxis.set_minor_locator(locmin)
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


ax.set_xlabel("Time (yr)",fontsize=fs)
ax.set_ylabel("Angular frequency (rad/yr)",fontsize=fs)

fig.savefig("KOI1925_tide_evol.png",bbox_inches='tight',dpi=300)
plt.close()
