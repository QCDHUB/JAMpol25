#!/usr/bin/env python
import sys, os
import numpy as np
import copy

sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import pylab as py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from tools.config import conf, load_config
from tools.tools import save, load

from analysis.qpdlib import ppdf
from analysis.obslib import pstf

#regen = True
regen = False

#ext = '.png'
ext = '.pdf'

Q2 = 10

def plot_moments():

    
    W2 = np.array([3.5,4,5,6,10])
    #W2 = np.array([3.5,4,5,10])

    flavs = ['up','dp','sp','g','gA','a8','Sigma']
    d2X = ['d2p','d2n']
    
    means, stds = {}, {}
    
    mom = 1
    #--generate data if it doesn't exist
    for i in range(len(W2)):
        if W2[i] == 3.5:  wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT3p5/pos_g'
        elif W2[i] == 10: wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noHT10/pos_g'
        else:             wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT%d/pos_g'%W2[i]
    
        load_config(wdir + '/input.py')
   
        #--fixed xmin and xmax
        xmin = 5e-3 #--fixed by minimum of data
        xmax = 0.52 #--fixed by Q2/(W2 - M2 + Q2) with W2 = 10
        #xmax = Q2/(W2[i] - 0.938**2 + Q2)
        filename = '%s/data/ppdf-moment-trunc-%d-Q2=%3.5f-xmin=%3.5f-xmax=%3.5f.dat'%(wdir,mom,Q2,xmin,xmax)
    
        if regen: ppdf.gen_moments_trunc(wdir, Q2 = Q2, flavors = flavs, mom=mom, xmin = xmin, xmax = xmax)
        try:
            moments = load(filename)
        except:
            ppdf.gen_moments_trunc(wdir, Q2 = Q2, flavors = flavs, mom=mom, xmin = xmin, xmax = xmax)
            moments = load(filename)
   
 
        for flav in flavs:
            if flav not in means: means[flav] = np.zeros(len(W2))
            if flav not in stds:  stds [flav] = np.zeros(len(W2))
            if flav=='Sigma': moment = np.array(moments['moments']['up'])[:,0] + np.array(moments['moments']['dp'])[:,0]
            else:             moment = np.array(moments['moments'][flav])[:,0]
            means[flav][i] = np.mean(moment)
            stds [flav][i] = np.std (moment)
       

        #--get d2
        #for tar in ['p','n']:
        #    if 'd2'+tar not in means: means['d2'+tar] = np.zeros(len(W2))
        #    if 'd2'+tar not in stds:  stds ['d2'+tar] = np.zeros(len(W2))
        #    filename = '%s/data/d2-trunc-%s-Q2=%3.5f-xmin=%3.5f-xmax=%3.5f.dat'%(wdir,tar,Q2,xmin,xmax)
        #    if regen: pstf.gen_d2_trunc(wdir, tar=tar, Q2 = Q2, xmin = xmin, xmax = xmax)
        #    try:
        #        moments = load(filename)
        #    except:
        #        pstf.gen_d2_trunc(wdir, tar=tar, Q2 = Q2, xmin = xmin, xmax = xmax)
        #        moments = load(filename)
    
        #    d2 = np.array(moments['moments'])[:,0]
        #    means['d2'+tar][i] = np.mean(d2,axis=0) 
        #    stds ['d2'+tar][i] = np.std (d2,axis=0)
     
    # Enter moments from ppdf_run.py here
    up_mean    = np.array(means['up'])
    dp_mean    = np.array(means['dp'])
    sp_mean    = np.array(means['sp'])
    g_mean     = np.array(means['g'])
    #d2p_mean   = np.array(means['d2p'])
    #d2n_mean   = np.array(means['d2n'])
    gA_mean    = np.array(means['gA'])
    a8_mean    = np.array(means['a8'])
    Sigma_mean = np.array(means['up'])
    
    up_std     = np.array(stds['up'])
    dp_std     = np.array(stds['dp'])
    sp_std     = np.array(stds['sp'])
    g_std      = np.array(stds['g'])
    #d2p_std    = np.array(stds['d2p'])
    #d2n_std    = np.array(stds['d2n'])
    gA_std     = np.array(stds['gA'])
    a8_std     = np.array(stds['a8'])
    Sigma_std  = np.array(stds['Sigma'])

    ncols = 1
    nrows = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*10,nrows*5))
    
    colors = ['r','b','forestgreen','magenta','darkgoldenrod']

    markersizes = [5,5,5,5]
    fmts = ['o','o','o','o']  

    dW2 = 0.10 
    capsize = 2 
    axs[0].errorbar(W2    ,up_mean,yerr=up_std,        color = colors[0],fmt = fmts[0],linestyle = 'None',capsize = capsize,markersize=markersizes[0])
    axs[0].errorbar(W2    ,dp_mean,yerr = dp_std,      color = colors[1],fmt = fmts[1],linestyle = 'None',capsize = capsize,markersize=markersizes[1])
    #axs[0].errorbar(W2    ,sp_mean,yerr = sp_std,      color = colors[2],fmt = fmts[2],linestyle = 'None',capsize = capsize,markersize=markersizes[2])
    axs[0].errorbar(W2+dW2,g_mean, yerr = g_std,       color = colors[2],fmt = fmts[3],linestyle = 'None',capsize = capsize,markersize=markersizes[3])
                                                                                                                                            
    axs[1].errorbar(W2,    gA_mean,   yerr = gA_std,   color = colors[3],fmt = fmts[0],linestyle = 'None',capsize = capsize,markersize=markersizes[0])
    axs[1].errorbar(W2,    Sigma_mean,yerr = Sigma_std,color = colors[4],fmt = fmts[1],linestyle = 'None',capsize = capsize,markersize=markersizes[1])
    #axs[1].errorbar(W2+dW2,a8_mean,   yerr = a8_std,   color = colors[2],fmt = fmts[2],linestyle = 'None',capsize = capsize,markersize=markersizes[2])
                                                                                                                                             
    #axs[2].errorbar(W2,    d2p_mean,yerr = d2p_std,    color = colors[0],fmt = fmts[0],linestyle = 'None',capsize = capsize,markersize=markersizes[0])
    #axs[2].errorbar(W2+dW2,d2n_mean,yerr = d2n_std,    color = colors[1],fmt = fmts[1],linestyle = 'None',capsize = capsize,markersize=markersizes[1])
   

    j = 1
    alpha = 0.03
    axs[0].fill_between(np.linspace(0,100,100),up_mean[j]-up_std[j],up_mean[j]+up_std[j],color=colors[0],alpha=alpha)
    axs[0].fill_between(np.linspace(0,100,100),dp_mean[j]-dp_std[j],dp_mean[j]+dp_std[j],color=colors[1],alpha=alpha)
    #axs[0].fill_between(np.linspace(0,100,100),sp_mean[j]-sp_std[j],sp_mean[j]+sp_std[j],color=colors[2],alpha=alpha)
    axs[0].fill_between(np.linspace(0,100,100),g_mean [j]-g_std [j],g_mean [j]+g_std [j],color=colors[2],alpha=alpha)

    axs[1].fill_between(np.linspace(0,100,100),gA_mean   [j]-gA_std   [j],gA_mean   [j]+gA_std   [j],color=colors[3],alpha=alpha)
    axs[1].fill_between(np.linspace(0,100,100),Sigma_mean[j]-Sigma_std[j],Sigma_mean[j]+Sigma_std[j],color=colors[4],alpha=alpha)
    #axs[1].fill_between(np.linspace(0,100,100),a8_mean   [j]-a8_std   [j],a8_mean   [j]+a8_std   [j],color=colors[2],alpha=alpha)

    #axs[2].fill_between(np.linspace(0,100,100),d2p_mean[j]-d2p_std[j],d2p_mean[j]+d2p_std[j],color=colors[0],alpha=alpha)
    #axs[2].fill_between(np.linspace(0,100,100),d2n_mean[j]-d2n_std[j],d2n_mean[j]+d2n_std[j],color=colors[1],alpha=alpha)

    for i in range(2):
        axs[i].tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=35,length=8)
        axs[i].tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=35,length=4)
        axs[i].set_xlim(3.3,10.5)
        axs[i].set_xticks([4,6,8,10])
        axs[i].xaxis.set_minor_locator(MultipleLocator(1))
        axs[i].axhline(0,0,1,color='black',alpha=0.3)

    axs[0].tick_params(labelbottom=False)

    axs[1].set_xticklabels([r'$4$',r'$6$',r'$8$',r'$10$'])
 
    axs[0].set_ylim(-0.6,0.9)
    axs[0].set_yticks([-0.4,0,0.4,0.8])
    axs[0].set_yticklabels([r'$-0.4$',r'$0$',r'$0.4$',r'$0.8$'])
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.1))
   
    axs[1].set_ylim(0.4,1.3)
    axs[1].set_yticks([0.5,1.0])
    axs[1].set_yticklabels([r'$0.5$',r'$1.0$'])
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.1))
 
    axs[1].set_xlabel(r'\boldmath{$W^2_{\rm{min}}$}~ \rm \bf [GeV\boldmath{$^2$}]',fontsize=35)

    x = 0.60
    size = 40
    axs[0].text(x     ,0.85, r'\boldmath$\Delta u^+$'   ,transform=axs[0].transAxes,size = size,color = colors[0])
    axs[0].text(x     ,0.65, r'\boldmath$\Delta g$'     ,transform=axs[0].transAxes,size = size,color = colors[2])
    #axs[0].text(x     ,0.38, r'\boldmath$\Delta s^+$'   ,transform=axs[0].transAxes,size = size,color = colors[2])
    axs[0].text(x     ,0.11, r'\boldmath$\Delta d^+$'   ,transform=axs[0].transAxes,size = size,color = colors[1])
   
    size = 40
    axs[1].text(0.45,0.78, r'\boldmath$\Delta u^+ - \Delta d^+$'          ,transform=axs[1].transAxes,size = size,color = colors[3])
    axs[1].text(0.45,0.35, r'\boldmath$\Delta u^+ + \Delta d^+$',transform=axs[1].transAxes,size = size,color = colors[4])
    #axs[1].text(x+0.03,0.19, r'\boldmath$a_8$'          ,transform=axs[1].transAxes,size = size,color = colors[2])

    size = 50
    #axs[2].text(x     ,0.55, r'\boldmath$d_2^p$'        ,transform=axs[2].transAxes,size = size,color = colors[0])
    #axs[2].text(x     ,0.20, r'\boldmath$d_2^n$'        ,transform=axs[2].transAxes,size = size,color = colors[1])

    #axs[2].text(0.60,0.80, r'\boldmath$\int_{0.005}^{0.52} {\rm d}x f$'        ,transform=axs[2].transAxes,size = 50)
    #axs[1].text(0.04,0.05, r'\textbf{\textrm{JAMpol25}}'      ,transform=axs[1].transAxes,size = 40)
    axs[1].text(0.04,0.05, r'$Q^2 = 10 ~ {\rm GeV}^2$'        ,transform=axs[1].transAxes,size = 25)

    #axs[1].set_ylabel(r'\boldmath${\rm Truncated Moments}~(0.005 < x < 0.52)$',size=40,labelpad=40)
    axs[0].set_title(r'\boldmath${\rm Truncated ~ Moments}~(0.005 < x < 0.52)$',size=32)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.01,wspace=0.01,top=0.95,right=0.99,left=0.13)
   
    axs[0].set_rasterized(True)
    axs[1].set_rasterized(True)
 
    filename = 'plots/fig_moments%s'%ext
    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

def print_moments():

    #--setup scenarios
    N    = 4
    xmin = [0.005,0.005,0.005,0.0]
    xmax = [0.53, 0.53, 0.76, 1.0]
    W2   = [10, 4, 4, 4]

    Q2 = 10
    mom = 1

    MOM = {}

    for i in range(N):

        if W2[i] == 3.5:  wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT3p5/pos_g'
        elif W2[i] == 10: wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noHT10/pos_g'
        else:             wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT%d/pos_g'%W2[i]
        #--ppdf moments
        flavs = ['up','dp','sp','g','Sigma','gA','a8']
        #--full moment
        if i==3: filename = '%s/data/ppdf-moment-%d-Q2=%3.5f.dat'%(wdir,mom,Q2)
        #--truncated moment
        else:    filename = '%s/data/ppdf-moment-trunc-%d-Q2=%3.5f-xmin=%3.5f-xmax=%3.5f.dat'%(wdir,mom,Q2,xmin[i],xmax[i])

        if regen:
            if i==3: ppdf.gen_moments      (wdir, Q2 = Q2, flavors = flavs, mom=mom)
            else:    ppdf.gen_moments_trunc(wdir, Q2 = Q2, flavors = flavs, mom=mom, xmin = xmin[i], xmax = xmax[i])
        try:
            moments = load(filename)
        except:
            if i==3: ppdf.gen_moments      (wdir, Q2 = Q2, flavors = flavs, mom=mom)
            else:    ppdf.gen_moments_trunc(wdir, Q2 = Q2, flavors = flavs, mom=mom, xmin = xmin[i], xmax = xmax[i])
            moments = load(filename)

        MOM[i] = {}
        for flav in flavs:
            if i==3: MOM[i][flav] = np.array(moments['moments'][flav])
            else:    MOM[i][flav] = np.array(moments['moments'][flav])[:,0]
        #--d2 moments
        #for tar in ['p','n']:
        #    #--full moment
        #    if i==3: filename = '%s/data/d2-%s.dat'%(wdir,tar)
        #    #--truncated moment
        #    else:    filename = '%s/data/d2-trunc-%s-Q2=%3.5f-xmin=%3.5f-xmax=%3.5f.dat'%(wdir,tar,Q2,xmin[i],xmax[i])
        #    if regen:
        #        if i==3: pstf.gen_d2(wdir, tar=tar)
        #        else:    pstf.gen_d2_trunc(wdir, tar=tar, Q2 = Q2, xmin = xmin[i], xmax = xmax[i])
        #    try:
        #        moments = load(filename)
        #    except:
        #        if i==3: pstf.gen_d2(wdir, tar=tar)
        #        else:    pstf.gen_d2_trunc(wdir, tar=tar, Q2 = Q2, xmin = xmin[i], xmax = xmax[i])
        #        moments = load(filename)

        #    if i==3:
        #        idx = np.asarray(moments['Q2']==10).nonzero()[0][0]
        #        d2 = np.array(moments['D2'])[:,idx]
        #    else:
        #        d2 = np.array(moments['moments'])[:,0]
        #    MOM[i]['d2'+tar] = d2

    #flavs = ['up','dp','sp','g','Sigma','gA','a8','d2p','d2n']
    flavs = ['up','dp','g','Sigma','gA']

    for i in range(N):
        print('SCENARIO %d: xmin = %3.5f, xmax = %3.5f, W2min = %d'%(i+1,xmin[i],xmax[i],W2[i]))
        for flav in flavs:
            if flav=='sigma': moment = MOM[i]['up'] + MOM[i]['dp']
            else: moment = MOM[i][flav]
            mean = np.mean(moment,axis=0) 
            std  = np.std (moment,axis=0)
            #print(mean,std)
            if flav in ['d2p','d2n']: std  = str(int(np.round(std*10000)))
            else:                     std  = str(int(np.round(std*100)))
            if flav in ['d2p','d2n']: moment = '%5.4f'%mean
            else:                     moment = '%5.2f'%mean
            if len(std)==1: moment +=' '
            moment += '('+std+')'
            print('%5s: %s'%(flav, moment))


if __name__=="__main__":

    plot_moments()
    print_moments()










