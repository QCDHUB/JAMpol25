#!/usr/bin/env python
import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

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

path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

#--from scipy stack 
from scipy.integrate import quad

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
# from analysis.corelib import classifier
from analysis.qpdlib import pdf,ppdf

import kmeanconf as kc

#regen = True
regen = False

#ext = '.png'
ext = '.pdf'

M2 = 0.938**2
Q2 = 10

W2min  = np.array([4,10,4])
xmax   = Q2/(Q2+W2min-M2)
colors = ['red','yellow','black']
alpha = [0.9,0.9,0.5]
zorders = [1.1,1.0,1.2]

def plot_ppdfs(WDIR):
 
    flavs = ['u','d','ub','db','ub-db','g','sigma']
    
    nrows = 4
    ncols = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*10,nrows*5))
 
    axs[3][1].axis("off")
 
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            if j ==0: axs[i][j].set_xlim(0.01,0.73)
            if j ==1: axs[i][j].set_xlim(0.01,0.60)
            axs [i][j].axhline(0,lw = 1,color = 'k',alpha = 0.3,zorder=10)

    k = 0
    hand = {}

    for wdir in WDIR:
        load_config(wdir + '/input.py')
        
        filename = wdir + '/data/ppdf-Q2=%3.5f.dat'%(Q2)
        if regen: ppdf.gen_xf(wdir,Q2=Q2) 
        try: data=load(filename)
        except:
            ppdf.gen_xf(wdir,Q2=Q2)
            data=load(filename)

        X = data['X']
        for flav in flavs:

            if   flav=='u':     i, j = 0, 0
            elif flav=='d':     i, j = 1, 0
            elif flav=='sigma': i, j = 2, 0
            elif flav=='g':     i, j = 3, 0
            elif flav=='ub':    i, j = 0, 1
            elif flav=='db':    i, j = 1, 1
            elif flav=='ub-db': i, j = 2, 1

            if flav=='ub-db':
                result = np.array(data['XF']['ub']) - np.array(data['XF']['db'])
            elif flav=='sigma':
                result = np.array(data['XF']['up']) + np.array(data['XF']['dp'])# + np.array(data_pos_1['XF']['sp'])
            else:
                result = np.array(data['XF'][flav])

            #--mask the result
            result[:,X>xmax[k]] = 'nan'

            mean = np.mean(result,axis=0)
            std  = np.std (result,axis=0)

            if k==2:
                hand[k] = axs [i][j].fill_between(X,mean-std,mean+std,color='none',edgecolor=colors[k],alpha=alpha[k],zorder=zorders[k],hatch='//')
            else:
                hand[k] = axs [i][j].fill_between(X,mean-std,mean+std,color=colors[k],alpha=alpha[k],zorder=zorders[k])

        k += 1       

 
    axs[0][0].set_ylim(0,0.4)
    axs[0][0].yaxis.set_minor_locator(MultipleLocator(0.05)) 
    axs[0][0].yaxis.set_major_locator(MultipleLocator(0.10)) 
    axs[0][0].set_yticks([0.1,0.2,0.3]) 
    axs[0][0].set_yticklabels([r'$0.1$',r'$0.2$',r'$0.3$']) 
 
    axs[0][1].set_ylim(-0.009,0.030)
    axs[0][1].yaxis.set_minor_locator(MultipleLocator(0.005)) 
    axs[0][1].set_yticks([0,0.01,0.02]) 
    axs[0][1].set_yticklabels([r'$0$',r'$0.01$',r'$0.02$']) 

    axs[1][0].set_ylim(-0.14,0.02)
    axs[1][0].yaxis.set_minor_locator(MultipleLocator(0.01)) 
    axs[1][0].set_yticks([-0.10,-0.05,0]) 
    axs[1][0].set_yticklabels([r'$-0.10$',r'$-0.05$',r'$0$']) 
        
    axs[1][1].set_ylim(-0.039,0.009)
    axs[1][1].yaxis.set_minor_locator(MultipleLocator(0.01)) 
    axs[1][1].set_yticks([-0.03,-0.02,-0.01,0]) 
    axs[1][1].set_yticklabels([r'$-0.03$',r'$-0.02$',r'$-0.01$',r'$0$']) 

    axs[2][0].set_ylim(-0.05,0.39)
    axs[2][0].yaxis.set_minor_locator(MultipleLocator(0.05)) 
    axs[2][0].set_yticks([0,0.1,0.2,0.3]) 
    axs[2][0].set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$']) 

    axs[2][1].set_ylim(-0.01,0.06)
    axs[2][1].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[2][1].set_yticks([0,0.02,0.04]) 
    axs[2][1].set_yticklabels([r'$0$',r'$0.02$',r'$0.04$']) 

    axs[3][0].set_ylim(-0.02,0.19)
    axs[3][0].yaxis.set_minor_locator(MultipleLocator(0.01)) 
    axs[3][0].set_yticks([0,0.05,0.10,0.15]) 
    axs[3][0].set_yticklabels([r'$0$',r'$0.05$',r'$0.10$',r'$0.15$']) 
 
    ls = 40
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            axs [i][j].tick_params(axis='both', which='major', top=True, direction='in',labelsize=ls,length=10)
            axs [i][j].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=ls,length=5)
            if j==0: axs[i][j].set_xticks([0.2,0.4,0.6])
            if j==1: axs[i][j].set_xticks([0.2,0.4])
            if i < 3 and (i!=3 and j!=1): axs[i][j].set_xticklabels([])
            else: 
                if j==0: axs[i][j].set_xticklabels([r'$0.2$',r'$0.4$',r'$0.6$'])
                if j==1: axs[i][j].set_xticklabels([r'$0.2$',r'$0.4$'])
    
            if i == 3 or (i==2 and j==1): axs[i][j].set_xlabel(r'\boldmath{$x$}',fontsize = 60)
    
            axs[i][j].xaxis.set_minor_locator(MultipleLocator(0.10)) 
       
    pdf_fs = 50
 
    #axs[0][1].text( 0.20, 0.85, r'\textrm{\textbf{JAMpol25}}',          transform=axs[0][1].transAxes,fontsize=pdf_fs)
    
    axs[0][0].text( 0.77, 0.85, r'\boldmath$x \Delta u$',          transform=axs[0][0].transAxes,fontsize=pdf_fs)
    axs[1][0].text( 0.77, 0.05, r'\boldmath$x \Delta d$',          transform=axs[1][0].transAxes,fontsize=pdf_fs)
    axs[0][1].text( 0.77, 0.85, r'\boldmath$x \Delta \bar{u}$',    transform=axs[0][1].transAxes,fontsize=pdf_fs)
    axs[1][1].text( 0.77, 0.05, r'\boldmath$x \Delta \bar{d}$',    transform=axs[1][1].transAxes,fontsize=pdf_fs)
    #axs[2][0].text( 0.77, 0.85, r'\boldmath$x \Delta \Sigma$',     transform=axs[2][0].transAxes,fontsize=pdf_fs)
    axs[2][0].text( 0.37, 0.82, r'\boldmath$x (\Delta u^+ + \Delta d^+)$',     transform=axs[2][0].transAxes,fontsize=pdf_fs)
    axs[2][1].text( 0.47, 0.82, r'\boldmath$x (\Delta \bar{u} - \Delta \bar{d})$',    transform=axs[2][1].transAxes,fontsize=pdf_fs)
    axs[3][0].text( 0.77, 0.85, r'\boldmath$x \Delta g$',          transform=axs[3][0].transAxes,fontsize=pdf_fs)

    axs[0][0].text(0.03, 0.05, r'$Q^2 = %d ~ {\rm GeV}^2$'%Q2,     transform=axs[0][0].transAxes,fontsize=35)
  

    handles, labels = [],[]
    handles.append(hand[0])
    handles.append(hand[2])
    handles.append(hand[1])
    labels.append(r'\boldmath$W^2 > 4 ~{\rm GeV}^2$')
    #labels.append(r'\boldmath$W^2 > 4 ~{\rm GeV}^2 ~ {\rm (no ~ HT)}$')
    labels.append(r'\boldmath$W^2 > 4 ~{\rm GeV}^2 ~ {\rm (LT ~ only)}$')
    labels.append(r'\boldmath$W^2 > 10~{\rm GeV}^2$')
    axs[3][1].legend(handles,labels,fontsize = 45, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.02,wspace=0.20,top=0.99,right=0.99)

    filename = 'plots/fig_ppdfs_W2' + ext 

    axs[0][0].set_rasterized(True)
    axs[0][1].set_rasterized(True)
    axs[1][0].set_rasterized(True)
    axs[1][1].set_rasterized(True)
    axs[2][0].set_rasterized(True)
    axs[2][1].set_rasterized(True)
    axs[3][0].set_rasterized(True)

    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

def plot_polarization(WDIR):
   
 
    ncols = 2
    nrows = 1
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*7,nrows*5))



    k = 0
    flavs = ['u','d']
    hand = {}

    for wdir in WDIR:
        load_config(wdir + '/input.py')
        
        filename = wdir + '/data/ppdf-Q2=%3.5f.dat'%(Q2)
        if regen: ppdf.gen_xf(wdir,Q2=Q2) 
        try: data=load(filename)
        except:
            ppdf.gen_xf(wdir,Q2=Q2)
            data=load(filename)

        filename = wdir + '/data/pdf-Q2=%3.5f.dat'%(Q2)
        if regen: pdf.gen_xf(wdir,Q2=Q2) 
        try: udata=load(filename)
        except:
            pdf.gen_xf(wdir,Q2=Q2)
            udata=load(filename)

        X = data['X']
        for flav in flavs:

            if   flav=='u':     i = 0
            elif flav=='d':     i = 1

            if flav=='u':
                result = np.array(data['XF']['u'])/np.array(udata['XF']['u'])
            elif flav=='d':
                result = np.array(data['XF']['d'])/np.array(udata['XF']['d'])

            result[:,X>xmax[k]] = 'nan'

            mean = np.mean(result,axis=0)
            std  = np.std (result,axis=0)

            if k==2:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color='none',edgecolor=colors[k],alpha=alpha[k],zorder=zorders[k],hatch='//')
            else:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color=colors[k],alpha=alpha[k],zorder=zorders[k])

        k+=1 
      
    ls = 25 

    for i in range(len(axs)):
        axs[i].set_xlim(0.01,0.73)
        axs [i].axhline(0,lw = 1,color = 'k',alpha = 0.3)
        axs [i].tick_params(axis='both', which='major', top=False, direction='in',labelsize=ls,length=10)
        axs [i].tick_params(axis='both', which='minor', top=False, direction='in',labelsize=ls,length=5)
        axs [i].set_xticks([0.2,0.4,0.6])
        axs[i].set_xlabel(r'\boldmath{$x$}',fontsize = 35)
        axs[i].xaxis.set_minor_locator(MultipleLocator(0.10)) 

    axs[0].set_ylim(0.001,0.99)
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.10)) 
    axs[0].set_yticks([0.2,0.4,0.6,0.8,1.0]) 
    axs[0].set_yticklabels([r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1.0$']) 
   
    axs[1].set_ylim(-1.10,1.00)
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.10)) 
    axs[1].set_yticks([-1.0,-0.5,0,0.5,1.0]) 
    axs[1].set_yticklabels([r'$-1.0$',r'$-0.5$',r'$0$',r'$0.5$',r'1.0']) 
  

    #axs[0].text(0.05, 0.05, r'\textrm{\textbf{JAMpol25}}',    transform=axs[0].transAxes,fontsize=40)
    axs[0].text(0.03, 0.80, r'\boldmath$\frac{\Delta u}{u}$', transform=axs[0].transAxes,fontsize=50)
    axs[1].text(0.03, 0.10, r'\boldmath$\frac{\Delta d}{d}$', transform=axs[1].transAxes,fontsize=50)
    
    axs[0].text(0.70, 0.05, r'$Q^2 = %d ~ {\rm GeV}^2$'%Q2,     transform=axs[0].transAxes,fontsize=20)


    handles, labels = [],[]
    handles.append(hand[0])
    handles.append(hand[2])
    handles.append(hand[1])
    labels.append(r'\boldmath$W^2 > 4 ~{\rm GeV}^2$')
    labels.append(r'\boldmath$W^2 > 4 ~{\rm GeV}^2 ~ {\rm (LT ~ only)}$')
    labels.append(r'\boldmath$W^2 > 10~{\rm GeV}^2$')
    axs[1].legend(handles,labels,fontsize = 25, loc=(0.00,0.55), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.01,wspace=0.15,top=0.97,right=0.99)

    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
 
    filename = 'plots/fig_polarization' + ext 
    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()
    

if __name__=="__main__":

    #wdir1 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/simul/pos_g'
    #wdir1 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noboot_HT4/pos_g'
    wdir1 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g'
    wdir2 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noHT10/pos_g'
    #wdir3 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noHT4/pos_g'
    wdir3 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/LT4/pos_g'


    WDIR = [wdir1,wdir2,wdir3]
    plot_ppdfs       (WDIR)
    plot_polarization(WDIR)

    

