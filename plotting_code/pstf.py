#!/usr/bin/env python
import sys,os
import numpy as np
import time
import argparse
import copy

#--matplotlib
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
from matplotlib.ticker import MultipleLocator
import pylab as py

from scipy.integrate import cumulative_trapezoid

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
# from analysis.corelib import classifier
from analysis.obslib import pstf

#--force data to be regenerated
#regen = True
regen = False

ext = '.png'
#ext = '.pdf'

path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

#msr_dir      = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g/'
#msr_dir      = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/jlab12/pos_g/'
msr_dir      = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/jlab12_doublestat/'

msr_dir_noHT = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noHT4/pos_g/'
msr_dir_LT   = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/LT4/pos_g/'
W2 = 4

def plot_pstfs(Q2=4):

    pdf_file_name = 'plots/fig_pstfs' + ext
    load_config('%s/input.py' % msr_dir)
    istep = core.get_istep()

    if regen:
        pstf.gen_pstf(msr_dir     ,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir     ,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir     ,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir     ,Q2=Q2,tar='n',stf='g2')
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='n',stf='g2')
        pstf.gen_pstf(msr_dir_LT  ,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir_LT  ,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir_LT  ,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir_LT  ,Q2=Q2,tar='n',stf='g2')
    #--load full results
    try:
        sf_data_p_g1    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
  
    #--load LT with TMC results 
    try:
        sf_data_p_g1_LT    = load(msr_dir_noHT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1_LT    = load(msr_dir_noHT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1_LT    = load(msr_dir_noHT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1_LT    = load(msr_dir_noHT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2_LT    = load(msr_dir_noHT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2_LT    = load(msr_dir_noHT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2_LT    = load(msr_dir_noHT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_noHT,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2_LT    = load(msr_dir_noHT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)

    #--load LT with no TMC results 
    try:
        sf_data_p_g1_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2_LT_noTMC    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    
    
    data = {}
    data['g1p'] = sf_data_p_g1
    data['g1n'] = sf_data_n_g1
    data['g2p'] = sf_data_p_g2
    data['g2n'] = sf_data_n_g2

    data_LT = {}
    data_LT['g1p'] = sf_data_p_g1_LT
    data_LT['g1n'] = sf_data_n_g1_LT
    data_LT['g2p'] = sf_data_p_g2_LT
    data_LT['g2n'] = sf_data_n_g2_LT

    data_LT_noTMC = {}
    data_LT_noTMC['g1p'] = sf_data_p_g1_LT_noTMC
    data_LT_noTMC['g1n'] = sf_data_n_g1_LT_noTMC
    data_LT_noTMC['g2p'] = sf_data_p_g2_LT_noTMC
    data_LT_noTMC['g2n'] = sf_data_n_g2_LT_noTMC
    
    ncols = 1
    nrows = 3
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*8,nrows*5))
    
    axsL = {}
    for i in range(3):
        divider = make_axes_locatable(axs[i])
        axsL[i] = divider.append_axes("right", size=3.50, pad = 0.005,sharey=axs[i])
        axs[i].spines['right'].set_visible(False)
        axsL[i].spines['left'].set_visible(False)
        axsL[i].yaxis.set_ticks_position('right')
    
    
    i = 1
    hand = {} 
    for stf in ['g1','g2','gT']:
        for tar in ['p', 'n']:
            if stf=='g1': j = 0
            if stf=='g2': j = 1
            if stf=='gT': j = 2

            if tar=='p': zorder, color = 2,'red'
            if tar=='n': zorder, color = 1,'dodgerblue'
            if stf=='gT': X = data['g1'+tar]['X']
            else:         X = data[stf+tar]['X']
            if stf=='gT':
                mean = X*np.mean(data['g1'+tar]['STF']+data['g2'+tar]['STF'],axis=0)
                std  = X*np.std (data['g1'+tar]['STF']+data['g2'+tar]['STF'],axis=0)
            else:
                mean = X*np.mean(data[stf+tar]['STF'],axis=0)
                std  = X*np.std (data[stf+tar]['STF'],axis=0)
            axs[j] .fill_between(X,mean - std,mean + std, color = color,alpha = 0.9,zorder=zorder)
            axsL[j].fill_between(X,mean - std,mean + std, color = color,alpha = 0.9,zorder=zorder)

            #--plot leading twist with AOT
            #if stf=='gT': X = data_LT['g1'+tar]['X']
            #else:         X = data_LT[stf+tar]['X']
            #if stf=='gT':
            #    mean = X*np.mean(data_LT['g1'+tar]['STF']+data_LT['g2'+tar]['STF'],axis=0)
            #    std  = X*np.std (data_LT['g1'+tar]['STF']+data_LT['g2'+tar]['STF'],axis=0)
            #else:
            #    mean = X*np.mean(data_LT[stf+tar]['STF'],axis=0)
            #    std  = X*np.std (data_LT[stf+tar]['STF'],axis=0)


            #hand['noHT'] = axs[j] .fill_between(X,mean - std,mean + std, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//',lw=2)
            #hand['noHT'] = axsL[j].fill_between(X,mean - std,mean + std, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//',lw=2)

            #--plot leading twist with no TMCs
            if tar=='p': zorder,color = 3,'firebrick'
            if tar=='n': zorder,color = 3,'darkblue'
            if stf=='gT': X = data_LT_noTMC['g1'+tar]['X']
            else:         X = data_LT_noTMC[stf+tar]['X']
            if stf=='gT':
                mean = X*np.mean(data_LT_noTMC['g1'+tar]['STF']+data_LT_noTMC['g2'+tar]['STF'],axis=0)
                std  = X*np.std (data_LT_noTMC['g1'+tar]['STF']+data_LT_noTMC['g2'+tar]['STF'],axis=0)
            else:
                mean = X*np.mean(data_LT_noTMC[stf+tar]['STF'],axis=0)
                std  = X*np.std (data_LT_noTMC[stf+tar]['STF'],axis=0)
            hand['LT'] = axs[j] .fill_between(X,mean - std,mean + std, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//')
            hand['LT'] = axsL[j].fill_between(X,mean - std,mean + std, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//')
    
    for i in range(3):
        axs[i].set_xlim(1e-2,0.1)
        axs[i].semilogx()
    
        axs[i].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
        axs[i].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=0,length=4)
        axs[i].set_xticks([0.01,0.1])
        axs[i].set_xticklabels([r'$0.01$',r'$0.1$'])
        axs[i].axhline(0,lw = 1,color = 'k',alpha = 0.3)
        
        
        axsL[i].tick_params(axis='both', which='major', top=True, right = True, left = False, labelright=False, direction='in',labelsize=30,length=8)
        axsL[i].tick_params(axis='both', which='minor', top=True, right = True, left = False, labelright=False, direction='in',labelsize=30,length=4)
        axsL[i].set_xlim(0.1,0.73)
        axsL[i].set_xticks([0.1,0.3,0.5,0.7])
        axsL[i].set_xticklabels([r'',r'$0.3$',r'$0.5$',r'$0.7$'])
        axsL[i].set_xlabel(r'\boldmath{$x$}',fontsize = 45,loc = 'left')
        axsL[i].axhline(0,lw = 1,color = 'k',alpha = 0.3)

        axsL[i].xaxis.set_minor_locator(MultipleLocator(0.10))

    for i in [0,1]:
        axs [i].tick_params(labelbottom=False)
        axsL[i].tick_params(labelbottom=False)
    
    axs[0].set_ylim(-0.019,0.065)
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[0].set_yticks([0,0.02,0.04,0.06])
    axs[0].set_yticklabels([r'$0$',r'$0.02$',r'$0.04$',r'$0.06$'])

    axs[1].set_ylim(-0.049,0.025)
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[1].set_yticks([-0.04,-0.02,0,0.02])
    axs[1].set_yticklabels([r'$-0.04$',r'$-0.02$',r'$0$',r'$0.02$'])
    
    axs[2].set_ylim(-0.025,0.045)
    axs[2].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[2].set_yticks([-0.02,0,0.02,0.04])
    axs[2].set_yticklabels([r'$-0.02$',r'$0$',r'$0.02$',r'$0.04$'])

    
    # ##########################################################
        
    axs[1].text(0.03,0.30,r'\rm \bf proton', transform=axs[1].transAxes,size = 40, color = 'red')
    axs[1].text(0.03,0.20,r'\rm \bf neutron',transform=axs[1].transAxes,size = 40, color = 'dodgerblue')
    
    #axs[0] .text(0.08, 0.87, r'\textrm{\textbf{JAMpol25}}', transform=axs [0].transAxes,fontsize=40,color = 'k')
    
    axs[0] .text(0.08, 0.87, r'\boldmath$x g_1$', transform=axs [0].transAxes,fontsize=50,color = 'k')
    axsL[1].text(0.60, 0.87, r'\boldmath$x g_2$', transform=axsL[1].transAxes,fontsize=50,color = 'k')
    axs[2] .text(0.08, 0.87, r'\boldmath$x g_T$', transform=axs [2].transAxes,fontsize=50,color = 'k')
   
    axs[0].text(0.08,0.75,r'$Q^2=%d~{\rm GeV}^2$'%Q2,transform=axs[0].transAxes,size=25)


 
    #hand['LT'] = axs[j] .fill_between(X,-100,-99, color = 'None', edgecolor='black',alpha = 1.0,zorder=3,hatch='//',lw=2)
    handles, labels = [],[]
    #handles.append(hand['noHT'])
    handles.append(hand['LT'])
    #labels.append(r'\textrm{\textbf{no HT}}')
    labels.append(r'\textrm{\textbf{LT only}}')
    axs[1].legend(handles,labels,fontsize = 30, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)
 
    py.tight_layout()
    py.subplots_adjust(top=0.99,right=0.99,hspace=0.01,wspace=0.01)
    
    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
    axs [2].set_rasterized(True)
    axsL[0].set_rasterized(True)
    axsL[1].set_rasterized(True)
    axsL[2].set_rasterized(True)
    
    py.savefig(pdf_file_name)
    py.close()
    print('Saving figure to %s'%pdf_file_name)

def plot_gXres():

    wdir = msr_dir
    Q2 = 10
    mode = 1

    nrows,ncols=1,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    TAR = ['p','n']
    hand = {}
    for stf in ['g1','g2']: 
        for tar in TAR:
            filename ='%s/data/gXres-%s-%s-Q2=%3.5f.dat'%(wdir,stf,tar,Q2)
            if regen: pstf.gen_gXres(wdir,Q2,tar,stf)
            try:
                data=load(filename)
            except:
                pstf.gen_gXres(wdir,Q2,tar,stf)
                data=load(filename)

            X    = data['X']
            data = data['STF'] 
            mean = X*np.mean(data,axis=0)
            std  = X*np.std (data,axis=0)

            if tar=='p' and stf =='g1': ax,color=ax11,'red'
            if tar=='n' and stf =='g1': ax,color=ax11,'b'
            if tar=='p' and stf =='g2': ax,color=ax12,'red'
            if tar=='n' and stf =='g2': ax,color=ax12,'b'

            label = None

            #--plot each replica
            if mode==0:
                for i in range(len(data)):
                    hand[tar] ,= ax .plot(X,data[i],color=color,alpha=0.1)
          
            #--plot average and standard deviation
            if mode==1:
                hand[tar] = ax .fill_between(X,(mean-std),(mean+std),color=color,alpha=0.7)


    for ax in [ax11,ax12]:
          ax.set_xlim(0.05,0.73)

          ax.set_xticks([0.2,0.4,0.6])
          ax.tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=25,length=8)
          ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=25,length=4)
          minorLocator = MultipleLocator(0.1)
          majorLocator = MultipleLocator(0.2)
          ax.xaxis.set_minor_locator(minorLocator)
          ax.xaxis.set_major_locator(majorLocator)
          ax.set_xlabel(r'\boldmath$x$' ,size=45)
    
          ax.axhline(0,lw = 1,color = 'k',alpha = 0.3)
    
          ax.set_ylim(-0.017,0.009)
          ax.yaxis.set_major_locator(MultipleLocator(0.005))
          ax.yaxis.set_minor_locator(MultipleLocator(0.001))
          

    ax12.tick_params(labelleft=False)

    ax11.text(0.03,0.06,r'\boldmath$xg_1^{\rm HT}$',transform=ax11.transAxes,size=45)
    ax12.text(0.03,0.06,r'\boldmath$xg_2^{\rm HT}$',transform=ax12.transAxes,size=45)

    ax11.text(0.60,0.17,r'\rm \bf proton' ,transform=ax11.transAxes,size = 45,color='red')
    ax12.text(0.60,0.05,r'\rm \bf neutron',transform=ax11.transAxes,size = 45,color='blue')

    ax12.text(0.70,0.90,r'$Q^2=10~{\rm GeV}^2$',transform=ax12.transAxes,size=20)

    py.tight_layout()
    py.subplots_adjust(hspace=0,wspace=0.01,top=0.99,right=0.99)

    filename = 'plots/fig_gXres'
    filename+=ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)

    py.savefig(filename)
    print ('Saving figure to %s'%filename)
    py.clf()

def plot_gXres_W2():

    wdir = msr_dir

    nrows,ncols=1,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5.5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()

    TAR = ['p','n']
    hand = {}
    for stf in ['g1','g2']: 
        for tar in TAR:
            W2 = 4
            filename ='%s/data/gXres-%s-%s-W2=%3.5f.dat'%(wdir,stf,tar,W2)
            if regen: pstf.gen_gXres_W2(wdir,W2,tar,stf)
            try:
                data=load(filename)
            except:
                pstf.gen_gXres_W2(wdir,W2,tar,stf)
                data=load(filename)

            X    = data['X']
            data = data['STF'] 
            mean = X*np.mean(data,axis=0)
            std  = X*np.std (data,axis=0)

            if tar=='p' and stf =='g1': ax,zorder,color=ax11,2,'red'
            if tar=='n' and stf =='g1': ax,zorder,color=ax11,1,'b'
            if tar=='p' and stf =='g2': ax,zorder,color=ax12,2,'red'
            if tar=='n' and stf =='g2': ax,zorder,color=ax12,1,'b'

            hand['4'] = ax .fill_between(X,(mean-std),(mean+std),color=color,alpha=0.90,zorder=zorder)

            #--also plot W2 = 10 for g1
            if stf=='g1':
                W2 = 10
                filename ='%s/data/gXres-%s-%s-W2=%3.5f.dat'%(wdir,stf,tar,W2)
                if regen: pstf.gen_gXres_W2(wdir,W2,tar,stf)
                try:
                    data=load(filename)
                except:
                    pstf.gen_gXres_W2(wdir,W2,tar,stf)
                    data=load(filename)

                X    = data['X']
                data = data['STF'] 
                mean = X*np.mean(data,axis=0)
                std  = X*np.std (data,axis=0)

                #hand['10'] = ax .fill_between(X,(mean-std),(mean+std),facecolor='None',edgecolor=color,hatch='//',alpha=1.0,zorder=2)

    xmin = 0.35
    xmax = 0.90

    Q2min = 1.3**2
    Q2max = xmax/(1-xmax)*(W2-0.938**2)

    for ax in [ax11,ax12]:
          ax.set_xlim(xmin,xmax)

          ax.set_xticks([0.4,0.6,0.8])
          ax.tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=30,length=8)
          ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=30,length=4)
          minorLocator = MultipleLocator(0.1)
          majorLocator = MultipleLocator(0.2)
          ax.xaxis.set_minor_locator(minorLocator)
          ax.xaxis.set_major_locator(majorLocator)
          ax.set_xlabel(r'\boldmath$x$' ,size=45)
    
          ax.axhline(0,lw = 1,color = 'k',alpha = 0.3)
    
          ax.set_ylim(-0.019,0.012)
          ax.yaxis.set_minor_locator(MultipleLocator(0.005))
          ax.set_yticks([-0.01,0,0.01])
          ax.set_yticklabels([r'$-0.01$',r'$0$',r'$0.01$'])
          

    ax12.tick_params(labelleft=False)

    ax11.text(0.05,0.08,r'\boldmath$xg_1^{\rm HT}$',transform=ax11.transAxes,size=45)
    ax12.text(0.05,0.08,r'\boldmath$xg_2^{\rm HT}$',transform=ax12.transAxes,size=45)

    ax11.text(0.60,0.17,r'\rm \bf proton' ,transform=ax11.transAxes,size = 45,color='red')
    ax11.text(0.60,0.05,r'\rm \bf neutron',transform=ax11.transAxes,size = 45,color='blue')

    ax11.text(0.03,0.88,r'\boldmath$W^2 = 4 ~ {\rm GeV}^2$',transform=ax11.transAxes,size=25)
    ax12.text(0.45,0.05,r'$%3.2f < Q^2 < %d ~ {\rm GeV}^2$'%(Q2min,Q2max),transform=ax12.transAxes,size=25)

    #ax12.text(0.70,0.05,r'$W^2=4~{\rm GeV}^2$',transform=ax12.transAxes,size=20)
    #ax12.text(0.70,0.90,r'$Q^2=10~{\rm GeV}^2$',transform=ax12.transAxes,size=20)

    #handles, labels = [],[]
    #handles.append(hand['4'])
    #handles.append(hand['10'])
    #labels.append(r'$W^2=4 ~{\rm GeV}^2$')
    #labels.append(r'$W^2=10~{\rm GeV}^2$')
    #ax11.legend(handles,labels,fontsize = 22, loc=(0.60,0.00), frameon=0,handlelength=1,handletextpad=0.6)

    py.tight_layout()
    py.subplots_adjust(hspace=0,wspace=0.01,top=0.99,right=0.99)

    filename = 'plots/fig_gXres_W2'
    filename+=ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)

    py.savefig(filename)
    print ('Saving figure to %s'%filename)
    py.clf()

def plot_d2():

  nrows,ncols=1,2
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*7,nrows*5.5))
  ax11 = py.subplot(nrows,ncols,1)
  ax12 = py.subplot(nrows,ncols,2)

  load_config('%s/input.py'%msr_dir)
  istep=core.get_istep()


  TAR = ['p','n']
  hand = {}
  for tar in TAR:
      filename ='%s/data/d2-%s.dat'%(msr_dir,tar)
      if regen: pstf.gen_d2(msr_dir,tar=tar)
      #--load data if it exists
      try:
          data=load(filename)
      #--generate data and then load it if it does not exist
      except:
          pstf.gen_d2(msr_dir,tar=tar)
          data=load(filename)

      Q2 = data['Q2']
      data = data['D2']
      mean = np.mean(data,axis=0)
      std  = np.std (data,axis=0)

      if tar=='p': ax,alpha,color,zorder=ax11,0.9,'red' , 1
      if tar=='n': ax,alpha,color,zorder=ax12,0.9,'blue', 1

      label = None

      hand[tar] = ax.fill_between(Q2,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)
      hand[tar] = ax.fill_between(Q2,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)

      #--no HT result
      filename ='%s/data/d2-%s-LT.dat'%(msr_dir,tar)
      #--load data if it exists
      try:
          data=load(filename)
      #--generate data and then load it if it does not exist
      except:
          pstf.gen_d2(msr_dir,tar=tar,LT_only=True)
          data=load(filename)

      Q2 = data['Q2']
      data = data['D2']
      mean = np.mean(data,axis=0)
      std  = np.std (data,axis=0)

      if tar=='p': ax,alpha,color,zorder=ax11,1.0,'firebrick', 3
      if tar=='n': ax,alpha,color,zorder=ax12,1.0,'darkblue', 3

      ax .fill_between(Q2,(mean-std),(mean+std),facecolor="None",edgecolor=color,hatch='//',alpha=alpha,zorder=zorder,lw=2)

  hand['LT'] = ax .fill_between(Q2,-100,-100, facecolor="None",edgecolor='black',hatch='//',alpha=alpha,zorder=zorder,lw=2)

  #--load JAM15 data 
  F = open('data/d2_JAM15.csv','r')
  L = F.readlines()
  F.close()

  L = [l.strip() for l in L]
  L = [[x for x in l.split()] for l in L]
  L = np.transpose(L)[0]

  Q2, pmax,pmin,pmean,nmax,nmin,nmean = [],[],[],[],[],[],[]
  for i in range(len(L)):
      if i==0: continue
      Q2.append(float(L[i].split(',')[0]))
      pmax .append(float(L[i].split(',')[1]))
      pmin .append(float(L[i].split(',')[2]))
      pmean.append(float(L[i].split(',')[3]))
      nmax .append(float(L[i].split(',')[4]))
      nmin .append(float(L[i].split(',')[5]))
      nmean.append(float(L[i].split(',')[6]))
  Q2 = np.array(Q2)
  pmax  = np.array(pmax)
  pmin  = np.array(pmin)
  pmean = np.array(pmean)
  nmax  = np.array(nmax)
  nmin  = np.array(nmin)
  nmean = np.array(nmean)
 
  alpha = 0.6
  #ax11.fill_between(Q2,pmin,pmax,facecolor="None",hatch='//',edgecolor='black', alpha=alpha,zorder=3)
  #ax12.fill_between(Q2,nmin,nmax,facecolor="None",hatch='//',edgecolor='black',alpha=alpha,zorder=3)
  hand['JAM15'] = ax11.fill_between(Q2,-100,-100,facecolor="None",hatch='//',edgecolor='black',alpha=alpha,zorder=3)

  for ax in [ax11,ax12]:
      ax.set_xlim(1.7,5.5)
      ax.set_xticks([2,3,4,5])
      ax.xaxis.set_minor_locator(MultipleLocator(0.5))

      ax.tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=30,length=8)
      ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=30,length=4)
      ax.set_xlabel(r'\boldmath$Q^2~[{\rm GeV}^2]$' ,size=30)

      ax.axhline(0,lw = 1,color = 'k',alpha = 0.5,zorder=6)
      ax.set_ylim(-0.029,0.029)
      ax.yaxis.set_minor_locator(MultipleLocator(0.005))
      ax.set_yticks([-0.02,-0.01,0,0.01,0.02])
      ax.set_yticklabels([r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$'])

    

  ax11.text(0.03,0.83,r'\boldmath$d_2$',transform=ax11.transAxes,size=60)

  ax11.text(0.60,0.85,r'\rm \bf proton', transform=ax11.transAxes,size = 45, color = 'red')
  ax11.text(0.60,0.73,r'\rm \bf neutron',transform=ax11.transAxes,size = 45, color = 'blue')

  ax12.tick_params(labelleft=False)

  handles, labels = [],[]
  handles.append(hand['JAM15'])
  labels.append(r'\textrm{\textbf{JAM15}}')
  #ax11.legend(handles,labels,fontsize = 22, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)

  handles, labels = [],[]
  handles.append(hand['LT'])
  labels.append(r'\textrm{\textbf{no HT}}')
  ax11.legend(handles,labels,fontsize = 25, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)

  py.tight_layout()
  py.subplots_adjust(hspace=0,wspace=0.01,top=0.99,right=0.99)

  filename = 'plots/fig_d2'

  filename+=ext

  ax11.set_rasterized(True)

  py.savefig(filename)
  print ('Saving figure to %s'%filename)
  py.clf()

def plot_d2_trunc(xmin=0.005,xmax=0.76):

  nrows,ncols=1,2
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*7,nrows*5.5))
  ax11 = py.subplot(nrows,ncols,1)
  ax12 = py.subplot(nrows,ncols,2)

  load_config('%s/input.py'%msr_dir)
  istep=core.get_istep()


  TAR = ['p','n']
  hand = {}
  for tar in TAR:
      filename ='%s/data/d2-trunc-funcQ2-%s-xmin=%3.5f-xmax=%3.5f.dat'%(msr_dir,tar,xmin,xmax)
      if regen: pstf.gen_d2_trunc_funcQ2(msr_dir,tar,xmin,xmax)
      #--load data if it exists
      try:
          data=load(filename)
      #--generate data and then load it if it does not exist
      except:
          pstf.gen_d2_trunc_funcQ2(msr_dir,tar,xmin,xmax)
          data=load(filename)

      Q2 = data['Q2']
      data = data['D2']
      mean = np.mean(data,axis=0)
      std  = np.std (data,axis=0)

      if tar=='p': ax,alpha,color,zorder=ax11,0.7,'red' , 2
      if tar=='n': ax,alpha,color,zorder=ax12,0.5,'blue', 1

      label = None

      hand[tar] = ax.fill_between(Q2,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)
      hand[tar] = ax.fill_between(Q2,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)

      #--no HT result
      filename ='%s/data/d2-trunc-funcQ2-%s-xmin=%3.5f-xmax=%3.5f-LT.dat'%(msr_dir,tar,xmin,xmax)
      if regen: pstf.gen_d2_trunc_funcQ2(msr_dir,tar,xmin,xmax,LT_only=True)
      #--load data if it exists
      try:
          data=load(filename)
      #--generate data and then load it if it does not exist
      except:
          pstf.gen_d2_trunc_funcQ2(msr_dir,tar,xmin,xmax,LT_only=True)
          data=load(filename)

      Q2 = data['Q2']
      data = data['D2']
      mean = np.mean(data,axis=0)
      std  = np.std (data,axis=0)

      if tar=='p': ax,alpha,color,zorder=ax11,0.8,'black', 4
      if tar=='n': ax,alpha,color,zorder=ax12,0.8,'black', 3

      ax .fill_between(Q2,(mean-std),(mean+std),facecolor="None",edgecolor=color,hatch='|',alpha=alpha,zorder=zorder)

  hand['LT'] = ax .fill_between(Q2,-100,-100, facecolor="None",edgecolor='black',hatch='|',alpha=alpha,zorder=zorder)

  #--load JAM15 data 
  F = open('data/d2_JAM15.csv','r')
  L = F.readlines()
  F.close()

  L = [l.strip() for l in L]
  L = [[x for x in l.split()] for l in L]
  L = np.transpose(L)[0]

  Q2, pmax,pmin,pmean,nmax,nmin,nmean = [],[],[],[],[],[],[]
  for i in range(len(L)):
      if i==0: continue
      Q2.append(float(L[i].split(',')[0]))
      pmax .append(float(L[i].split(',')[1]))
      pmin .append(float(L[i].split(',')[2]))
      pmean.append(float(L[i].split(',')[3]))
      nmax .append(float(L[i].split(',')[4]))
      nmin .append(float(L[i].split(',')[5]))
      nmean.append(float(L[i].split(',')[6]))
  Q2 = np.array(Q2)
  pmax  = np.array(pmax)
  pmin  = np.array(pmin)
  pmean = np.array(pmean)
  nmax  = np.array(nmax)
  nmin  = np.array(nmin)
  nmean = np.array(nmean)
 
  alpha = 0.6
  ax11.fill_between(Q2,pmin,pmax,facecolor="None",hatch='//',edgecolor='black', alpha=alpha,zorder=3)
  ax12.fill_between(Q2,nmin,nmax,facecolor="None",hatch='//',edgecolor='black',alpha=alpha,zorder=3)
  hand['JAM15'] = ax11.fill_between(Q2,-100,-100,facecolor="None",hatch='//',edgecolor='black',alpha=alpha,zorder=3)

  for ax in [ax11,ax12]:
      ax.set_xlim(1.7,5.5)
      ax.set_xticks([2,3,4,5])
      ax.xaxis.set_minor_locator(MultipleLocator(0.5))

      ax.tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=25,length=8)
      ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=25,length=4)
      ax.set_xlabel(r'\boldmath$Q^2~({\rm GeV}^2)$' ,size=30)

      ax.axhline(0,lw = 1,color = 'k',alpha = 0.3)
      ax.set_ylim(-0.029,0.029)
      #ax.yaxis.set_minor_locator(MultipleLocator(0.001))
      #ax.yaxis.set_major_locator(MultipleLocator(0.005))


    

  ax11.text(0.03,0.83,r'\boldmath$d_2$',transform=ax11.transAxes,size=60)

  ax11.text(0.55,0.85,r'\rm \bf proton', transform=ax11.transAxes,size = 45, color = 'red')
  ax11.text(0.55,0.73,r'\rm \bf neutron',transform=ax11.transAxes,size = 45, color = 'blue')

  ax12.tick_params(labelleft=False)

  handles, labels = [],[]
  handles.append(hand['JAM15'])
  labels.append(r'\textrm{\textbf{JAM15}}')
  ax11.legend(handles,labels,fontsize = 22, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)

  handles, labels = [],[]
  handles.append(hand['LT'])
  labels.append(r'\textrm{\textbf{no HT}}')
  ax12.legend(handles,labels,fontsize = 22, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)

  py.tight_layout()
  py.subplots_adjust(hspace=0,wspace=0.01,top=0.99,right=0.99)

  filename = 'plots/fig_d2_trunc'

  filename+=ext

  ax11.set_rasterized(True)

  py.savefig(filename)
  print ('Saving figure to %s'%filename)
  py.clf()

def plot_BC_SR():

  Nick = False

  nrows,ncols=1,1
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*7,nrows*5))
  ax11 = py.subplot(nrows,ncols,1)

  load_config('%s/input.py'%msr_dir)
  istep=core.get_istep()

  TAR = ['p','n']
  hand = {}
  for tar in TAR:
      if Nick: filename = '%s/data/pstf-BC-SR-LT-TMC-HT-%s-full.dat'%(msr_dir,tar)
      else:    filename = '%s/data/pstf-BCSR-%s.dat'%(msr_dir,tar)
      #--load data if it exists
      try:
          data=load(filename)
      #--generate data and then load it if it does not exist
      except:
          if Nick==False: pstf.gen_BCSR(msr_dir,tar)
          data=load(filename)

      if Nick: Q2 = data['Q2'][:,0]
      else:    Q2 = data['Q2']
      data = data['SR'] 
      mean = np.mean(data,axis=0)
      std  = np.std (data,axis=0)

      if tar=='p': ax,color,zorder=ax11,'red' ,2
      if tar=='n': ax,color,zorder=ax11,'blue',1

      hand[tar] = ax .fill_between(Q2,(mean-std),(mean+std),color=color,alpha=0.7,zorder=zorder)


  for ax in [ax11]:
      ax.set_xlim(1.2,10)
      ax.set_xticks([2,4,6,8])
      ax.xaxis.set_minor_locator(MultipleLocator(1))

      ax.tick_params(axis='both', which='both', top=True, direction='in',labelsize=20)
      ax.tick_params(axis='both', which='major', top=False, direction='in',labelsize=25,length=10)
      ax.tick_params(axis='both', which='minor', top=False, direction='in',labelsize=0,length=5)
      ax.set_xlabel(r'\boldmath$Q^2~({\rm GeV}^2)$' ,size=30)

  ax11.axhline(0,lw = 1,color = 'k',alpha = 0.3)
    
  ax11.set_ylim(-0.49,0.39)
  ax11.yaxis.set_minor_locator(MultipleLocator(0.05))
  ax11.yaxis.set_major_locator(MultipleLocator(0.10))

  #ax11.text(0.03,0.83,r'\textrm{\textbf{BCSR}}',transform=ax11.transAxes,size=60)
  #ax11.text(0.03,0.83,r'\textrm{\textbf{Burkhardt-Cottingham sum rule}}',transform=ax11.transAxes,size=40)
  ax11.text(0.03,0.10,r'\boldmath$\int_{0}^{1} {\rm d} x g_2$',transform=ax11.transAxes,size=45)

  ax11.text(0.55,0.15,r'\rm \bf proton', transform=ax11.transAxes,size = 45, color = 'red')
  ax11.text(0.55,0.03,r'\rm \bf neutron',transform=ax11.transAxes,size = 45, color = 'blue')

  py.tight_layout()
  py.subplots_adjust(hspace=0,right=0.99,top=0.99)

  if Nick: filename = 'plots/fig_BCSR_Nick'
  else:    filename = 'plots/fig_BCSR'

  filename+=ext

  py.savefig(filename)
  print ('Saving figure to %s'%filename)
  py.clf()

def plot_BC_SR_trunc():

    xmin = 5e-3
    W2 = 4

    nrows,ncols=1,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5.5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    load_config('%s/input.py'%msr_dir)
    istep=core.get_istep()

    TAR = ['p','n']
    hand = {}
    for tar in TAR:
        filename = '%s/data/pstf-BCSR-trunc-%s-xmin=%3.5f-W2=%3.5f.dat'%(msr_dir,tar,xmin,W2)
        if regen: pstf.gen_BCSR_trunc(msr_dir,tar,xmin,W2)
        #--load data if it exists
        try:
            data=load(filename)
        #--generate data and then load it if it does not exist
        except:
            pstf.gen_BCSR_trunc(msr_dir,tar,xmin,W2)
            data=load(filename)

        Q2 = data['Q2']
        data = data['SR'] 
        mean = np.mean(data,axis=0)
        std  = np.std (data,axis=0)

        if tar=='p': ax,color,zorder=ax11,'red' ,1
        if tar=='n': ax,color,zorder=ax12,'blue',1

        label = None

        hand[tar] = ax .fill_between(Q2,(mean-std),(mean+std),color=color,alpha=0.9,zorder=zorder)

        filename = '%s/data/pstf-BCSR-trunc-%s-xmin=%3.5f-W2=%3.5f-LT.dat'%(msr_dir,tar,xmin,W2)
        if regen: pstf.gen_BCSR_trunc(msr_dir,tar,xmin,W2,LT_only=True)
        #--load data if it exists
        try:
            data=load(filename)
        #--generate data and then load it if it does not exist
        except:
            pstf.gen_BCSR_trunc(msr_dir,tar,xmin,W2,LT_only=True)
            data=load(filename)

        Q2 = data['Q2']
        data = data['SR'] 
        mean = np.mean(data,axis=0)
        std  = np.std (data,axis=0)

        if tar=='p': ax,color,zorder=ax11,'firebrick' ,3
        if tar=='n': ax,color,zorder=ax12,'darkblue',3

        label = None

        ax .fill_between(Q2,(mean-std),(mean+std),edgecolor=color,facecolor="None",alpha=1.0,zorder=zorder,hatch='//',lw=2)

    hand['LT'] = ax .fill_between(Q2,-100,-100,edgecolor='black', facecolor="None",alpha=1.0,zorder=zorder,hatch='//',lw=2)

    for ax in [ax11,ax12]:
        ax.set_xlim(1.7,10)
        ax.set_xticks([2,4,6,8])
        ax.xaxis.set_minor_locator(MultipleLocator(1))

        ax.tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=30,length=8)
        ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=30,length=4)
        ax.set_xlabel(r'\boldmath$Q^2~[{\rm GeV}^2]$' ,size=30)

        ax.axhline(0,lw = 1,color = 'k',alpha = 0.3)
      
        ax.set_ylim(-0.06,0.19)
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.set_yticks([-0.05,0,0.05,0.10,0.15])
        ax.set_yticklabels([r'$-0.05$',r'$0$',r'$0.05$',r'$0.10$',r'$0.15$'])

    ax12.tick_params(labelleft=False)

    #ax11.text(0.03,0.83,r'\textrm{\textbf{BCSR}}',transform=ax11.transAxes,size=60)
    #ax11.text(0.03,0.83,r'\textrm{\textbf{Burkhardt-Cottingham sum rule}}',transform=ax11.transAxes,size=40)
    #ax11.text(0.03,0.85,r'\boldmath$\int_{x_{\rm min}}^{x_{\rm max}} {\rm d} x g_2$',transform=ax11.transAxes,size=45)
    ax11.text(0.03,0.85,r'\boldmath$\int_{0.005}^{x_{\rm max}} {\rm d} x g_2$',transform=ax11.transAxes,size=45)

    ax11.text(0.60,0.85,r'\rm \bf proton', transform=ax11.transAxes,size = 45, color = 'red')
    ax11.text(0.60,0.73,r'\rm \bf neutron',transform=ax11.transAxes,size = 45, color = 'blue')

    handles, labels = [],[]
    handles.append(hand['LT'])
    labels.append(r'\textrm{\textbf{no HT}}')
    ax12.legend(handles,labels,fontsize = 22, loc=(0.70,0.85), frameon=0,handlelength=1.5,handletextpad=0.6)

    py.tight_layout()
    py.subplots_adjust(hspace=0,wspace=0.01,right=0.99,top=0.99)

    filename = 'plots/fig_BCSR_trunc'

    filename+=ext

    ax11.set_rasterized(True)
    ax12.set_rasterized(True)

    py.savefig(filename)
    print ('Saving figure to %s'%filename)
    py.clf()

def plot_d2_integrand(Q2=4,xmin=0.005,xmax=0.999):

    load_config('%s/input.py' % msr_dir)
    istep = core.get_istep()

    if regen:
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='n',stf='g2')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
    #--load full results
    try:
        sf_data_p_g1    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
  
    #--load LT with TMC results 
    #if regen:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1',LT_only=True)
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1',LT_only=True)
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2',LT_only=True)
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2',LT_only=True)
    #try:
    #    sf_data_p_g1_LT    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1',LT_only=True)
    #    sf_data_p_g1_LT    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f-LT.dat'%Q2)
    #try:
    #    sf_data_n_g1_LT    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1',LT_only=True)
    #    sf_data_n_g1_LT    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f-LT.dat'%Q2)
    #try:
    #    sf_data_p_g2_LT    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2',LT_only=True)
    #    sf_data_p_g2_LT    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f-LT.dat'%Q2)
    #try:
    #    sf_data_n_g2_LT    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2',LT_only=True)
    #    sf_data_n_g2_LT    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f-LT.dat'%Q2)

    #--load LT with no TMC results 
    try:
        sf_data_p_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    
    
    data = {}
    data['g1p'] = sf_data_p_g1
    data['g1n'] = sf_data_n_g1
    data['g2p'] = sf_data_p_g2
    data['g2n'] = sf_data_n_g2

    data_LT = {}
    data_LT['g1p'] = sf_data_p_g1_LT
    data_LT['g1n'] = sf_data_n_g1_LT
    data_LT['g2p'] = sf_data_p_g2_LT
    data_LT['g2n'] = sf_data_n_g2_LT

    #data_LT_noTMC = {}
    #data_LT_noTMC['g1p'] = sf_data_p_g1_LT_noTMC
    #data_LT_noTMC['g1n'] = sf_data_n_g1_LT_noTMC
    #data_LT_noTMC['g2p'] = sf_data_p_g2_LT_noTMC
    #data_LT_noTMC['g2n'] = sf_data_n_g2_LT_noTMC
    
    ncols = 2
    nrows = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*8,nrows*5))
    
    
    i = 1
    hand = {} 
    for tar in ['p', 'n']:
        j = 0

        if tar=='p': ax, zorder, color = axs[0][0],1,'red'
        if tar=='n': ax, zorder, color = axs[0][1],1,'dodgerblue'
        X = data['g1p']['X']
        g1 = data['g1'+tar]['STF']
        g2 = data['g2'+tar]['STF']
        func = X**2 * (2*g1 + 3*g2) 
        mean = np.mean(func,axis=0)
        std  = np.std (func,axis=0)
        hand['JAM %s'%tar] = ax.fill_between(X,mean - std,mean + std, color = color,alpha = 0.9,zorder=zorder)

        
        #--plot leading twist with or without AOT
        g1_LT = data_LT['g1'+tar]['STF']
        g2_LT = data_LT['g2'+tar]['STF']
        func_LT = X**2 * (2*g1_LT + 3*g2_LT) 
        mean_LT = np.mean(func_LT,axis=0)
        std_LT  = np.std (func_LT,axis=0)

        if tar=='p': ax, zorder,color = axs[0][0], 1.1,'firebrick'
        if tar=='n': ax, zorder,color = axs[0][1], 1.1,'blue'

        hand['LT %s'%tar] = ax .fill_between(X,mean_LT - std_LT,mean_LT + std_LT, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//',lw=2)

        #--plot truncated moment as function of xmax
        j = 1
        if tar=='p': ax, zorder, color = axs[1][0], 1,'red'
        if tar=='n': ax, zorder, color = axs[1][1], 1,'dodgerblue'
        d2 = cumulative_trapezoid(func,X,initial=0.0)
        d2_mean = np.mean(d2,axis=0)
        d2_std  = np.std (d2,axis=0) 
        hand['JAM %s'%tar] = ax.fill_between(X,d2_mean - d2_std,d2_mean + d2_std, color = color,alpha = 0.9,zorder=zorder)

        if tar=='p': ax, zorder,color = axs[1][0], 1.1,'firebrick'
        if tar=='n': ax, zorder,color = axs[1][1], 1.1,'blue'
        d2_LT = cumulative_trapezoid(func_LT,X,initial=0.0)
        d2_mean_LT = np.mean(d2_LT,axis=0)
        d2_std_LT  = np.std (d2_LT,axis=0) 

        hand['LT %s'%tar] = ax.fill_between(X,d2_mean_LT - d2_std_LT,d2_mean_LT + d2_std_LT, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//',lw=2)
    
    for i in range(2):
        for j in range(2):
            axs[i][j].tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=30,length=8)
            axs[i][j].tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=0,length=4)
            axs[i][j].axhline(0,lw = 1,color = 'k',alpha = 0.3)
            
            
            axs[i][j].set_xlim(0.01,0.99)
            axs[i][j].set_xticks([0.1,0.3,0.5,0.7,0.9])
            axs[i][j].set_xticklabels([r'$0.1$',r'$0.3$',r'$0.5$',r'$0.7$',r'$0.9$'])
            axs[i][j].axhline(0,lw = 1,color = 'k',alpha = 0.3)

            axs[i][j].xaxis.set_minor_locator(MultipleLocator(0.10))

    #--plot extrapolation region
    M2 = 0.938**3
    W2min = 4
    xlim = Q2/(Q2 + W2min - M2)
    x = np.linspace(xlim,1,10)
    do = -100*np.ones(len(x))
    up =  100*np.ones(len(x))
    axs[0][0].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)
    axs[0][1].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)
    axs[1][0].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)
    axs[1][1].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)

    #axs[1][0].text(0.60, 0.90, r'\textrm{extrapolation}', transform=axs[1][0].transAxes,fontsize=30,alpha=0.5,rotation=0)


    axs[0][0].set_xlabel(r'\boldmath{$x$}',fontsize = 45)
    axs[0][1].set_xlabel(r'\boldmath{$x$}',fontsize = 45)
    axs[1][0].set_xlabel(r'\boldmath{$x_{\rm max}$}',fontsize = 45)
    axs[1][1].set_xlabel(r'\boldmath{$x_{\rm max}$}',fontsize = 45)
 
    axs[0][0].set_ylim(-0.025,0.019)
    axs[0][0].yaxis.set_minor_locator(MultipleLocator(0.005))
    axs[0][0].set_yticks([-0.02,-0.01,0,0.01])
    axs[0][0].set_yticklabels([r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$'])

    axs[0][1].set_ylim(-0.04,0.04)
    axs[0][1].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[0][1].set_yticks([-0.02,0,0.02])
    axs[0][1].set_yticklabels([r'$-0.02$',r'$0$',r'$0.02$'])

    axs[1][0].set_ylim(-0.001,0.005)
    axs[1][0].yaxis.set_minor_locator(MultipleLocator(0.001))
    axs[1][0].set_yticks([0,0.002,0.004])
    axs[1][0].set_yticklabels([r'$0$',r'$0.002$',r'$0.004$'])
    
    axs[1][1].set_ylim(-0.005,0.009)
    axs[1][1].yaxis.set_minor_locator(MultipleLocator(0.001))
    axs[1][1].set_yticks([-0.004,0,0.004,0.008])
    axs[1][1].set_yticklabels([r'$-0.004$',r'$0$',r'$0.004$',r'$0.008$'])
    
    # ##########################################################

    #axs[0][1].text(0.03,0.85,r'\rm \bf JAMpol25', transform=axs[0][1].transAxes,size = 40)
       
    axs[0][0].text(0.60,0.85,r'\rm \bf proton', transform=axs[0][0].transAxes,size = 40, color = 'red')
    axs[0][1].text(0.05,0.85,r'\rm \bf neutron',transform=axs[0][1].transAxes,size = 40, color = 'blue')
 
    #axs[0][0].text(0.03,0.05,r'\rm \bf proton', transform=axs[0][0].transAxes,size = 40)
    #axs[0][1].text(0.03,0.05,r'\rm \bf neutron',transform=axs[0][1].transAxes,size = 40)
    #axs[1][0].text(0.03,0.85,r'\rm \bf proton', transform=axs[1][0].transAxes,size = 40)
    #axs[1][1].text(0.03,0.85,r'\rm \bf neutron',transform=axs[1][1].transAxes,size = 40)
    
    axs[0][0].set_ylabel(r'\boldmath$x^2 (2 g_1 + 3 g_2)$',size=40)   
    axs[1][0].set_ylabel(r'\boldmath$\int_{%s}^{x_{\rm max}} {\rm d} x x^2 (2 g_1 + 3 g_2)$'%xmin,size=30)   
 
    #axs[0][0].text(0.03, 0.80, r'\boldmath$x^2 (2 g_1 + 3 g_2)$', transform=axs [0][0].transAxes,fontsize=50,color = 'k')
    #axs[1][1].text(0.03, 0.15, r'\boldmath$\int_{0.005}^{x_{\rm max}} {\rm d} x x^2 (2 g_1 + 3 g_2)$', transform=axs [1][1].transAxes,fontsize=40,color = 'k')
   
    axs[1][1].text(0.03,0.05,r'$Q^2=%d~{\rm GeV}^2$'%Q2,transform=axs[1][1].transAxes,size=25)


    hand['LT'] = axs[1][1] .fill_between(X,-100,-99, color = 'None', edgecolor='black',alpha = 1.0,zorder=3,hatch='//',lw=2)
 
    handles, labels = [],[]
    handles.append(hand['JAM p'])
    handles.append(hand['LT p'])
    labels.append(r'\textrm{\textbf{JAMpol25}}')
    labels.append(r'\textrm{\textbf{LT only}}')
    axs[0][0].legend(handles,labels,fontsize = 25, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)

    handles, labels = [],[]
    handles.append(hand['JAM n'])
    handles.append(hand['LT n'])
    labels.append(r'\textrm{\textbf{JAMpol25}}')
    labels.append(r'\textrm{\textbf{LT only}}')
    axs[0][1].legend(handles,labels,fontsize = 25, loc=(0.00,0.00), frameon=0,handlelength=1,handletextpad=0.6)
 
    py.tight_layout()
    py.subplots_adjust(top=0.99,right=0.99)#,wspace=0.01)
    
    axs[0][0].set_rasterized(True)
    axs[0][1].set_rasterized(True)
    axs[1][0].set_rasterized(True)
    axs[1][1].set_rasterized(True)
    
    #filename = 'plots/fig_d2_integrand_Q2=%3.5f'%Q2 + ext
    filename = 'plots/fig_d2_integrand' + ext
    py.savefig(filename)
    py.close()
    print('Saving figure to %s'%filename)

def plot_Bjorken_sumrule(Q2=4,xmin=0.005,xmax=0.999):

    load_config('%s/input.py' % msr_dir)
    istep = core.get_istep()

    if regen:
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir   ,Q2=Q2,tar='n',stf='g2')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
    #--load full results
    try:
        sf_data_p_g1    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
  
    #--load LT with TMC results 
    #if regen:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1',LT_only=True)
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1',LT_only=True)
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2',LT_only=True)
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2',LT_only=True)
    #try:
    #    sf_data_p_g1_LT    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g1',LT_only=True)
    #    sf_data_p_g1_LT    = load(msr_dir    + 'data/pstf-g1-p-Q2=%3.5f-LT.dat'%Q2)
    #try:
    #    sf_data_n_g1_LT    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g1',LT_only=True)
    #    sf_data_n_g1_LT    = load(msr_dir    + 'data/pstf-g1-n-Q2=%3.5f-LT.dat'%Q2)
    #try:
    #    sf_data_p_g2_LT    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='p',stf='g2',LT_only=True)
    #    sf_data_p_g2_LT    = load(msr_dir    + 'data/pstf-g2-p-Q2=%3.5f-LT.dat'%Q2)
    #try:
    #    sf_data_n_g2_LT    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f-LT.dat'%Q2)
    #except:
    #    pstf.gen_pstf(msr_dir,Q2=Q2,tar='n',stf='g2',LT_only=True)
    #    sf_data_n_g2_LT    = load(msr_dir    + 'data/pstf-g2-n-Q2=%3.5f-LT.dat'%Q2)

    #--load LT with no TMC results 
    try:
        sf_data_p_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    
    
    data = {}
    data['g1p'] = sf_data_p_g1
    data['g1n'] = sf_data_n_g1

    data_LT = {}
    data_LT['g1p'] = sf_data_p_g1_LT
    data_LT['g1n'] = sf_data_n_g1_LT

    #data_LT_noTMC = {}
    #data_LT_noTMC['g1p'] = sf_data_p_g1_LT_noTMC
    #data_LT_noTMC['g1n'] = sf_data_n_g1_LT_noTMC
    
    ncols = 1
    nrows = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*8,nrows*5))
    
    i = 1
    hand = {} 
    j = 0

    zorder, color = 1.0,'red'
    X = data['g1p']['X']
    if xmin < X[0]: xmin = X[0]
    g1p = data['g1p']['STF']
    g1n = data['g1n']['STF']
    func = g1p - g1n 
    mean = np.mean(func,axis=0)
    std  = np.std (func,axis=0)
    hand['JAM'] = axs[j] .fill_between(X,mean - std,mean + std, color = color,alpha = 0.9,zorder=zorder)

    #--plot leading twist
    g1p_LT = data_LT['g1p']['STF']
    g1n_LT = data_LT['g1n']['STF']
    func_LT = g1p_LT - g1n_LT 
    mean_LT = np.mean(func_LT,axis=0)
    std_LT  = np.std (func_LT,axis=0)

    zorder,color = 1.1,'black'

    hand['LT'] = axs[j] .fill_between(X,mean_LT - std_LT,mean_LT + std_LT, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//',lw=2)

    #--plot truncated moment as function of xmax
    #--set results below xmin = 0.005 to zero
    func   [:,X < xmin] = 0
    func_LT[:,X < xmin] = 0
    j = 1
    zorder, color = 1.0,'red'
    cum = cumulative_trapezoid(func,X,initial=0.0)
    cum_mean = np.mean(cum,axis=0)
    cum_std  = np.std (cum,axis=0) 
    axs[j] .fill_between(X,cum_mean - cum_std,cum_mean + cum_std, color = color,alpha = 0.9,zorder=zorder)

    zorder,color = 1.1,'black'
    cum_LT = cumulative_trapezoid(func_LT,X,initial=0.0)
    cum_mean_LT = np.mean(cum_LT,axis=0)
    cum_std_LT  = np.std (cum_LT,axis=0) 

    hand['LT'] = axs[j] .fill_between(X,cum_mean_LT - cum_std_LT,cum_mean_LT + cum_std_LT, color = 'None', edgecolor=color,alpha = 1.0,zorder=zorder,hatch='//',lw=2)
    
    for i in range(2):
        axs[i].tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=30,length=8)
        axs[i].tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=0,length=4)
        axs[i].axhline(0,lw = 1,color = 'k',alpha = 0.3)
        
        axs[i].semilogx() 
        axs[i].set_xlim(0.005,1)
        #axs[i].set_xlim(0.01,1.00)
        axs[i].set_xticks([1e-2,1e-1,1])
        axs[i].set_xticklabels([r'$0.01$',r'$0.1$',r'$1$'])
        #axs[i].axhline(0,lw = 1,color = 'k',alpha = 0.3)

        #axs[i].xaxis.set_minor_locator(MultipleLocator(0.10))

    axs[0].set_xlabel(r'\boldmath{$x$}',fontsize = 45)
    axs[1].set_xlabel(r'\boldmath{$x_{\rm max}$}',fontsize = 45)
 
    axs[0].set_ylim(-0.1,1.9)
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.1))
    axs[0].set_yticks([0,0.5,1.0,1.5])
    axs[0].set_yticklabels([r'$0$',r'$0.5$',r'$1.0$',r'$1.5$'])

    axs[1].set_ylim(0,0.20)
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[1].set_yticks([0.05,0.10,0.15,0.20])
    axs[1].set_yticklabels([r'$0.05$',r'$0.10$',r'$0.15$',r'$0.20$'])
  
    #--get and plot result from theory
    gA     = 1.269
    gA_err = 0.003
    RESMAN(nworkers=1,parallel=False,datasets=False) 
    alphaS = conf['alphaS'].get_alphaS(Q2)
    theory_LO    = gA/6
    theory_NLO   = gA/6 * (1 - alphaS/np.pi)
    theory_NNLO  = gA/6 * (1 - alphaS/np.pi - 3.5833 * (alphaS/np.pi)**2)
    theory_NNNLO = gA/6 * (1 - alphaS/np.pi - 3.5833 * (alphaS/np.pi)**2 - 20.2153 * (alphaS/np.pi)**3)
    print('LO: %4.3f'%theory_LO)
    print('NLO: %4.3f'%theory_NLO)
    print('NNLO: %4.3f'%theory_NNLO)
    print('NNNLO: %4.3f'%theory_NNNLO)
    hand['theory'] = axs[1].axhline(theory_NLO,0.85,1.0,color='blue',ls='--',alpha=1.0,zorder=1.1,lw=3)
   
    #--plot extrapolation region
    M2 = 0.938**3
    W2min = 4
    xlim = Q2/(Q2 + W2min - M2)
    x = np.linspace(xlim,1,10)
    do = -100*np.ones(len(x))
    up =  100*np.ones(len(x))
    axs[0].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)
    axs[1].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)

    #axs[1].text(0.95, 0.05, r'\textrm{extrapolation}', transform=axs [1].transAxes,fontsize=30,alpha=0.5,rotation=90)
 
    # ##########################################################
        
    #axs[1].text(0.03,0.85,r'\rm \bf proton', transform=axs[1].transAxes,size = 40, color = 'red')
    #axs[1].text(0.03,0.75,r'\rm \bf neutron',transform=axs[1].transAxes,size = 40, color = 'blue')
    
    
    axs[0].text(0.50, 0.80, r'\boldmath$g_1^p - g_1^n$', transform=axs [0].transAxes,fontsize=50,color = 'k')
    axs[1].text(0.03, 0.80, r'\boldmath$\int_{%s}^{x_{\rm max}} {\rm d}x (g_1^p - g_1^n)$'%xmin, transform=axs [1].transAxes,fontsize=40,color = 'k')
   
    axs[0].text(0.03,0.10,r'$Q^2=%d~{\rm GeV}^2$'%Q2,transform=axs[0].transAxes,size=25)


 
    handles, labels = [],[]
    handles.append(hand['JAM'])
    handles.append(hand['LT'])
    handles.append(hand['theory'])
    labels.append(r'\textrm{\textbf{JAMpol25}}')
    labels.append(r'\textrm{\textbf{JAMpol25 LT}}')
    labels.append(r'\textrm{\textbf{LT theory}}')
    axs[1].legend(handles,labels,fontsize = 25, loc=(0.00,0.25), frameon=0,handlelength=1,handletextpad=0.6)
 
    py.tight_layout()
    py.subplots_adjust(top=0.99,right=0.99,wspace=0.01)
    
    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
    
    filename = 'plots/fig_bjorken_sumrule' + ext
    py.savefig(filename)
    py.close()
    print('Saving figure to %s'%filename)

def plot_Bjorken_sumrule_nobuo(Q2=4):

    xmax = 0.999
    XMIN = [0.005, 1e-5]
    load_config('%s/input.py' % msr_dir)
    istep = core.get_istep()

    #--load LT only results? 
    if regen:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
    try:
        sf_data_p_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g1')
        sf_data_p_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g1')
        sf_data_n_g1_LT    = load(msr_dir_LT    + 'data/pstf-g1-n-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_p_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='p',stf='g2')
        sf_data_p_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-p-Q2=%3.5f.dat'%Q2)
    try:
        sf_data_n_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    except:
        pstf.gen_pstf(msr_dir_LT,Q2=Q2,tar='n',stf='g2')
        sf_data_n_g2_LT    = load(msr_dir_LT    + 'data/pstf-g2-n-Q2=%3.5f.dat'%Q2)
    
    

    data_LT = {}
    data_LT['g1p'] = sf_data_p_g1_LT
    data_LT['g1n'] = sf_data_n_g1_LT
    
    ncols = 1
    nrows = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*10,nrows*5))
    
    i = 1
    hand = {} 
    j = 0

    zorder, color = 1.0,'red'
    X = data_LT['g1p']['X']
    g1p_LT = data_LT['g1p']['STF']
    g1n_LT = data_LT['g1n']['STF']
    func_LT = g1p_LT - g1n_LT 
    mean_LT = np.mean(func_LT,axis=0)
    std_LT  = np.std (func_LT,axis=0)

    zorder,color = 1.1,'firebrick'

    hand['LT'] = axs[j] .fill_between(X,mean_LT - std_LT,mean_LT + std_LT, color = 'firebrick',alpha = 1.0,zorder=zorder)

    #--plot truncated moment as function of xmax
    j = 1
    k = 0
    #--set results below xmin = 0.005 to zero
    for xmin in XMIN:
        func_LT_trunc = copy.copy(func_LT)
        func_LT_trunc[:,X < xmin] = 0

        if k==0: zorder,color = 1.2,'firebrick'
        if k==1: zorder,color = 1.1,'darkgreen'
        if k==2: zorder,color = 1.0,'blue'
        cum_LT = cumulative_trapezoid(func_LT_trunc,X,initial=0.0)
        cum_mean_LT = np.mean(cum_LT,axis=0)
        cum_std_LT  = np.std (cum_LT,axis=0) 

        hand[k] = axs[j] .fill_between(X,cum_mean_LT - cum_std_LT,cum_mean_LT + cum_std_LT, color =color,alpha = 0.7,zorder=zorder)
    
        k+=1

    for i in range(2):
        axs[i].tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=30,length=8)
        axs[i].tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=0,length=4)
        #axs[i].axhline(0,lw = 1,color = 'k',alpha = 0.3)
        
        axs[i].semilogx() 
        axs[i].set_xlim(1e-5,1)
        #axs[i].set_xlim(0.01,1.00)
        axs[i].set_xticks([1e-4,1e-3,1e-2,1e-1,1])
        axs[i].set_xticklabels([r'$0.0001$',r'$0.001$',r'$0.01$',r'$0.1$',r'$1$'])
        #axs[i].axhline(0,lw = 1,color = 'k',alpha = 0.3)

        #axs[i].xaxis.set_minor_locator(MultipleLocator(0.10))

    axs[0].set_xlabel(r'\boldmath{$x$}',fontsize = 45)
    axs[1].set_xlabel(r'\boldmath{$x_{\rm max}$}',fontsize = 45)
 
    axs[0].set_ylim(0,10)
    axs[0].yaxis.set_minor_locator(MultipleLocator(1))
    axs[0].set_yticks([2,4,6,8])
    axs[0].set_yticklabels([r'$2$',r'$4$',r'$6$',r'$8$'])

    axs[1].set_ylim(0,0.22)
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.01))
    axs[1].set_yticks([0.05,0.10,0.15,0.20])
    axs[1].set_yticklabels([r'$0.05$',r'$0.10$',r'$0.15$',r'$0.20$'])
  
    #--get and plot result from theory
    gA     = 1.269
    gA_err = 0.003
    RESMAN(nworkers=1,parallel=False,datasets=False) 
    alphaS = conf['alphaS'].get_alphaS(Q2)
    theory_LO    = gA/6
    theory_NLO   = gA/6 * (1 - alphaS/np.pi)
    theory_NNLO  = gA/6 * (1 - alphaS/np.pi - 3.5833 * (alphaS/np.pi)**2)
    theory_NNNLO = gA/6 * (1 - alphaS/np.pi - 3.5833 * (alphaS/np.pi)**2 - 20.2153 * (alphaS/np.pi)**3)
    #print('LO: %4.3f'%theory_LO)
    #print('NLO: %4.3f'%theory_NLO)
    #print('NNLO: %4.3f'%theory_NNLO)
    #print('NNNLO: %4.3f'%theory_NNNLO)
    hand['theory'] = axs[1].axhline(theory_NLO,0.85,1.0,color='blue',ls='--',alpha=1.0,zorder=1.1,lw=3)
   
    #--plot extrapolation region
    M2 = 0.938**3
    W2min = 4
    xlim = Q2/(Q2 + W2min - M2)
    x = np.linspace(xlim,1,10)
    do = -100*np.ones(len(x))
    up =  100*np.ones(len(x))
    axs[0].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)
    axs[1].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)

    xlim = 0.005
    x = np.linspace(0,xlim,10)
    do = -100*np.ones(len(x))
    up =  100*np.ones(len(x))
    axs[0].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)
    axs[1].fill_between(x,do,up,color='gray',alpha=0.1,zorder=0.9)

    #axs[1].text(0.95, 0.05, r'\textrm{extrapolation}', transform=axs [1].transAxes,fontsize=30,alpha=0.5,rotation=90)
 
    # ##########################################################
        
    #axs[1].text(0.03,0.85,r'\rm \bf proton', transform=axs[1].transAxes,size = 40, color = 'red')
    #axs[1].text(0.03,0.75,r'\rm \bf neutron',transform=axs[1].transAxes,size = 40, color = 'blue')
    
    
    #axs[0].text(0.65, 0.80, r'\boldmath$g_1^p - g_1^n$', transform=axs [0].transAxes,fontsize=50,color = 'k')
    #axs[0].text(0.65, 0.55, r'\textrm{\textbf{LT only}}', transform=axs [0].transAxes,fontsize=50,color = 'k')
    
    #axs[1].text(0.03, 0.80, r'\boldmath$\int_{x_{\rm min}}^{x_{\rm max}} {\rm d}x (g_1^p - g_1^n)$', transform=axs [1].transAxes,fontsize=40,color = 'k')
   
    axs[0].text(0.75,0.90,r'$Q^2=%d~{\rm GeV}^2$'%Q2,transform=axs[0].transAxes,size=20)

    axs[0].set_ylabel(r'\boldmath$g_1^p - g_1^n ~ ({\rm LT ~ only})$',size=30,labelpad=40)   
    axs[1].set_ylabel(r'\boldmath$\int_{x_{\rm min}}^{x_{\rm max}} {\rm d}x (g_1^p - g_1^n)$',size=30,labelpad=20)   

 
    handles, labels = [],[]
    handles.append(hand[0])
    handles.append(hand[1])
    #handles.append(hand[2])
    handles.append(hand['theory'])
    labels.append(r'\boldmath$x_{\rm min} = 0.005$')
    labels.append(r'\boldmath$x_{\rm min} = 0.00001$')
    #labels.append(r'\boldmath$x_{\rm min} = 0.00001$')
    labels.append(r'\textrm{\textbf{LT theory}}')
    axs[1].legend(handles,labels,fontsize = 30, loc=(0.00,0.48), frameon=0,handlelength=1,handletextpad=0.6)
 
    py.tight_layout()
    py.subplots_adjust(top=0.99,right=0.99,wspace=0.01)
    
    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
    
    filename = 'plots/fig_bjorken_sumrule' + ext

    py.savefig(filename)
    py.close()
    print('Saving figure to %s'%filename)

if __name__=="__main__":

    plot_pstfs       (Q2=4)
    plot_d2_integrand(Q2=4)

    plot_Bjorken_sumrule_nobuo(Q2=4)


    plot_gXres_W2()
    plot_d2()
    plot_BC_SR_trunc()
    





