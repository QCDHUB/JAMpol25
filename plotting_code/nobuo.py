#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd

# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

#--matplotlib
import matplotlib
import matplotlib as plt
matplotlib.use('Agg')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
import pylab as py
from matplotlib.ticker import MultipleLocator
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

#--from scipy stack 
from scipy.integrate import fixed_quad
from scipy import interpolate

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#-- from qcdlib
from qcdlib import aux

#--from local
from analysis.corelib import core
# from analysis.corelib import classifier
from analysis.qpdlib import pdf,ppdf


#--from obslib
from obslib.psidis.reader import READER

import kmeanconf as kc

ext = '.png'
#ext = '.pdf'

cwd = os.getcwd()

def plot_solid(wdir,kc,had='pi+'):

    print('\ngenerating %s PSIDIS plots from %s'%(had,wdir))
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))

    nrows,ncols=2,4
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax13 = py.subplot(nrows,ncols,3)
    ax14 = py.subplot(nrows,ncols,4)
    ax21 = py.subplot(nrows,ncols,5)
    ax22 = py.subplot(nrows,ncols,6)
    ax23 = py.subplot(nrows,ncols,7)
    ax24 = py.subplot(nrows,ncols,8)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][30000]='psidis/expdata/30000.xlsx'
    conf['datasets']['psidis']['xlsx'][30001]='psidis/expdata/30001.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    data = predictions['reactions']['psidis']

    #--get theory by seperating solutions and taking mean
    for idx in tabs:
        predictions = copy.copy(data[idx]['prediction-rep'])
        del data[idx]['prediction-rep']
        del data[idx]['residuals-rep']
        del data[idx]['shift-rep']
        del data[idx]['rres-rep']
        del data[idx]['r-residuals']
        del data[idx]['n-residuals']
        predictions_ic = [predictions[i] for i in range(len(predictions))]
        data[idx]['thy']  = np.mean(predictions_ic,axis=0)
        data[idx]['dthy'] = np.std(predictions_ic,axis=0)


    if had=='pi+': idx = 30000
    if had=='pi-': idx = 30001

    Q2bins = [1.5,2.5,3.5,4.5,5.5,7]
    zbins  = [0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575]

    Q2bins = np.unique(data[idx]['Q2'])
    zbins  = np.unique(data[idx]['z'])

    #--set up bins
    DATA = {}
    for i in range(len(zbins)):
        DATA[i] = {}
        for j in range(len(Q2bins)):
            query = "z==%s and Q2==%s"%(zbins[i],Q2bins[j])
            DATA[i][j] = pd.DataFrame(data[idx]).query(query)
    
    #######################
    #--plot absolute values
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in tabs:
        for i in range(len(zbins)):
            for j in range(len(Q2bins)):
                tab = DATA[i][j]
                X = tab['X']
                values = tab['value']
                alpha  = tab['alpha']
                if i==0: ax = ax11
                if i==1: ax = ax12
                if i==2: ax = ax13
                if i==3: ax = ax14
                if i==4: ax = ax21
                if i==5: ax = ax22
                if i==6: ax = ax23
                if i==7: ax = ax24
                if j==0: color,fmt,ms = 'firebrick','o',4
                if j==1: color,fmt,ms = 'darkgreen','^',4
                if j==2: color,fmt,ms = 'black'    ,'*',4
                if j==3: color,fmt,ms = 'cyan'     ,'s',4
                if j==4: color,fmt,ms = 'magenta'  ,'v',4
                if j==5: color,fmt,ms = 'orange'   ,'D',4

                hand[j] = ax.errorbar(X,values,yerr=alpha,color=color,fmt=fmt,ms=ms,capsize=0.0)

                thy = tab['thy']
                std = tab['dthy']
                down = thy - std
                up   = thy + std
                thy_plot[j] ,= ax.plot(X,thy,color=color)
                thy_band[j]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

                if j==0 and had=='pi+': ax.text(0.02, 0.88, r'\boldmath$z = %3.3f$'%(zbins[i]),transform=ax.transAxes,size=30)
                if j==0 and had=='pi-': ax.text(0.02, 0.04, r'\boldmath$z = %3.3f$'%(zbins[i]),transform=ax.transAxes,size=30)
  
    for ax in [ax11,ax12,ax13,ax14,ax21,ax22,ax23,ax24]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_xlim(0.03,0.57)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,color='black',alpha=0.2)

    for ax in [ax12,ax13,ax14,ax22,ax23,ax24]:
        ax.tick_params(labelleft=False)

    for ax in [ax11,ax12,ax13,ax14]:
        ax.tick_params(labelbottom=False)
        if had=='pi+':
            ax.set_ylim(-0.06,0.06)
            ax.set_yticks([-0.05,0,0.05])
            ax.set_yticklabels([r'$-0.05$',r'$0$',r'$0.05$'])
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        if had=='pi-':
            ax.set_ylim(-0.09,-0.01)
            ax.set_yticks([-0.08,-0.06,-0.04,-0.02])
            #ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))

    for ax in [ax21,ax22,ax23,ax24]:
        # ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=40)
        ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=40)
        # ax.xaxis.set_label_coords(0.94,0.00)
        if had=='pi+': 
            ax.set_ylim(-0.04,0.10)
            ax.set_yticks([0,0.04,0.08])
            ax.set_yticklabels([r'$0$',r'$0.04$',r'$0.08$'])
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        if had=='pi-': 
            ax.set_ylim(-0.09,-0.01)
            ax.set_yticks([-0.08,-0.06,-0.04,-0.02])
            #ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
            ax.yaxis.set_minor_locator(MultipleLocator(0.01))

    #for ax in [ax12,ax22]:
    #    ax.tick_params(axis='both',which='both',labelleft=False)



    if had=='pi+': ax11.text(0.60, 0.13, r'\boldmath$A^{\pi^+,^3{\rm He}}_{\parallel}$',transform=ax11.transAxes,size=50)
    if had=='pi-': ax11.text(0.02, 0.75, r'\boldmath$A^{\pi^-,^3{\rm He}}_{\parallel}$',transform=ax11.transAxes,size=50)

    handles, labels = [],[]
    for j in [0,1]:
        handles.append((thy_plot[j],thy_band[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    if had=='pi+': ax12.legend(handles,labels,frameon=False,fontsize=20,loc=(0.72,0.00),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)
    if had=='pi-': ax12.legend(handles,labels,frameon=False,fontsize=20,loc=(0.00,0.70),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)


    handles, labels = [],[]
    for j in [2,3]:
        handles.append((thy_plot[j],thy_band[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    if had=='pi+': ax13.legend(handles,labels,frameon=False,fontsize=20,loc=(0.72,0.00),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)
    if had=='pi-': ax13.legend(handles,labels,frameon=False,fontsize=20,loc=(0.00,0.70),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)

    handles, labels = [],[]
    for j in [4,5]:
        handles.append((thy_plot[j],thy_band[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])

    if had=='pi+': ax14.legend(handles,labels,frameon=False,fontsize=20,loc=(0.70,0.00),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)
    if had=='pi-': ax14.legend(handles,labels,frameon=False,fontsize=20,loc=(0.00,0.70),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)
    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/nobuo_solid_%s'%(cwd,had)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving SoLID PSIDIS %s plot to %s'%(had,filename))

def asym_impact(WDIR,Q2=10):

    nrows,ncols=2,1
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    thy  = {}
    hand = {}

    j = 0
    for wdir in WDIR:
        load_config('%s/input.py'%wdir)

        filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
        data=load(filename)
        print('Loading %s'%filename)

        X=data['X']

        asym = np.array(data['XF']['ub']) - np.array(data['XF']['db'])
        add  = np.array(data['XF']['ub']) + np.array(data['XF']['db'])

        if j==0: color,alpha,zorder='red'   ,0.9,1.1
        if j==1: color,alpha,zorder='cyan'  ,0.9,1.0
        #--plot average and standard deviation
        mean = np.mean(asym,axis=0)
        std  = np.std (asym,axis=0)
        thy[j] = ax11.fill_between(X,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)
        mean = np.mean(add,axis=0)
        std  = np.std (add,axis=0)
        thy[j] = ax12.fill_between(X,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)
        j+=1


    for ax in [ax11,ax12]:
        ax.set_xlim(0.00,0.60)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))

        ax.tick_params(axis='both', which='major', top=True, direction='in',labelsize=25,length=8)
        ax.tick_params(axis='both', which='minor', top=True, direction='in',labelsize=25,length=4)
        ax.axhline(0  ,color='k',alpha=0.5,zorder=6)


    ax11.set_ylim(-0.01,0.055)
    ax11.set_yticks([0,0.02,0.04])
    ax11.set_yticklabels([r'$0$',r'$0.02$',r'$0.04$'])
    ax11.yaxis.set_minor_locator(MultipleLocator(0.01)) 

    ax12.set_ylim(-0.05,0.05)
    ax12.set_yticks([-0.04,-0.02,0,0.02,0.04])
    ax12.set_yticklabels([r'$-0.04$',r'$-0.02$',r'$0$',r'$0.02$',r'$0.04$'])
    ax12.yaxis.set_minor_locator(MultipleLocator(0.01)) 

    ax11.tick_params(labelbottom=False)

    ax12.set_xlabel(r'\boldmath$x$',size=45,loc = 'left')
    ax12.xaxis.set_label_coords(0.92,0.00)

    ax11.text(0.35,0.85,r'\boldmath{$x (\Delta \bar{u} - \Delta \bar{d})$}',                    transform=ax11.transAxes,size=45)
    ax12.text(0.35,0.85,r'\boldmath{$x (\Delta \bar{u} + \Delta \bar{d})$}',                    transform=ax12.transAxes,size=45)

    ax11.text(0.02,0.05,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=ax11.transAxes,size=25)

    handles,labels = [],[]
    handles.append(thy[1])
    handles.append(thy[0])
    labels.append(r'\textrm{\textbf{JAM}}')
    labels.append(r'\textrm{\textbf{+SoLID}}')
    ax12.legend(handles,labels,loc=(0.50,0.00),fontsize=35,frameon=0,handletextpad=0.3,handlelength=1.0)
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.2,hspace=0.02,top=0.99,right=0.99)

    filename = 'plots/nobuo_asym_impact'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def sea_impact(WDIR,Q2=10):

    nrows,ncols=1,1
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax11 = py.subplot(nrows,ncols,1)

    thy  = {}
    hand = {}

    j = 0
    for wdir in WDIR:
        load_config('%s/input.py'%wdir)

        filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
        data=load(filename)
        print('Loading %s'%filename)

        X=data['X']

        sea  = np.array(data['XF']['s']) #--note that s = sea in this parameterization
        mean = np.mean(sea,axis=0)
        std  = np.std (sea,axis=0)

        if j==0: color,alpha,zorder='red'   ,0.9,1.1
        if j==1: color,alpha,zorder='cyan'  ,0.9,1.0
        #--plot average and standard deviation
        thy[j] = ax11.fill_between(X,(mean-std),(mean+std),color=color,alpha=alpha,zorder=zorder)
        j+=1


    ax11.set_xlim(0.15,0.60)
    ax11.set_xticks([0.2,0.3,0.4,0.5])
    ax11.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax11.tick_params(axis='both', which='major', top=True, direction='in',labelsize=25,length=8)
    ax11.tick_params(axis='both', which='minor', top=True, direction='in',labelsize=25,length=4)


    #for i in [1,2,3,4]:
    #      axs[i] .tick_params(labelbottom=False)
    #      axLs[i].tick_params(labelbottom=False)

    #

    ax11.set_ylim(-0.09,0.09)
    ax11.set_yticks([-0.08,-0.04,0,0.04,0.08])
    ax11.set_yticklabels([r'$-0.08$',r'$-0.04$',r'$0$',r'$0.04$',r'$0.08$'])
    ax11.yaxis.set_minor_locator(MultipleLocator(0.01)) 

    ax11.axhline(0  ,color='k',alpha=0.5,zorder=6)

    ax11.set_xlabel(r'\boldmath$x$',size=45,loc = 'left')
    ax11.xaxis.set_label_coords(0.92,0.00)

    ax11.text(0.20,0.85,r'\textrm{\textbf{Symmetric Sea}}',                    transform=ax11.transAxes,size=45)

    ax11.text(0.02,0.05,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=ax11.transAxes,size=25)

    handles,labels = [],[]
    handles.append(thy[1])
    handles.append(thy[0])
    labels.append(r'\textrm{\textbf{JAM}}')
    labels.append(r'\textrm{\textbf{+SoLID}}')
    ax11.legend(handles,labels,loc=(0.60,0.00),fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.2,hspace=0.02,top=0.99,right=0.99)

    filename = 'plots/nobuo_sea_impact'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def ppdf_impact(WDIR,Q2=10):

    flavs = ['u','d','ub','db','g','sigma']
    
    nrows = 3
    ncols = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*10,nrows*5))
 
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            if j ==0: axs[i][j].set_xlim(0.01,0.73)
            if j ==1: axs[i][j].set_xlim(0.01,0.73)
            axs [i][j].axhline(0,lw = 1,color = 'k',alpha = 0.3,zorder=10)

    k = 0
    hand = {}
    for wdir in WDIR:
        filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
        data=load(filename)
        print('Loading %s'%filename)

        X = data['X']
        for flav in flavs:

            if   flav=='u':     i, j = 0, 0
            elif flav=='d':     i, j = 1, 0
            elif flav=='sigma': i, j = 2, 0
            elif flav=='g':     i, j = 2, 1
            elif flav=='ub':    i, j = 0, 1
            elif flav=='db':    i, j = 1, 1
            elif flav=='ub-db': i, j = 2, 1

            if k==0: alpha,color,zorder = 0.9, 'red' , 1.1
            if k==1: alpha,color,zorder = 0.7, 'cyan', 1.0

            if flav=='ub-db':
                result = np.array(data['XF']['ub']) - np.array(data['XF']['db'])
            elif flav=='sigma':
                result = np.array(data['XF']['up']) + np.array(data['XF']['dp']) + np.array(data['XF']['sp'])
            else:
                result = np.array(data['XF'][flav])

            mean = np.mean(result,axis=0)
            std  = np.std (result,axis=0)
           
            hand[k] = axs [i][j].fill_between(X,mean-std,mean+std,color=color,alpha=alpha,zorder=zorder)

        k+=1
        
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

    axs[2][0].set_ylim(-0.05,0.35)
    axs[2][0].yaxis.set_minor_locator(MultipleLocator(0.05)) 
    axs[2][0].set_yticks([0,0.1,0.2,0.3]) 
    axs[2][0].set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$']) 

    axs[2][1].set_ylim(0.00,0.16)
    axs[2][1].yaxis.set_minor_locator(MultipleLocator(0.01)) 
    axs[2][1].set_yticks([0,0.05,0.10,0.15]) 
    axs[2][1].set_yticklabels([r'$0$',r'$0.05$',r'$0.10$',r'$0.15$']) 
 
    ls = 40
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            axs [i][j].tick_params(axis='both', which='major', top=True, direction='in',labelsize=ls,length=10)
            axs [i][j].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=ls,length=5)
            axs[i][j].set_xticks([0.2,0.4,0.6])
            if i < 2: axs[i][j].set_xticklabels([])
            else: 
                axs[i][j].set_xticklabels([r'$0.2$',r'$0.4$',r'$0.6$'])
    
            if i == 2: axs[i][j].set_xlabel(r'\boldmath{$x$}',fontsize = 60)
    
            axs[i][j].xaxis.set_minor_locator(MultipleLocator(0.10)) 
       
    pdf_fs = 50
 
    axs[0][0].text( 0.77, 0.85, r'\boldmath$x \Delta u$',          transform=axs[0][0].transAxes,fontsize=pdf_fs)
    axs[1][0].text( 0.77, 0.05, r'\boldmath$x \Delta d$',          transform=axs[1][0].transAxes,fontsize=pdf_fs)
    axs[0][1].text( 0.77, 0.85, r'\boldmath$x \Delta \bar{u}$',    transform=axs[0][1].transAxes,fontsize=pdf_fs)
    axs[1][1].text( 0.77, 0.05, r'\boldmath$x \Delta \bar{d}$',    transform=axs[1][1].transAxes,fontsize=pdf_fs)
    axs[2][0].text( 0.77, 0.85, r'\boldmath$x \Delta \Sigma$',     transform=axs[2][0].transAxes,fontsize=pdf_fs)
    axs[2][1].text( 0.77, 0.85, r'\boldmath$x \Delta g$',          transform=axs[2][1].transAxes,fontsize=pdf_fs)

    axs[2][1].text(0.03, 0.05, r'$Q^2 = 10 ~ {\rm GeV}^2$',     transform=axs[2][1].transAxes,fontsize=35)
  

    handles, labels = [],[]
    handles.append(hand[1])
    handles.append(hand[0])
    labels.append(r'\textbf{\textbf{JAM}}')
    labels.append(r'\textbf{\textbf{+SoLID}}')
    axs[0][1].legend(handles,labels,fontsize = 45, loc=(0.52,0.35), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.02,wspace=0.20,top=0.99,right=0.99)

    filename = 'plots/nobuo_ppdfs' + ext 

    axs[0][0].set_rasterized(True)
    axs[0][1].set_rasterized(True)
    axs[1][0].set_rasterized(True)
    axs[1][1].set_rasterized(True)
    axs[2][0].set_rasterized(True)
    axs[2][1].set_rasterized(True)

    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

def plot_ppdfs_rel_uncert(WDIR, Q2=10):

    flavs = ['u','d','ub','db','ub-db','sp','g','sigma']
    
    nrows = 4
    ncols = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*7,nrows*5))
  
    axsL = {} 
    for i in range(nrows):
        axsL[i] = {}
        for j in range(ncols):
            divider = make_axes_locatable(axs[i][j])
            axsL[i][j] = divider.append_axes("right", size=3.25, pad=0.005, sharey=axs[i][j])
            axs[i][j].spines['right'].set_visible(False)
            axsL[i][j].spines['left'].set_visible(False)
            axsL[i][j].yaxis.set_ticks_position('right')
            if j ==0: axsL[i][j].set_xlim(0.1,0.73)
            if j ==1: axsL[i][j].set_xlim(0.1,0.53)
            axs[i][j].semilogx()
            axs[i][j].set_xlim(5e-3,0.1)
            if i==3 and j==1: continue
            axs [i][j].axhline(1,lw = 1,color = 'k',alpha = 0.3)
            axsL[i][j].axhline(1,lw = 1,color = 'k',alpha = 0.3)

    axs [3][1].axis("off")
    axsL[3][1].axis("off")

    filename = '%s/data/ppdf-Q2=%3.5f.dat'%(WDIR[0],Q2)
    data1=load(filename)
    print('Loading %s'%filename)
    filename = '%s/data/ppdf-Q2=%3.5f.dat'%(WDIR[1],Q2)
    data2=load(filename)
    print('Loading %s'%filename)

    hand = {}
    for flav in flavs:

        if   flav=='u':     i, j = 0, 0
        elif flav=='d':     i, j = 1, 0
        elif flav=='sigma': i, j = 2, 0
        elif flav=='g':     i, j = 3, 0
        elif flav=='ub':    i, j = 0, 1
        elif flav=='db':    i, j = 1, 1
        elif flav=='ub-db': i, j = 2, 1
        #elif flav=='sp':    i, j = 3, 1
        else: continue

        if flav=='ub-db':
            result1 = np.array(data1['XF']['ub']) - np.array(data1['XF']['db'])
            result2 = np.array(data2['XF']['ub']) - np.array(data2['XF']['db'])
        elif flav=='sigma':
            result1 = np.array(data1['XF']['up']) + np.array(data1['XF']['dp']) + np.array(data1['XF']['sp'])
            result2 = np.array(data2['XF']['up']) + np.array(data2['XF']['dp']) + np.array(data2['XF']['sp'])
        else:
            result1 = np.array(data1['XF'][flav])
            result2 = np.array(data2['XF'][flav])


        std1 = np.std(result1,axis=0)
        std2 = np.std(result2,axis=0)
  
        X = data1['X']
        ratio = std1/std2
        axs [i][j].plot(X,ratio,color='red',alpha=1.0,lw=3)
        axsL[i][j].plot(X,ratio,color='red',alpha=1.0,lw=3)

 
    ls = 30
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            axs [i][j].tick_params(axis='both', which='major', top=True, direction='in',labelsize=ls,length=10)
            axs [i][j].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=ls,length=5)
            axsL[i][j].tick_params(axis='both', which='major', top=True, right=True, labelright=False, direction='in',labelsize=ls,length=10)
            axsL[i][j].tick_params(axis='both', which='minor', top=True, right=True, labelright=False, direction='in',labelsize=ls,length=5)
            axs[i][j].set_xticks([0.01])
            if i < 3: axs[i][j].set_xticklabels([])
            else: axs[i][j].set_xticklabels([r'$0.01$'])
            if j==0: axsL[i][j].set_xticks([0.1,0.3,0.5,0.7])
            if j==1: axsL[i][j].set_xticks([0.1,0.3,0.5])
            if i < 3: axsL[i][j].set_xticklabels([])
            else: 
                if j==0: axsL[i][j].set_xticklabels([r'$0.1$',r'$0.3$',r'$0.5$',r'$0.7$'])
                if j==1: axsL[i][j].set_xticklabels([r'$0.1$',r'$0.3$',r'$0.5$'])
    
            if i == 3: axsL[i][j].set_xlabel(r'\boldmath{$x$}',fontsize = 45,loc = 'left')
   
            if j==1: axs[i][j].tick_params(labelleft=False) 

            axsL[i][j].xaxis.set_minor_locator(MultipleLocator(0.10)) 
      
            axs[i][j].set_ylim(0,1.49)
            axs[i][j].yaxis.set_minor_locator(MultipleLocator(0.10)) 
            axs[i][j].yaxis.set_major_locator(MultipleLocator(0.20)) 

 
    pdf_fs = 45
 
    axsL[0][0].text( 0.45, 0.85, r'\boldmath$x \Delta u$',          transform=axsL[0][0].transAxes,fontsize=pdf_fs)
    axs [1][0].text( 0.05, 0.05, r'\boldmath$x \Delta d$',          transform=axs[1][0].transAxes,fontsize=pdf_fs)
    axs [0][1].text( 0.10, 0.85, r'\boldmath$x \Delta \bar{u}$',    transform=axs[0][1].transAxes,fontsize=pdf_fs)
    axsL[1][1].text( 0.40, 0.05, r'\boldmath$x \Delta \bar{d}$',    transform=axsL[1][1].transAxes,fontsize=pdf_fs)
    axs [2][0].text( 0.10, 0.85, r'\boldmath$x \Delta \Sigma$',     transform=axs[2][0].transAxes,fontsize=pdf_fs)
    axsL[2][1].text(-1.00, 0.85, r'\boldmath$x (\Delta \bar{u} - \Delta \bar{d})$',    transform=axsL[2][1].transAxes,fontsize=pdf_fs)
    axs [3][0].text( 0.10, 0.05, r'\boldmath$x \Delta g$',          transform=axs[3][0].transAxes,fontsize=pdf_fs)
    #axs [3][1].text( 0.10, 0.85, r'\boldmath$x \Delta s^+$',        transform=axs[3][1].transAxes,fontsize=pdf_fs)

    axs[3][1].text(0.03, 0.05, r'$Q^2 = 10 ~ {\rm GeV}^2$',     transform=axs[3][1].transAxes,fontsize=35)

    #axs[0][0].text( 0.02, 0.70, r'\boldmath$\delta_{W^2 > 4}/\delta_{W^2 > 10}$',          transform=axs[0][0].transAxes,fontsize=45)
    axsL[0][0].text(-0.90, 0.15, r'\boldmath$\frac{\delta~({\rm +SoLID})}{\delta~({\rm JAM})}$',          transform=axsL[0][0].transAxes,fontsize=60)
  

    #handles, labels = [],[]
    #handles.append(hand[2])
    #handles.append(hand[1])
    #labels.append(r'\boldmath{$W^2 > 4$}  \rm \bf GeV\boldmath{$^2$}')
    #labels.append(r'\boldmath{$W^2 > 10$} \rm \bf GeV\boldmath{$^2$}')
    #axsL[0][0].legend(handles,labels,fontsize = 26, loc=(-1.10,0.72), frameon=0,handlelength=1,handletextpad=0.6)
  

    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.99,right=0.99,left=0.05)

    filename = 'plots/nobuo_ppdfs_rel_uncert' + ext 

    axs[0][0].set_rasterized(True)
    axs[0][1].set_rasterized(True)
    axs[1][0].set_rasterized(True)
    axs[1][1].set_rasterized(True)
    axs[2][0].set_rasterized(True)
    axs[2][1].set_rasterized(True)
    axs[3][0].set_rasterized(True)
    axs[3][1].set_rasterized(True)
    axsL[0][0].set_rasterized(True)
    axsL[0][1].set_rasterized(True)
    axsL[1][0].set_rasterized(True)
    axsL[1][1].set_rasterized(True)
    axsL[2][0].set_rasterized(True)
    axsL[2][1].set_rasterized(True)
    axsL[3][0].set_rasterized(True)
    axsL[3][1].set_rasterized(True)

    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

def polarization_impact(WDIR,Q2=10):

    M2 = 0.938**2 
    regen = False 
    W2min  = np.array([4,4])
    xmax   = Q2/(Q2+W2min-M2)
 
    ncols = 2
    nrows = 1
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*7,nrows*5))


    k = 0
    flavs = ['u','d']
    hand = {}

    colors = ['red','cyan']
    alpha = [0.9,0.9]
    zorders = [1.1,1.0]

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

            #result[:,X>xmax[k]] = 'nan'

            mean = np.mean(result,axis=0)
            std  = np.std (result,axis=0)

            if k==2:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color='none',edgecolor=colors[k],alpha=alpha[k],zorder=zorders[k],hatch='//')
            else:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color=colors[k],alpha=alpha[k],zorder=zorders[k])

        k+=1 
      
    ls = 25 

    for i in range(len(axs)):
        axs[i].set_xlim(0.01,0.80)
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
   
    #axs[1].set_ylim(-1.10,1.00)
    axs[1].set_ylim(-1.50,1.50)
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.10)) 
    axs[1].set_yticks([-1.0,-0.5,0,0.5,1.0]) 
    axs[1].set_yticklabels([r'$-1.0$',r'$-0.5$',r'$0$',r'$0.5$',r'1.0']) 
  

    axs[0].text(0.03, 0.80, r'\boldmath$\frac{\Delta u}{u}$', transform=axs[0].transAxes,fontsize=50)
    axs[1].text(0.03, 0.10, r'\boldmath$\frac{\Delta d}{d}$', transform=axs[1].transAxes,fontsize=50)
    
    axs[0].text(0.65, 0.05, r'$Q^2 = %d ~ {\rm GeV}^2$'%Q2,     transform=axs[0].transAxes,fontsize=20)


    handles,labels = [],[]
    handles.append(hand[1])
    handles.append(hand[0])
    labels.append(r'\textrm{\textbf{JAM}}')
    labels.append(r'\textrm{\textbf{+SoLID}}')
    axs[1].legend(handles,labels,fontsize = 25, loc=(0.00,0.70), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.01,wspace=0.15,top=0.97,right=0.99)

    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
 
    filename = 'plots/nobuo_polarization_impact' + ext 
    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

def plot_nuc_ratio(wdir1,wdir2,had='pi+'):

    load_config('%s/input.py'%wdir1)
    istep=core.get_istep()
    predictions1 = load('%s/data/predictions-%d.dat'%(wdir1,istep))
    predictions2 = load('%s/data/predictions-%d.dat'%(wdir2,istep))

    nrows,ncols=2,4
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax13 = py.subplot(nrows,ncols,3)
    ax14 = py.subplot(nrows,ncols,4)
    ax21 = py.subplot(nrows,ncols,5)
    ax22 = py.subplot(nrows,ncols,6)
    ax23 = py.subplot(nrows,ncols,7)
    ax24 = py.subplot(nrows,ncols,8)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][30000]='psidis/expdata/30000.xlsx'
    conf['datasets']['psidis']['xlsx'][30001]='psidis/expdata/30001.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    data1 = predictions1['reactions']['psidis']
    data2 = predictions2['reactions']['psidis']

    #--get theory by seperating solutions and taking mean
    for idx in tabs:
        predictions1 = copy.copy(data1[idx]['prediction-rep'])
        predictions2 = copy.copy(data2[idx]['prediction-rep'])
        del data1[idx]['prediction-rep']
        del data1[idx]['residuals-rep']
        del data1[idx]['shift-rep']
        del data1[idx]['rres-rep']
        del data1[idx]['r-residuals']
        del data1[idx]['n-residuals']
        del data2[idx]['prediction-rep']
        del data2[idx]['residuals-rep']
        del data2[idx]['shift-rep']
        del data2[idx]['rres-rep']
        del data2[idx]['r-residuals']
        del data2[idx]['n-residuals']
        data1[idx]['thy']  = np.mean(predictions1,axis=0)
        data1[idx]['dthy'] = np.std (predictions1,axis=0)
        data2[idx]['thy']  = np.mean(predictions2,axis=0)
        data2[idx]['dthy'] = np.std (predictions2,axis=0)


    if had=='pi+': idx = 30000
    if had=='pi-': idx = 30001

    Q2bins = [1.5,2.5,3.5,4.5,5.5,7]
    zbins  = [0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575]

    Q2bins = np.unique(data1[idx]['Q2'])
    zbins  = np.unique(data1[idx]['z'])

    #--set up bins
    DATA1 = {}
    DATA2 = {}
    for i in range(len(zbins)):
        DATA1[i] = {}
        DATA2[i] = {}
        for j in range(len(Q2bins)):
            query = "z==%s and Q2==%s"%(zbins[i],Q2bins[j])
            DATA1[i][j] = pd.DataFrame(data1[idx]).query(query)
            DATA2[i][j] = pd.DataFrame(data2[idx]).query(query)
    
    #######################
    #--plot (abs(nuc - eff)/err)
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in tabs:
        for i in range(len(zbins)):
            for j in range(len(Q2bins)):
                tab1 = DATA1[i][j]
                tab2 = DATA2[i][j]
                X = tab1['X']
                values = tab1['value']
                alpha  = tab1['alpha']
                if i==0: ax = ax11
                if i==1: ax = ax12
                if i==2: ax = ax13
                if i==3: ax = ax14
                if i==4: ax = ax21
                if i==5: ax = ax22
                if i==6: ax = ax23
                if i==7: ax = ax24
                if j==0: color,fmt,ms = 'firebrick','o',20
                if j==1: color,fmt,ms = 'darkgreen','^',20
                if j==2: color,fmt,ms = 'black'    ,'*',20
                if j==3: color,fmt,ms = 'cyan'     ,'s',20
                if j==4: color,fmt,ms = 'magenta'  ,'v',20
                if j==5: color,fmt,ms = 'orange'   ,'D',20


                thy1 = tab1['thy']
                thy2 = tab2['thy']

                ratio = np.abs(thy1-thy2)/alpha
                thy_plot[j] ,= ax.plot(X,ratio,color=color,ls=':',alpha=0.5)
                #thy_band[j]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

                hand[j] = ax.scatter(X,ratio,color=color,marker=fmt,s=ms)

                if j==0: ax.text(0.02, 0.90, r'\boldmath$z = %3.3f$'%(zbins[i]),transform=ax.transAxes,size=30)
  
    for ax in [ax11,ax12,ax13,ax14,ax21,ax22,ax23,ax24]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_xlim(0.03,0.57)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(1,0,1,color='black',alpha=0.2)

    for ax in [ax12,ax13,ax14,ax22,ax23,ax24]:
        ax.tick_params(labelleft=False)

    for ax in [ax11,ax12,ax13,ax14]:
        ax.tick_params(labelbottom=False)

    for ax in [ax21,ax22,ax23,ax24]:
        ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=40)

    for ax in [ax11,ax12,ax13,ax14,ax21,ax22,ax23,ax24]:
        ax.semilogy()
        ax.set_ylim(1e-1,50)
        #ax.set_yticks([0,0.2,0.4,0.6,0.8])
        #ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
        #ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    #for ax in [ax12,ax22]:
    #    ax.tick_params(axis='both',which='both',labelleft=False)



    #if had=='pi+': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^+,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    #if had=='pi-': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^-,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    if had=='pi+': ax11.text(0.80, 0.05, r'\boldmath$\pi^+$',transform=ax11.transAxes,size=60)
    if had=='pi-': ax11.text(0.80, 0.05, r'\boldmath$\pi^-$',transform=ax11.transAxes,size=60)

    ax11.set_ylabel(r'\boldmath$| {\rm A_{LL}^{^3{\rm He}}(KPSV) - A_{LL}^{^3{\rm He}}(SS)} |/\delta A_{LL}^{^3{\rm He}}$',size=30)
    ax11.yaxis.label.set_position((-0.2,0.0))

    fs = 25
    x,y = 0.70, 0.00
    handles, labels = [],[]
    for j in [0,1]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax12.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)


    handles, labels = [],[]
    for j in [2,3]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax13.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)

    handles, labels = [],[]
    for j in [4,5]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])

    ax14.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)
    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/nobuo_nuc_ratio_%s'%(cwd,had)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving SoLID PSIDIS %s plot to %s'%(had,filename))

def plot_nuc_ratio2(wdir1,wdir2,had='pi+'):

    load_config('%s/input.py'%wdir1)
    istep=core.get_istep()
    predictions1 = load('%s/data/predictions-%d.dat'%(wdir1,istep))
    predictions2 = load('%s/data/predictions-%d.dat'%(wdir2,istep))

    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax13 = py.subplot(nrows,ncols,3)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][30000]='psidis/expdata/30000.xlsx'
    conf['datasets']['psidis']['xlsx'][30001]='psidis/expdata/30001.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    data1 = predictions1['reactions']['psidis']
    data2 = predictions2['reactions']['psidis']

    #--get theory by seperating solutions and taking mean
    for idx in tabs:
        predictions1 = copy.copy(data1[idx]['prediction-rep'])
        predictions2 = copy.copy(data2[idx]['prediction-rep'])
        del data1[idx]['prediction-rep']
        del data1[idx]['residuals-rep']
        del data1[idx]['shift-rep']
        del data1[idx]['rres-rep']
        del data1[idx]['r-residuals']
        del data1[idx]['n-residuals']
        del data2[idx]['prediction-rep']
        del data2[idx]['residuals-rep']
        del data2[idx]['shift-rep']
        del data2[idx]['rres-rep']
        del data2[idx]['r-residuals']
        del data2[idx]['n-residuals']
        data1[idx]['thy']  = np.mean(predictions1,axis=0)
        data1[idx]['dthy'] = np.std (predictions1,axis=0)
        data2[idx]['thy']  = np.mean(predictions2,axis=0)
        data2[idx]['dthy'] = np.std (predictions2,axis=0)


    if had=='pi+': idx = 30000
    if had=='pi-': idx = 30001

    Q2bins = [1.5,2.5,3.5,4.5,5.5,7]
    zbins  = [0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575]

    Q2bins = np.unique(data1[idx]['Q2'])
    zbins  = np.unique(data1[idx]['z'])

    #--set up bins
    DATA1 = {}
    DATA2 = {}
    for i in range(len(zbins)):
        DATA1[i] = {}
        DATA2[i] = {}
        for j in range(len(Q2bins)):
            query = "z==%s and Q2==%s"%(zbins[i],Q2bins[j])
            DATA1[i][j] = pd.DataFrame(data1[idx]).query(query)
            DATA2[i][j] = pd.DataFrame(data2[idx]).query(query)
    
    #######################
    #--plot (abs(nuc - eff)/err)
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in tabs:
        cnt = 0
        for i in range(len(zbins)):
            if i not in [0,2,6]: continue
            for j in range(len(Q2bins)):
                tab1 = DATA1[cnt][j]
                tab2 = DATA2[cnt][j]
                X = tab1['X']
                values = tab1['value']
                alpha  = tab1['alpha']
                if cnt==0: ax = ax11
                if cnt==1: ax = ax12
                if cnt==2: ax = ax13
                if j==0: color,fmt,ms = 'firebrick','o',20
                if j==1: color,fmt,ms = 'darkgreen','^',20
                if j==2: color,fmt,ms = 'black'    ,'*',20
                if j==3: color,fmt,ms = 'cyan'     ,'s',20
                if j==4: color,fmt,ms = 'magenta'  ,'v',20
                if j==5: color,fmt,ms = 'orange'   ,'D',20


                thy1 = tab1['thy']
                thy2 = tab2['thy']

                ratio = np.abs(thy1-thy2)/alpha
                thy_plot[j] ,= ax.plot(X,ratio,color=color,ls=':',alpha=0.5)
                #thy_band[j]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

                hand[j] = ax.scatter(X,ratio,color=color,marker=fmt,s=ms)

                if j==0: ax.text(0.02, 0.87, r'\boldmath$z = %3.3f$'%(zbins[i]),transform=ax.transAxes,size=30)
            cnt += 1
  
    for ax in [ax11,ax12,ax13]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=25)
        ax.set_xlim(0.03,0.57)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(1,0,1,color='black',alpha=0.2)

    for ax in [ax12,ax13]:
        ax.tick_params(labelleft=False)

    for ax in [ax11,ax12,ax13]:
        ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=35)

    for ax in [ax11,ax12,ax13]:
        ax.semilogy()
        ax.set_ylim(2e-2,50)
        #ax.set_yticks([0,0.2,0.4,0.6,0.8])
        #ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
        #ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    #for ax in [ax12,ax22]:
    #    ax.tick_params(axis='both',which='both',labelleft=False)



    #if had=='pi+': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^+,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    #if had=='pi-': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^-,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    if had=='pi+': ax11.text(0.15, 0.05, r'\boldmath$\pi^+$',transform=ax11.transAxes,size=60)
    if had=='pi-': ax11.text(0.15, 0.05, r'\boldmath$\pi^-$',transform=ax11.transAxes,size=60)

    ax11.set_ylabel(r'\boldmath$| {\rm A_{LL}^{^3{\rm He}}(KPSV) - A_{LL}^{^3{\rm He}}(SS)} |/\delta A_{LL}^{^3{\rm He}}$',size=16)
    ax11.yaxis.label.set_position((-0.2,0.35))

    fs = 25
    x,y = 0.70, 0.00
    handles, labels = [],[]
    for j in [0,1]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax11.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)


    handles, labels = [],[]
    for j in [2,3]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax12.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)

    handles, labels = [],[]
    for j in [4,5]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])

    ax13.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)
    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/nobuo_nuc_ratio2_%s'%(cwd,had)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving SoLID PSIDIS %s plot to %s'%(had,filename))

def plot_nuc_ratio3(wdir1,wdir2):

    load_config('%s/input.py'%wdir1)
    istep=core.get_istep()
    predictions1 = load('%s/data/predictions-%d.dat'%(wdir1,istep))
    predictions2 = load('%s/data/predictions-%d.dat'%(wdir2,istep))

    nrows,ncols=1,2
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][30000]='psidis/expdata/30000.xlsx'
    conf['datasets']['psidis']['xlsx'][30001]='psidis/expdata/30001.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    data1 = predictions1['reactions']['psidis']
    data2 = predictions2['reactions']['psidis']

    #--get theory by seperating solutions and taking mean
    for idx in tabs:
        predictions1 = copy.copy(data1[idx]['prediction-rep'])
        predictions2 = copy.copy(data2[idx]['prediction-rep'])
        del data1[idx]['prediction-rep']
        del data1[idx]['residuals-rep']
        del data1[idx]['shift-rep']
        del data1[idx]['rres-rep']
        del data1[idx]['r-residuals']
        del data1[idx]['n-residuals']
        del data2[idx]['prediction-rep']
        del data2[idx]['residuals-rep']
        del data2[idx]['shift-rep']
        del data2[idx]['rres-rep']
        del data2[idx]['r-residuals']
        del data2[idx]['n-residuals']
        data1[idx]['thy']  = np.mean(predictions1,axis=0)
        data1[idx]['dthy'] = np.std (predictions1,axis=0)
        data2[idx]['thy']  = np.mean(predictions2,axis=0)
        data2[idx]['dthy'] = np.std (predictions2,axis=0)


    idxs = [30000, 30001]

    Q2bins = [1.5,2.5,3.5,4.5,5.5,7]
    zbins  = [0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575]


    #--set up bins
    DATA1 = {}
    DATA2 = {}
    for idx in idxs:
        Q2bins = np.unique(data1[idx]['Q2'])
        zbins  = np.unique(data1[idx]['z'])
        DATA1[idx] = {}
        DATA2[idx] = {}
        for i in range(len(zbins)):
            DATA1[idx][i] = {}
            DATA2[idx][i] = {}
            for j in range(len(Q2bins)):
                query = "z==%s and Q2==%s"%(zbins[i],Q2bins[j])
                DATA1[idx][i][j] = pd.DataFrame(data1[idx]).query(query)
                DATA2[idx][i][j] = pd.DataFrame(data2[idx]).query(query)
    
    #######################
    #--plot (abs(nuc - eff)/err)
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in tabs:
        cnt = 0
        for i in range(len(zbins)):
            if i not in [0]: continue
            for j in range(len(Q2bins)):
                tab1 = DATA1[idx][cnt][j]
                tab2 = DATA2[idx][cnt][j]
                X = tab1['X']
                values = tab1['value']
                alpha  = tab1['alpha']
                if idx==30000: ax = ax11
                if idx==30001: ax = ax12
                if j==0: color,fmt,ms = 'firebrick','o',20
                if j==1: color,fmt,ms = 'darkgreen','^',20
                if j==2: color,fmt,ms = 'black'    ,'*',20
                if j==3: color,fmt,ms = 'cyan'     ,'s',20
                if j==4: color,fmt,ms = 'magenta'  ,'v',20
                if j==5: color,fmt,ms = 'orange'   ,'D',20


                thy1 = tab1['thy']
                thy2 = tab2['thy']

                ratio = np.abs(thy1-thy2)/0.03#/alpha
                thy_plot[j] ,= ax.plot(X,ratio,color=color,ls=':',alpha=0.5)
                #thy_band[j]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

                hand[j] = ax.scatter(X,ratio,color=color,marker=fmt,s=ms)

                if idx==30000 and j==0: ax.text(0.02, 0.87, r'\boldmath$z = %3.3f$'%(zbins[i]),transform=ax.transAxes,size=30)
            cnt += 1
  
    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        ax.set_xlim(0.03,0.57)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(1,0,1,color='black',alpha=0.2)

    for ax in [ax12]:
        ax.tick_params(labelleft=False)

    for ax in [ax11,ax12]:
        ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=35)

    for ax in [ax11,ax12]:
        ax.semilogy()
        #ax.set_ylim(2e-2,50)
        ax.set_ylim(1e-4,9e-2)
        #ax.set_yticks([0,0.2,0.4,0.6,0.8])
        #ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
        #ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    #for ax in [ax12,ax22]:
    #    ax.tick_params(axis='both',which='both',labelleft=False)



    #if had=='pi+': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^+,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    #if had=='pi-': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^-,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    ax11.text(0.02, 0.05, r'\boldmath$\pi^+$',transform=ax11.transAxes,size=60)
    ax12.text(0.02, 0.05, r'\boldmath$\pi^-$',transform=ax12.transAxes,size=60)

    ax11.set_ylabel(r'\boldmath$| {\rm A_{LL}^{^3{\rm He}}({\rm eff}) - A_{LL}^{^3{\rm He}}({\rm nuc})} |/\langle A \rangle$',size=20)
    ax11.yaxis.label.set_position((-0.2,0.40))

    fs = 20
    x,y = 0.70, 0.00
    handles, labels = [],[]
    for j in [0,1,2]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax11.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)


    handles, labels = [],[]
    for j in [3,4,5]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax12.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)

    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/nobuo_nuc_ratio3'%(cwd)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving SoLID PSIDIS plot to %s'%(filename))

def plot_nuc_ratio4(wdir1,wdir2):

    load_config('%s/input.py'%wdir1)
    istep=core.get_istep()
    predictions1 = load('%s/data/predictions-%d.dat'%(wdir1,istep))
    predictions2 = load('%s/data/predictions-%d.dat'%(wdir2,istep))

    nrows,ncols=1,2
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][30000]='psidis/expdata/30000.xlsx'
    conf['datasets']['psidis']['xlsx'][30001]='psidis/expdata/30001.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    data1 = predictions1['reactions']['psidis']
    data2 = predictions2['reactions']['psidis']

    #--get theory by seperating solutions and taking mean
    for idx in tabs:
        predictions1 = copy.copy(data1[idx]['prediction-rep'])
        predictions2 = copy.copy(data2[idx]['prediction-rep'])
        del data1[idx]['prediction-rep']
        del data1[idx]['residuals-rep']
        del data1[idx]['shift-rep']
        del data1[idx]['rres-rep']
        del data1[idx]['r-residuals']
        del data1[idx]['n-residuals']
        del data2[idx]['prediction-rep']
        del data2[idx]['residuals-rep']
        del data2[idx]['shift-rep']
        del data2[idx]['rres-rep']
        del data2[idx]['r-residuals']
        del data2[idx]['n-residuals']
        data1[idx]['thy']  = np.mean(predictions1,axis=0)
        data1[idx]['dthy'] = np.std (predictions1,axis=0)
        data2[idx]['thy']  = np.mean(predictions2,axis=0)
        data2[idx]['dthy'] = np.std (predictions2,axis=0)


    idxs = [30000, 30001]

    Q2bins = [1.5,2.5,3.5,4.5,5.5,7]
    zbins  = [0.225,0.275,0.325,0.375,0.425,0.475,0.525,0.575]


    #--set up bins
    DATA1 = {}
    DATA2 = {}
    for idx in idxs:
        Q2bins = np.unique(data1[idx]['Q2'])
        zbins  = np.unique(data1[idx]['z'])
        DATA1[idx] = {}
        DATA2[idx] = {}
        for i in range(len(zbins)):
            DATA1[idx][i] = {}
            DATA2[idx][i] = {}
            for j in range(len(Q2bins)):
                query = "z==%s and Q2==%s"%(zbins[i],Q2bins[j])
                DATA1[idx][i][j] = pd.DataFrame(data1[idx]).query(query)
                DATA2[idx][i][j] = pd.DataFrame(data2[idx]).query(query)
    
    #######################
    #--plot (abs(nuc - eff)/err)
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in tabs:
        cnt = 0
        for i in range(len(zbins)):
            if i not in [0]: continue
            for j in range(len(Q2bins)):
                tab1 = DATA1[idx][cnt][j]
                tab2 = DATA2[idx][cnt][j]
                X = tab1['X']
                values = tab1['value']
                alpha  = tab1['alpha']
                if idx==30000: ax = ax11
                if idx==30001: ax = ax12
                if j==0: color,fmt,ms = 'firebrick','o',20
                if j==1: color,fmt,ms = 'darkgreen','^',20
                if j==2: color,fmt,ms = 'black'    ,'*',20
                if j==3: color,fmt,ms = 'cyan'     ,'s',20
                if j==4: color,fmt,ms = 'magenta'  ,'v',20
                if j==5: color,fmt,ms = 'orange'   ,'D',20


                thy1 = tab1['thy']
                thy2 = tab2['thy']

                ratio = np.abs(thy1-thy2)/alpha
                thy_plot[j] ,= ax.plot(X,ratio,color=color,ls=':',alpha=0.5)
                #thy_band[j]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

                hand[j] = ax.scatter(X,ratio,color=color,marker=fmt,s=ms)

                if idx==30000 and j==0: ax.text(0.02, 0.87, r'\boldmath$z = %3.3f$'%(zbins[i]),transform=ax.transAxes,size=30)
            cnt += 1
  
    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=20)
        ax.set_xlim(0.03,0.57)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(1,0,1,color='black',alpha=0.2)

    for ax in [ax12]:
        ax.tick_params(labelleft=False)

    for ax in [ax11,ax12]:
        ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=35)

    for ax in [ax11,ax12]:
        ax.semilogy()
        ax.set_ylim(6e-3,20)
        #ax.set_ylim(1e-4,9e-2)
        #ax.set_yticks([0,0.2,0.4,0.6,0.8])
        #ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
        #ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    #for ax in [ax12,ax22]:
    #    ax.tick_params(axis='both',which='both',labelleft=False)



    #if had=='pi+': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^+,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    #if had=='pi-': ax11.text(0.60, 0.15, r'\boldmath$A^{\pi^-,^3{\rm He}}_{LL}$',transform=ax11.transAxes,size=50)
    ax11.text(0.02, 0.05, r'\boldmath$\pi^+$',transform=ax11.transAxes,size=60)
    ax12.text(0.02, 0.05, r'\boldmath$\pi^-$',transform=ax12.transAxes,size=60)

    ax11.set_ylabel(r'\boldmath$| {\rm A_{LL}^{^3{\rm He}}({\rm eff}) - A_{LL}^{^3{\rm He}}({\rm nuc})} |/\Delta A$',size=20)
    ax11.yaxis.label.set_position((-0.2,0.40))

    fs = 20
    x,y = 0.70, 0.00
    handles, labels = [],[]
    for j in [0,1,2]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax11.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)


    handles, labels = [],[]
    for j in [3,4,5]:
        #handles.append((thy_plot[j],thy_band[j],hand[j]))
        handles.append((thy_plot[j],hand[j]))
        labels.append(r'$Q^2 = %2.1f$'%Q2bins[j])
    
    ax12.legend(handles,labels,frameon=False,fontsize=fs,loc=(x,y),handletextpad = 0.5, handlelength = 1.0,ncol=1,columnspacing=0.5)

    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/nobuo_nuc_ratio4'%(cwd)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving SoLID PSIDIS plot to %s'%(filename))

if __name__ == "__main__":

    wdir1 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/nobuo/fixedFF'#/pos_g'
    wdir2 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/star_redo/pos_g'

    WDIR = [wdir1,wdir2]

    #plot_ppdfs_rel_uncert(WDIR, Q2=10)
    #ppdf_impact(WDIR, Q2=10)
    #asym_impact(WDIR, Q2=10)
    #sea_impact(WDIR, Q2=10)
    #polarization_impact(WDIR,Q2=10)

    wdir =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/nobuo/fixedFF'

    #plot_solid(wdir, kc, had='pi+')
    #plot_solid(wdir, kc, had='pi-')

    wdir1 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/nobuo/fixedFF_KPSV'
    wdir2 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/nobuo/fixedFF_SS'

    wdir1 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/nobuo/fixedFF'
    wdir2 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/nobuo/fixedFF_KPSV'

    #plot_nuc_ratio3(wdir1,wdir2)
    plot_nuc_ratio4(wdir1,wdir2)


