#!/usr/bin/env python
import sys, os
import numpy as np
import copy

path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

os.environ["LHAPDF_DATA_PATH"] = '/work/JAM/ccocuzza/lhapdf/python3/sets'

import lhapdf
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib import cm
import pylab as py

from qcdlib.qpdcalc import QPDCALC

#--from scipy stack 
from scipy.integrate import quad
from scipy.integrate import cumulative_trapezoid

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
# from analysis.corelib import classifier

## from local
from analysis.qpdlib import pdf

import kmeanconf as kc

#regen = True
regen = False

#ext = '.png'
ext = '.pdf'

flavors = []
flavors.append('up')
flavors.append('dp')
flavors.append('sp')
flavors.append('g')
flavors.append('ub')
flavors.append('db')

# cmap = matplotlib.cm.get_cmap('plasma')

def plot_ppdfs_with_strange(wdir,Q2=10,SETS={}):

    nrows,ncols=3,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    axs,axLs = {},{}
    for i in range(N):
        axs[i+1] = py.subplot(nrows,ncols,i+1)
        divider = make_axes_locatable(axs[i+1])
        axLs[i+1] = divider.append_axes("right",size=2.75,pad=0,sharey=axs[i+1])
        axLs[i+1].set_xlim(0.1,0.9)
        axLs[i+1].spines['left'].set_visible(False)
        axLs[i+1].yaxis.set_ticks_position('right')
        py.setp(axLs[i+1].get_xticklabels(),visible=True)

        axs[i+1].spines['right'].set_visible(False)

    thy  = {}
    hand = {}
    load_config('%s/input.py'%wdir)

    filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
    if regen: pdf.gen_xf(wdir,Q2)
    try:
        data=load(filename)
    except:
        pdf.gen_xf(wdir,Q2)
        data=load(filename)
    print('Loading %s'%filename)

    replicas=core.get_replicas(wdir)
    names    = core.get_replicas_names(wdir)

    X=data['X']

    for flav in data['XF']:
        mean = np.mean(data['XF'][flav],axis=0)
        std  = np.std (data['XF'][flav],axis=0)

        if flav=='up'      : ax,axL,zorder = axs[1],axLs[1], 1.4 
        elif flav=='ub'    : ax,axL,zorder = axs[2],axLs[2], 1.4 
        elif flav=='dp'    : ax,axL,zorder = axs[3],axLs[3], 1.4 
        elif flav=='db'    : ax,axL,zorder = axs[4],axLs[4], 1.4 
        elif flav=='g'     : ax,axL,zorder = axs[5],axLs[5], 1.4 
        elif flav=='sp'    : ax,axL,zorder = axs[6],axLs[6], 1
        
        else: continue


        #--plot each replica
        #--plot average and standard deviation
        thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=zorder)
        axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=zorder)

        #--plot other PPDF sets
        for SET in SETS:
            _SET, _color, _alpha, _zorder = SETS[SET][0], SETS[SET][1], SETS[SET][2], SETS[SET][4]

            if SET=='DSSV08':
                _X=10**np.linspace(-4,-1,200)
                _X=np.append(_X,np.linspace(0.1,0.99,200))
                ub = _SET.xfxQ2(-2,_X,Q2)
                db = _SET.xfxQ2(-1,_X,Q2)
                if flav=='up'   : ppdf = _SET.xfxQ2(-2 ,_X,Q2) + _SET.xfxQ2( 2,_X,Q2) 
                elif flav=='dp' : ppdf = _SET.xfxQ2(-1 ,_X,Q2) + _SET.xfxQ2( 1,_X,Q2) 
                elif flav=='g'  : ppdf = _SET.xfxQ2( 21,_X,Q2)
                elif flav=='sp' : ppdf = _SET.xfxQ2(-3 ,_X,Q2) + _SET.xfxQ2( 3,_X,Q2) 
                elif flav=='ub' : ppdf = _SET.xfxQ2(-2 ,_X,Q2)
                elif flav=='db' : ppdf = _SET.xfxQ2(-1 ,_X,Q2)
                mean = ppdf[0]
                std  = 0
                for i in range(1,20):
                    std += (ppdf[i] - ppdf[-i])**2
                std = np.sqrt(std)/2.0
                hand[SET] = ax.fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=_zorder)
                axL.           fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=_zorder)

            if SET=='NNPDF' or SET=='JAM17' or SET=='DSSV14':
                ppdf = _SET.get_xpdf(flav,X,Q2)
                hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=_zorder)
                axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=_zorder)



    for i in range(N):
          axs[i+1].set_xlim(8e-3,0.1)
          axs[i+1].semilogx()

          axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
          axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
          axs[i+1].set_xticks([0.01,0.1])
          axs[i+1].set_xticklabels([r'$0.01$',r'$0.1$'])

          if i in [0,2,4]: axLs[i+1].set_xlim(0.1,0.71)
          if i in [1,3,5]: axLs[i+1].set_xlim(0.1,0.51)

          axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=8)
          axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=4)
          if i in [0,2,4]: axLs[i+1].set_xticks([0.1,0.3,0.5,0.7])
          if i in [0,2,4]: axLs[i+1].set_xticklabels([r'',r'$0.3$',r'$0.5$',r'$0.7$'])
          if i in [1,3,5]: axLs[i+1].set_xticks([0.1,0.3,0.5])
          if i in [1,3,5]: axLs[i+1].set_xticklabels([r'',r'$0.3$',r'$0.5$'])

          axLs[i+1].xaxis.set_minor_locator(MultipleLocator(0.10)) 

    for i in [1,2,3,4]:
          axs[i] .tick_params(labelbottom=False)
          axLs[i].tick_params(labelbottom=False)

    

    axs[1].set_ylim(-0.01,0.39)
    axs[1].set_yticks([0,0.1,0.2,0.3])
    axs[1].set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$'])
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.05)) 

    axs[2].set_ylim(-0.025,0.035)
    axs[2].set_yticks([-0.02,-0.01,0,0.01,0.02,0.03])
    axs[2].set_yticklabels([r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$',r'$0.03$'])
    axs[2].yaxis.set_minor_locator(MultipleLocator(0.005)) 


    axs[3].set_ylim(-0.18,0.04)
    axs[3].set_yticks([-0.15,-0.10,-0.05,0])
    axs[3].set_yticklabels([r'$-0.15$',r'$-0.10$',r'$-0.05$',r'$0$'])
    axs[3].yaxis.set_minor_locator(MultipleLocator(0.01)) 


    axs[4].set_ylim(-0.035,0.03)
    axs[4].set_yticks([-0.03,-0.02,-0.01,0,0.01,0.02])
    axs[4].set_yticklabels([r'$-0.03$',r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$'])
    axs[4].yaxis.set_minor_locator(MultipleLocator(0.01)) 


    axs[5].set_ylim(-0.15,0.23)  
    axs[5].set_yticks([-0.1,0.0,0.1,0.2])
    axs[5].set_yticklabels([r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])
    axs[5].yaxis.set_minor_locator(MultipleLocator(0.05)) 


    axs[6].set_ylim(-0.10,0.10)
    axs[6].set_yticks([-0.08,-0.04,0.00,0.04,0.08])
    axs[6].set_yticklabels([r'$-0.08$',r'$-0.04$',r'$0$',r'$0.04$',r'$0.08$'])
    axs[6].yaxis.set_minor_locator(MultipleLocator(0.02)) 

    for i in range(N):
        axs [i+1].axhline(0  ,color='k',alpha=0.5,zorder=6)
        axLs[i+1].axhline(0  ,color='k',alpha=0.5,zorder=6)
        # axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
        # axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

    axLs[5].set_xlabel(r'\boldmath$x$',size=45,loc = 'left')
    axLs[6].set_xlabel(r'\boldmath$x$',size=45,loc = 'left')   
    # axLs[5].xaxis.set_label_coords(0.95,0.00)
    # axLs[6].xaxis.set_label_coords(0.95,0.00)
      

    axs [1].text(0.13,0.85,r'\boldmath{$x \Delta u^+$}',                        transform=axs[1] .transAxes,size=45)
    axLs[2].text(0.35,0.85,r'\boldmath{$x \Delta \bar{u}$}',                    transform=axLs[2].transAxes,size=45)
    axs [3].text(0.13,0.85,r'\boldmath{$x \Delta d^+$}',                        transform=axs[3] .transAxes,size=45)
    axs [4].text(0.13,0.85,r'\boldmath{$x \Delta \bar{d}$}',                    transform=axs[4] .transAxes,size=45)
    axLs[5].text(0.35,0.85,r'\boldmath{$x \Delta g$}'  ,                        transform=axLs[5].transAxes,size=45)
    axLs[6].text(0.35,0.85,r'\boldmath{$x \Delta s^+$}',                        transform=axLs[6].transAxes,size=45)


    axLs[5].text(0.10,0.05,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axLs[5].transAxes,size=25)

    handles,labels = [thy],[r'\rm \bf JAM25']
    handles.append(hand['JAM17'])
    labels.append(SETS['JAM17'][3])
    #axLs[3].legend(handles,labels,loc=(0.05,0.68),fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
    axLs[1].legend(handles,labels,loc=(-1.05,0.50),fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
    
    handles = []
    labels = []
    handles.append(hand['NNPDF'])
    labels.append(SETS['NNPDF'][3])
    handles.append(hand['DSSV14'])
    labels.append(SETS['DSSV14'][3])
    #axLs[4].legend(handles,labels,loc='upper right',fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
    axLs[3].legend(handles,labels,loc=(-1.05,0.00),fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
      
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.2,hspace=0.02)

    # filename = '%s/gallery/ppdfs-Q2=%3.5f'%(wdir,Q2)
    filename = 'plots/fig_ppdfs_compare'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def plot_ppdfs_backup(wdir,Q2=10,SETS={}):

    nrows,ncols=3,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    axs,axLs = {},{}
    for i in range(6):
        axs[i+1] = py.subplot(nrows,ncols,i+1)
        divider = make_axes_locatable(axs[i+1])
        axLs[i+1] = divider.append_axes("right",size=2.75,pad=0,sharey=axs[i+1])
        axLs[i+1].set_xlim(0.1,0.9)
        axLs[i+1].spines['left'].set_visible(False)
        axLs[i+1].yaxis.set_ticks_position('right')
        py.setp(axLs[i+1].get_xticklabels(),visible=True)

        axs[i+1].spines['right'].set_visible(False)

    axs [6].axis("off")
    axLs[6].axis("off")

    thy  = {}
    hand = {}
    load_config('%s/input.py'%wdir)

    filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
    if regen: pdf.gen_xf(wdir,Q2)
    try:
        data=load(filename)
    except:
        pdf.gen_xf(wdir,Q2)
        data=load(filename)
    print('Loading %s'%filename)

    replicas=core.get_replicas(wdir)
    names    = core.get_replicas_names(wdir)

    X=data['X']

    for flav in data['XF']:
        mean = np.mean(data['XF'][flav],axis=0)
        std  = np.std (data['XF'][flav],axis=0)

        if flav=='up'      : ax,axL,zorder = axs[1],axLs[1], 1.4 
        elif flav=='ub'    : ax,axL,zorder = axs[2],axLs[2], 1.4 
        elif flav=='dp'    : ax,axL,zorder = axs[3],axLs[3], 1.4 
        elif flav=='db'    : ax,axL,zorder = axs[4],axLs[4], 1.4 
        elif flav=='g'     : ax,axL,zorder = axs[5],axLs[5], 1.4 
        
        else: continue


        #--plot each replica
        #--plot average and standard deviation
        thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=zorder)
        axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=zorder)

        #--plot other PPDF sets
        for SET in SETS:
            _SET, _color, _alpha, _zorder = SETS[SET][0], SETS[SET][1], SETS[SET][2], SETS[SET][4]

            if SET=='DSSV08':
                _X=10**np.linspace(-4,-1,200)
                _X=np.append(_X,np.linspace(0.1,0.99,200))
                ub = _SET.xfxQ2(-2,_X,Q2)
                db = _SET.xfxQ2(-1,_X,Q2)
                if flav=='up'   : ppdf = _SET.xfxQ2(-2 ,_X,Q2) + _SET.xfxQ2( 2,_X,Q2) 
                elif flav=='dp' : ppdf = _SET.xfxQ2(-1 ,_X,Q2) + _SET.xfxQ2( 1,_X,Q2) 
                elif flav=='g'  : ppdf = _SET.xfxQ2( 21,_X,Q2)
                elif flav=='ub' : ppdf = _SET.xfxQ2(-2 ,_X,Q2)
                elif flav=='db' : ppdf = _SET.xfxQ2(-1 ,_X,Q2)
                mean = ppdf[0]
                std  = 0
                for i in range(1,20):
                    std += (ppdf[i] - ppdf[-i])**2
                std = np.sqrt(std)/2.0
                hand[SET] = ax.fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=_zorder)
                axL.           fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=_zorder)

            if SET=='NNPDF' or SET=='JAM17' or SET=='DSSV14':
                ppdf = _SET.get_xpdf(flav,X,Q2)
                hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=_zorder)
                axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=_zorder)



    for i in range(5):
          axs[i+1].set_xlim(8e-3,0.1)
          axs[i+1].semilogx()

          axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
          axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
          axs[i+1].set_xticks([0.01,0.1])
          axs[i+1].set_xticklabels([r'$0.01$',r'$0.1$'])

          if i in [0,2,4]: axLs[i+1].set_xlim(0.1,0.71)
          if i in [1,3,5]: axLs[i+1].set_xlim(0.1,0.51)

          axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=8)
          axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=4)
          if i in [0,2,4]: axLs[i+1].set_xticks([0.1,0.3,0.5,0.7])
          if i in [0,2,4]: axLs[i+1].set_xticklabels([r'',r'$0.3$',r'$0.5$',r'$0.7$'])
          if i in [1,3,5]: axLs[i+1].set_xticks([0.1,0.3,0.5])
          if i in [1,3,5]: axLs[i+1].set_xticklabels([r'',r'$0.3$',r'$0.5$'])

          axLs[i+1].xaxis.set_minor_locator(MultipleLocator(0.10)) 

    for i in [1,2,3]:
          axs[i] .tick_params(labelbottom=False)
          axLs[i].tick_params(labelbottom=False)

    

    axs[1].set_ylim(-0.01,0.39)
    axs[1].set_yticks([0,0.1,0.2,0.3])
    axs[1].set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$'])
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.05)) 

    axs[2].set_ylim(-0.025,0.035)
    axs[2].set_yticks([-0.02,-0.01,0,0.01,0.02,0.03])
    axs[2].set_yticklabels([r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$',r'$0.03$'])
    axs[2].yaxis.set_minor_locator(MultipleLocator(0.005)) 


    axs[3].set_ylim(-0.18,0.04)
    axs[3].set_yticks([-0.15,-0.10,-0.05,0])
    axs[3].set_yticklabels([r'$-0.15$',r'$-0.10$',r'$-0.05$',r'$0$'])
    axs[3].yaxis.set_minor_locator(MultipleLocator(0.01)) 


    axs[4].set_ylim(-0.035,0.03)
    axs[4].set_yticks([-0.03,-0.02,-0.01,0,0.01,0.02])
    axs[4].set_yticklabels([r'$-0.03$',r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$'])
    axs[4].yaxis.set_minor_locator(MultipleLocator(0.01)) 


    axs[5].set_ylim(-0.15,0.23)  
    axs[5].set_yticks([-0.1,0.0,0.1,0.2])
    axs[5].set_yticklabels([r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])
    axs[5].yaxis.set_minor_locator(MultipleLocator(0.05)) 


    for i in range(5):
        axs [i+1].axhline(0  ,color='k',alpha=0.5,zorder=6)
        axLs[i+1].axhline(0  ,color='k',alpha=0.5,zorder=6)
        # axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
        # axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

    axLs[4].set_xlabel(r'\boldmath$x$',size=45,loc = 'left')
    axLs[5].set_xlabel(r'\boldmath$x$',size=45,loc = 'left')   
    # axLs[5].xaxis.set_label_coords(0.95,0.00)
    # axLs[6].xaxis.set_label_coords(0.95,0.00)
      

    axs [1].text(0.13,0.85,r'\boldmath{$x \Delta u^+$}',                        transform=axs[1] .transAxes,size=40)
    axLs[2].text(0.35,0.85,r'\boldmath{$x \Delta \bar{u}$}',                    transform=axLs[2].transAxes,size=40)
    axs [3].text(0.13,0.85,r'\boldmath{$x \Delta d^+$}',                        transform=axs[3] .transAxes,size=40)
    axs [4].text(0.13,0.85,r'\boldmath{$x \Delta \bar{d}$}',                    transform=axs[4] .transAxes,size=40)
    axLs[5].text(0.35,0.85,r'\boldmath{$x \Delta g$}'  ,                        transform=axLs[5].transAxes,size=40)


    axLs[5].text(0.10,0.05,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axLs[5].transAxes,size=25)

    handles,labels = [thy],[r'\rm \bf JAMpol25']
    handles.append(hand['JAM17'])
    handles.append(hand['NNPDF'])
    handles.append(hand['DSSV14'])
    labels.append(SETS['JAM17'][3])
    labels.append(SETS['NNPDF'][3])
    labels.append(SETS['DSSV14'][3])
    axs[6].legend(handles,labels,loc=(0.00,0.05),fontsize=35,frameon=0,handletextpad=0.3,handlelength=1.0)
      
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.2,hspace=0.02)

    # filename = '%s/gallery/ppdfs-Q2=%3.5f'%(wdir,Q2)
    filename = 'plots/fig_ppdfs_compare'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def plot_ppdfs(wdir,Q2=10,SETS={}):

    nrows,ncols=3,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    axs,axLs = {},{}
    for i in range(6):
        axs[i+1] = py.subplot(nrows,ncols,i+1)
        divider = make_axes_locatable(axs[i+1])
        axLs[i+1] = divider.append_axes("right",size=2.75,pad=0,sharey=axs[i+1])
        axLs[i+1].set_xlim(0.1,0.9)
        axLs[i+1].spines['left'].set_visible(False)
        axLs[i+1].yaxis.set_ticks_position('right')
        py.setp(axLs[i+1].get_xticklabels(),visible=True)

        axs[i+1].spines['right'].set_visible(False)

    axs [6].axis("off")
    axLs[6].axis("off")

    thy  = {}
    hand = {}
    load_config('%s/input.py'%wdir)

    filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
    if regen: pdf.gen_xf(wdir,Q2)
    try:
        data=load(filename)
    except:
        pdf.gen_xf(wdir,Q2)
        data=load(filename)
    print('Loading %s'%filename)

    replicas=core.get_replicas(wdir)
    names    = core.get_replicas_names(wdir)

    X=data['X']

    for flav in data['XF']:
        mean = np.mean(data['XF'][flav],axis=0)
        std  = np.std (data['XF'][flav],axis=0)

        if flav=='up'      : ax,axL,zorder = axs[1],axLs[1], 1.0 
        elif flav=='ub'    : ax,axL,zorder = axs[2],axLs[2], 1.0 
        elif flav=='dp'    : ax,axL,zorder = axs[3],axLs[3], 1.0 
        elif flav=='db'    : ax,axL,zorder = axs[4],axLs[4], 1.0 
        elif flav=='g'     : ax,axL,zorder = axs[5],axLs[5], 1.0 
        
        else: continue


        #--plot each replica
        #--plot average and standard deviation
        thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=zorder)
        axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=1.0,zorder=zorder)
        thy_mean ,= ax.plot(X,mean,color='firebrick',alpha=1.0,zorder=1.5)
        axL           .plot(X,mean,color='firebrick',alpha=1.0,zorder=1.5)

        #--plot other PPDF sets
        for SET in SETS:
            _SET, _color, _alpha, _zorder = SETS[SET][0], SETS[SET][1], SETS[SET][2], SETS[SET][4]

            if SET=='DSSV08':
                _X=10**np.linspace(-4,-1,200)
                _X=np.append(_X,np.linspace(0.1,0.99,200))
                ub = _SET.xfxQ2(-2,_X,Q2)
                db = _SET.xfxQ2(-1,_X,Q2)
                if flav=='up'   : ppdf = _SET.xfxQ2(-2 ,_X,Q2) + _SET.xfxQ2( 2,_X,Q2) 
                elif flav=='dp' : ppdf = _SET.xfxQ2(-1 ,_X,Q2) + _SET.xfxQ2( 1,_X,Q2) 
                elif flav=='g'  : ppdf = _SET.xfxQ2( 21,_X,Q2)
                elif flav=='ub' : ppdf = _SET.xfxQ2(-2 ,_X,Q2)
                elif flav=='db' : ppdf = _SET.xfxQ2(-1 ,_X,Q2)
                mean = ppdf[0]
                std  = 0
                for i in range(1,20):
                    std += (ppdf[i] - ppdf[-i])**2
                std = np.sqrt(std)/2.0
                #hand[SET] = ax.fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=_zorder)
                #axL.           fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=_zorder)

            if SET=='NNPDF':
                ppdf = _SET.get_xpdf(flav,X,Q2)
                color = 'gold'
                hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=color,alpha=0.4,zorder=0.9)#,facecolor='none',hatch='//')
                axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=color,alpha=0.4,zorder=0.9)#,facecolor='none',hatch='//')
                hand[SET + ' mean'] ,= ax .plot(X,ppdf['xf0'],color=color,alpha=1.0,zorder=1.5)
                axL                       .plot(X,ppdf['xf0'],color=color,alpha=1.0,zorder=1.5)

            if SET=='JAM17':
                color = 'green'
                ppdf = _SET.get_xpdf(flav,X,Q2)
                hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=color,alpha=0.5,zorder=1.2,facecolor='none',hatch='|')
                axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=color,alpha=0.5,zorder=1.2,facecolor='none',hatch='|')
                hand[SET + ' mean'] ,= ax .plot(X,ppdf['xf0'],color=color,alpha=1.0,zorder=1.5)
                axL                       .plot(X,ppdf['xf0'],color=color,alpha=1.0,zorder=1.5)

            if SET=='DSSV14':
                ppdf = _SET.get_xpdf(flav,X,Q2)
                color = 'magenta'
                hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=color,alpha=0.3,zorder=1.1)
                axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=color,alpha=0.3,zorder=1.1)
                hand[SET+ ' mean'] ,= ax .plot(X,ppdf['xf0'],color='purple',alpha=1.0,zorder=1.5)
                axL.                      plot(X,ppdf['xf0'],color='purple',alpha=1.0,zorder=1.5)


    for i in range(5):
          axs[i+1].set_xlim(8e-3,0.1)
          axs[i+1].semilogx()

          axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
          axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
          axs[i+1].set_xticks([0.01,0.1])
          axs[i+1].set_xticklabels([r'$0.01$',r'$0.1$'])

          if i in [0,2,4]: axLs[i+1].set_xlim(0.1,0.71)
          if i in [1,3,5]: axLs[i+1].set_xlim(0.1,0.51)

          axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=8)
          axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=4)
          if i in [0,2,4]: axLs[i+1].set_xticks([0.1,0.3,0.5,0.7])
          if i in [0,2,4]: axLs[i+1].set_xticklabels([r'',r'$0.3$',r'$0.5$',r'$0.7$'])
          if i in [1,3,5]: axLs[i+1].set_xticks([0.1,0.3,0.5])
          if i in [1,3,5]: axLs[i+1].set_xticklabels([r'',r'$0.3$',r'$0.5$'])

          axLs[i+1].xaxis.set_minor_locator(MultipleLocator(0.10)) 

    for i in [1,2,3]:
          axs[i] .tick_params(labelbottom=False)
          axLs[i].tick_params(labelbottom=False)

    

    axs[1].set_ylim(-0.01,0.39)
    axs[1].set_yticks([0,0.1,0.2,0.3])
    axs[1].set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$'])
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.05)) 

    axs[2].set_ylim(-0.025,0.035)
    axs[2].set_yticks([-0.02,-0.01,0,0.01,0.02,0.03])
    axs[2].set_yticklabels([r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$',r'$0.03$'])
    axs[2].yaxis.set_minor_locator(MultipleLocator(0.005)) 


    axs[3].set_ylim(-0.18,0.04)
    axs[3].set_yticks([-0.15,-0.10,-0.05,0])
    axs[3].set_yticklabels([r'$-0.15$',r'$-0.10$',r'$-0.05$',r'$0$'])
    axs[3].yaxis.set_minor_locator(MultipleLocator(0.01)) 


    axs[4].set_ylim(-0.035,0.03)
    axs[4].set_yticks([-0.03,-0.02,-0.01,0,0.01,0.02])
    axs[4].set_yticklabels([r'$-0.03$',r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$'])
    axs[4].yaxis.set_minor_locator(MultipleLocator(0.01)) 


    axs[5].set_ylim(-0.15,0.23)  
    axs[5].set_yticks([-0.1,0.0,0.1,0.2])
    axs[5].set_yticklabels([r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])
    axs[5].yaxis.set_minor_locator(MultipleLocator(0.05)) 


    for i in range(5):
        axs [i+1].axhline(0  ,color='k',alpha=0.5,zorder=6)
        axLs[i+1].axhline(0  ,color='k',alpha=0.5,zorder=6)
        # axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
        # axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

    axLs[4].set_xlabel(r'\boldmath$x$',size=45,loc = 'left')
    axLs[5].set_xlabel(r'\boldmath$x$',size=45,loc = 'left')   
    # axLs[5].xaxis.set_label_coords(0.95,0.00)
    # axLs[6].xaxis.set_label_coords(0.95,0.00)
      

    axs [1].text(0.13,0.85,r'\boldmath{$x \Delta u^+$}',                        transform=axs[1] .transAxes,size=40)
    axLs[2].text(0.35,0.85,r'\boldmath{$x \Delta \bar{u}$}',                    transform=axLs[2].transAxes,size=40)
    axs [3].text(0.13,0.85,r'\boldmath{$x \Delta d^+$}',                        transform=axs[3] .transAxes,size=40)
    axs [4].text(0.13,0.85,r'\boldmath{$x \Delta \bar{d}$}',                    transform=axs[4] .transAxes,size=40)
    axLs[5].text(0.35,0.85,r'\boldmath{$x \Delta g$}'  ,                        transform=axLs[5].transAxes,size=40)


    axLs[5].text(0.10,0.05,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axLs[5].transAxes,size=25)

    handles,labels = [(thy,thy_mean)],[r'\rm \bf JAMpol25']
    handles.append((hand['JAM17'] ,hand['JAM17 mean']))
    handles.append((hand['NNPDF'] ,hand['NNPDF mean']))
    handles.append((hand['DSSV14'],hand['DSSV14 mean']))
    labels.append(SETS['JAM17'][3])
    labels.append(SETS['NNPDF'][3])
    labels.append(SETS['DSSV14'][3])
    axs[6].legend(handles,labels,loc=(0.00,0.05),fontsize=35,frameon=0,handletextpad=0.3,handlelength=1.0)
      
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.2,hspace=0.02)

    # filename = '%s/gallery/ppdfs-Q2=%3.5f'%(wdir,Q2)
    filename = 'plots/fig_ppdfs_compare'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)

if __name__=="__main__":


    msr_dir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g'

    # get PPDF sets for comparison
    SETS = {}
    #SETS['NNPDF'] = (QPDCALC('NNPDFpol11_100',ismc=True),'gold'    , 0.8 ,r'\textrm{NNPDFpol1.1}', 1)
    SETS['NNPDF'] = (QPDCALC('NNPDFpol20_nnlo_as_01180',ismc=True),'gold'    , 0.8 ,r'\textrm{NNPDFpol2.0}', 1)
    SETS['JAM17'] = (QPDCALC('JAM17_PPDF_nlo',ismc=True),'limegreen', 0.8 ,r'\textrm{JAM17}', 1.1)
    SETS['DSSV14']  = (QPDCALC('DSSV_REP_LHAPDF6',ismc=True),'dodgerblue', 0.8 ,r'\textrm{DSSV14}', 1.2)
    # SETS['DSSV']  = (DSSV(),                             'darkblue', 0.6 ,r'\textrm{\textbf{DSSV08}}')
    plot_ppdfs(msr_dir, Q2 = 10,SETS = SETS)


