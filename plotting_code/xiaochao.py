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
from analysis.obslib import stf,pstf


#--from obslib
from obslib.psidis.reader import READER

import kmeanconf as kc

regen = False

ext = '.png'
#ext = '.pdf'

cwd = os.getcwd()

def plot_JLab_E12_06_110(WDIR):

    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    hand = {}
    thy_plot = {}
    thy_band = {}

    j = 0
    for wdir in WDIR:
        load_config('%s/input.py' % wdir)
        istep = core.get_istep()
        replicas = core.get_replicas(wdir)
        core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

        predictions = load('%s/data/predictions-%d.dat' % (wdir, istep))

        data = predictions['reactions']['pidis']

        for idx in data:
            predictions = copy.copy(data[idx]['prediction-rep'])
            del data[idx]['prediction-rep']
            del data[idx]['residuals-rep']
            del data[idx]['shift-rep']
            del data[idx]['rres-rep']
            del data[idx]['r-residuals']
            del data[idx]['n-residuals']
            data[idx]['thy']  = np.mean(predictions, axis = 0)
            data[idx]['dthy'] = np.std (predictions, axis = 0)
            if 'X' in data[idx]: data[idx]['x'] = data[idx]['X']


        DATA = {}
        DATA['JLab E12 Apa']  = pd.DataFrame(data[20005]) #--
        DATA['JLab E12 Ape']  = pd.DataFrame(data[20006]) #--

        combo1 = ('firebrick', '*', 10)
        combo2 = ('darkgreen', '^', 8)
        #--plot data points
        for exp in DATA:
            if exp=='JLab E12 Apa': ax = ax11 
            if exp=='JLab E12 Ape': ax = ax12 
            ax.axhline(0,0,1,color='black',alpha=0.1)
            X     = DATA[exp]['X']
            val   = DATA[exp]['value']
            alpha = DATA[exp]['alpha']
            color,marker,ms = combo1 
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none',zorder=1.5)

            #--plot theory interpolated between data points
            X_thy= DATA[exp]['X']
            mean = DATA[exp]['thy']
            std  = DATA[exp]['dthy']
            down = mean - std
            up   = mean + std
            if j==0: alpha,zorder,color=0.6,1.1,'firebrick'
            if j==1: alpha,zorder,color=0.3,1.0,'cyan'
            thy_plot[j] ,= ax.plot(X_thy,mean,linestyle='solid',color=color,zorder=1.5)
            thy_band[j]  = ax.fill_between(X_thy,down,up,color=color,alpha=alpha,zorder=zorder)
            #if len(X)==1:
            #    x = [X_thy[0]*0.95,X_thy[0]*1.05]
            #    down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
            #    thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
            #    thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)
        j += 1

    for ax in [ax11,ax12]:
        ax.tick_params(axis = 'both', labelsize = 40)

        ax.yaxis.set_tick_params(which = 'major', length = 10)
        ax.yaxis.set_tick_params(which = 'minor', length = 5)

        ax.xaxis.set_tick_params(which = 'major', length = 10)
        ax.xaxis.set_tick_params(which = 'minor', length = 5)

        ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=60)
        ax.xaxis.set_label_coords(0.95,0.00)
        ax.tick_params(axis='x',pad=8)

        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
        ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


        ax.set_xlim(0.41, 0.77)
        ax.set_xticks([0.5,0.6,0.7])
        #ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.2$',r'$0.5$'])
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax11.set_ylim(-0.023, 0.049)
    ax11.set_yticks([-0.02,-0.01,0,0.01,0.02,0.03,0.04])
    ax11.set_yticklabels([r'$-0.02$',r'$-0.01$',r'$0$',r'$0.01$',r'$0.02$',r'$0.03$',r'$0.04$'])
    ax11.yaxis.set_minor_locator(MultipleLocator(0.005))

    ax12.set_ylim(-0.059, 0.059)
    ax12.set_yticks([-0.04,-0.02,0,0.02,0.04])
    ax12.set_yticklabels([r'$-0.04$',r'$-0.02$',r'$0$',r'$0.02$',r'$0.04$'])
    ax12.yaxis.set_minor_locator(MultipleLocator(0.005))


    #ax11.text(0.03,0.93    ,r'\textrm{\textbf{JLab E12-06-110}}' , transform = ax11.transAxes, size = 50)
    ax12.text(0.03,0.80    ,r'\textrm{\textbf{PRELIMINARY}}'     , transform = ax12.transAxes, size = 50)

    ax11.text(0.03,0.90    ,r'\boldmath$A^{^3 {\rm He}}_{\parallel}$'     , transform = ax11.transAxes, size = 90)
    ax12.text(0.03,0.90    ,r'\boldmath$A^{^3 {\rm He}}_{\perp}$'         , transform = ax12.transAxes, size = 90)


    handles,labels = [],[]
    handles.append(hand['JLab E12 Apa'])
    handles.append((thy_band[1],thy_plot[1]))
    handles.append((thy_band[0],thy_plot[0]))
    labels.append(r'\textrm{\textbf{JLab E12-06-110}}')
    labels.append(r'\textrm{\textbf{JAMpol25}')
    labels.append(r'\textrm{\textbf{+JLab E12-06-110}}')
    ax11.legend(handles,labels,loc=(0.00,0.60), fontsize = 45, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #ax11.title(r'\boldmath{$^3 {\rm He}$} \rm \bf target',fontsize = 60)


    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/plots/xiaochao_JLab_E12_06_110'%cwd
    filename += ext
    ax11.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

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
        ax.set_xlim(0.00,0.76)
        ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6])
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
    labels.append(r'\textrm{\textbf{+E12-06-110}}')
    ax12.legend(handles,labels,loc=(0.30,0.00),fontsize=35,frameon=0,handletextpad=0.3,handlelength=1.0)
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.2,hspace=0.02,top=0.99,right=0.99)

    filename = 'plots/xiaochao_asym_impact'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def ppdf_impact(WDIR,Q2=10):

    flavs = ['u','d','ub','db','g','sigma','sp']
    
    nrows = 4
    ncols = 2
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*10,nrows*5))
 
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            if j ==0: axs[i][j].set_xlim(0.01,0.76)
            if j ==1: axs[i][j].set_xlim(0.01,0.76)
            axs [i][j].axhline(0,lw = 1,color = 'k',alpha = 0.3,zorder=10)

    axs[3][1].axis("off")

    k = 0
    hand = {}
    for wdir in WDIR:
        filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
        data=load(filename)
        print('Loading %s'%filename)

        X = data['X']
        for flav in flavs:

            if   flav=='u':     i, j = 0, 0
            elif flav=='ub':    i, j = 0, 1
            elif flav=='d':     i, j = 1, 0
            elif flav=='db':    i, j = 1, 1
            elif flav=='sp':    i, j = 2, 0
            elif flav=='sigma': i, j = 2, 1
            elif flav=='g':     i, j = 3, 0

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

    axs[2][0].set_ylim(-0.12,0.12)
    axs[2][0].yaxis.set_minor_locator(MultipleLocator(0.02)) 
    axs[2][0].set_yticks([-0.08,-0.04,0,0.04,0.08]) 
    axs[2][0].set_yticklabels([r'$-0.08$',r'$-0.04$',r'$0$',r'$0.04$',r'$0.08$']) 

    axs[2][1].set_ylim(-0.05,0.35)
    axs[2][1].yaxis.set_minor_locator(MultipleLocator(0.05)) 
    axs[2][1].set_yticks([0,0.1,0.2,0.3]) 
    axs[2][1].set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$']) 

    axs[3][0].set_ylim(0.00,0.16)
    axs[3][0].yaxis.set_minor_locator(MultipleLocator(0.01)) 
    axs[3][0].set_yticks([0,0.05,0.10,0.15]) 
    axs[3][0].set_yticklabels([r'$0$',r'$0.05$',r'$0.10$',r'$0.15$']) 
 
    ls = 40
    for i in range(nrows):
        for j in range(ncols):
            if i==3 and j==1: continue
            axs [i][j].tick_params(axis='both', which='major', top=True, direction='in',labelsize=ls,length=10)
            axs [i][j].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=ls,length=5)
            axs[i][j].set_xticks([0.2,0.4,0.6])
            if i in [1,2,3] and j==0: axs[i][j].set_xticklabels([])
            elif i in [2] and j==1:
                axs[i][j].set_xticklabels([r'$0.2$',r'$0.4$',r'$0.6$'])
            else: 
                axs[i][j].set_xticklabels([r'$0.2$',r'$0.4$',r'$0.6$'])
    
            if i == 2 and j==1: axs[i][j].set_xlabel(r'\boldmath{$x$}',fontsize = 60)
            if i == 3 and j==0: axs[i][j].set_xlabel(r'\boldmath{$x$}',fontsize = 60)
    
            axs[i][j].xaxis.set_minor_locator(MultipleLocator(0.10)) 
       
    pdf_fs = 50
 
    axs[0][0].text( 0.77, 0.85, r'\boldmath$x \Delta u$',          transform=axs[0][0].transAxes,fontsize=pdf_fs)
    axs[1][0].text( 0.77, 0.05, r'\boldmath$x \Delta d$',          transform=axs[1][0].transAxes,fontsize=pdf_fs)
    axs[0][1].text( 0.77, 0.85, r'\boldmath$x \Delta \bar{u}$',    transform=axs[0][1].transAxes,fontsize=pdf_fs)
    axs[1][1].text( 0.77, 0.05, r'\boldmath$x \Delta \bar{d}$',    transform=axs[1][1].transAxes,fontsize=pdf_fs)
    axs[2][0].text( 0.77, 0.05, r'\boldmath$x \Delta s^+$',        transform=axs[2][0].transAxes,fontsize=pdf_fs)
    axs[2][1].text( 0.77, 0.85, r'\boldmath$x \Delta \Sigma$',     transform=axs[2][1].transAxes,fontsize=pdf_fs)
    axs[3][0].text( 0.77, 0.85, r'\boldmath$x \Delta g$',          transform=axs[3][0].transAxes,fontsize=pdf_fs)

    axs[3][0].text(0.03, 0.05, r'$Q^2 = 10 ~ {\rm GeV}^2$',     transform=axs[3][0].transAxes,fontsize=35)
  

    handles, labels = [],[]
    handles.append(hand[1])
    handles.append(hand[0])
    labels.append(r'\textbf{\textbf{JAM}}')
    labels.append(r'\textbf{\textbf{+E12-06-110}}')
    axs[3][1].legend(handles,labels,fontsize = 40, loc=(0.10,0.00), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.02,wspace=0.20,top=0.99,right=0.99)

    filename = 'plots/xiaochao_ppdfs' + ext 

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

def polarization_impact(WDIR,Q2=10,order='NLO'):

    #M2 = 0.938**2 
    #regen = False 
    #W2min  = np.array([4,4])
    #xmax   = Q2/(Q2+W2min-M2)
 
    ncols = 2
    nrows = 1
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*7,nrows*5))


    k = 0
    flavs = ['up','dp']
    hand = {}

    colors = ['red','cyan']
    alpha = [0.9,0.9]
    zorders = [1.1,1.0]

    DATA = {}

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

            if   flav=='up':     i = 0
            elif flav=='dp':     i = 1

            if flav=='up':
                result = (np.array(data['XF']['u']) + np.array(data['XF']['ub']))/(np.array(udata['XF']['u']) + np.array(udata['XF']['ub']))
            elif flav=='dp':
                result = (np.array(data['XF']['d']) + np.array(data['XF']['db']))/(np.array(udata['XF']['d']) + np.array(udata['XF']['db']))

            #result[:,X>xmax[k]] = 'nan'

            mean = np.mean(result,axis=0)
            std  = np.std (result,axis=0)

            if k==2:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color='none',edgecolor=colors[k],alpha=alpha[k],zorder=zorders[k],hatch='//')
            else:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color=colors[k],alpha=alpha[k],zorder=zorders[k])

            if k==0:
                DATA['d%s/%s after mean'%(flav,flav)] = mean
                DATA['d%s/%s after std' %(flav,flav)] = std
            if k==1:
                DATA['d%s/%s before mean'%(flav,flav)] = mean
                DATA['d%s/%s before std' %(flav,flav)] = std
        k+=1 
   
    #--plot LO calculation
    wdir = WDIR[0]

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


    if order=='NLO': filename = wdir + '/data/stf-Q2=%3.5f.dat'%(Q2)
    if order=='LO':  filename = wdir + '/data/stf-Q2=%3.5f-LO.dat'%(Q2)
    if regen: stf.gen_stf(wdir,Q2=Q2,order=order) 
    try: STF=load(filename)
    except:
        stf.gen_stf(wdir,Q2=Q2,order=order)
        STF=load(filename)

    PSTF = {}
    for tar in ['p','n']:
        if order=='NLO': filename = wdir + '/data/pstf-%s-%s-Q2=%3.5f.dat'   %('g1',tar,Q2)
        if order=='LO':  filename = wdir + '/data/pstf-%s-%s-Q2=%3.5f-LO.dat'%('g1',tar,Q2)
        if regen:
            pstf.gen_pstf(wdir,Q2=Q2,tar=tar,stf='g1',order=order) 
        try: 
            PSTF[tar]=load(filename)
        except:
            pstf.gen_pstf(wdir,Q2=Q2,tar=tar,stf='g1',order=order)
            PSTF[tar]=load(filename)

    X  = STF['X']
    Q2 = STF['Q2']
    F2p = STF['XF']['p']['F2']/X
    F2n = STF['XF']['n']['F2']/X
    FLp = STF['XF']['p']['FL']/X
    FLn = STF['XF']['n']['FL']/X

    M2 = 0.938**2
    rho2 = 1 + 4*M2*X**2/Q2
    F1p = (rho2 * F2p - FLp)/(2*X) 
    F1n = (rho2 * F2n - FLn)/(2*X)

    g1p = PSTF['p']['STF']
    g1n = PSTF['n']['STF']

    u  = np.array(udata['XF']['u'])
    ub = np.array(udata['XF']['ub'])
    d  = np.array(udata['XF']['d'])
    db = np.array(udata['XF']['db'])
    rat = (d + db)/(u + ub) 

    test = np.std(F1p,axis=0)

    u_LO_eq = 4/15 * g1p/F1p * (4 + 1*rat) - 1/15 * g1n/F1n * (1 + 4*rat)
    d_LO_eq = 4/15 * g1n/F1n * (4 + 1/rat) - 1/15 * g1p/F1p * (1 + 4/rat)

    u_mean = np.mean(u_LO_eq,axis=0)
    u_std  = np.std (u_LO_eq,axis=0)
    d_mean = np.mean(d_LO_eq,axis=0)
    d_std  = np.std (d_LO_eq,axis=0)


    hand['LO eq'] = axs [0].fill_between(X,u_mean-u_std,u_mean+u_std,color='green',alpha=0.5,zorder=1.5)
    hand['LO eq'] = axs [1].fill_between(X,d_mean-d_std,d_mean+d_std,color='green',alpha=0.5,zorder=1.5)
 
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
  

    axs[0].text(0.10, 0.05, r'\textrm{\textbf{PRELIMINARY}}',     transform=axs[0].transAxes,fontsize=30)
    axs[0].text(0.03, 0.80, r'\boldmath$\frac{\Delta u^+}{u^+}$', transform=axs[0].transAxes,fontsize=50)
    axs[1].text(0.03, 0.10, r'\boldmath$\frac{\Delta d^+}{d^+}$', transform=axs[1].transAxes,fontsize=50)
    
    axs[0].text(0.70, 0.05, r'$Q^2 = %d ~ {\rm GeV}^2$'%Q2[0],     transform=axs[0].transAxes,fontsize=20)


    handles,labels = [],[]
    handles.append(hand[1])
    handles.append(hand[0])
    handles.append(hand['LO eq'])
    labels.append(r'\textrm{\textbf{JAM}}')
    labels.append(r'\textrm{\textbf{+E12-06-110}}')
    labels.append(r'\textrm{\textbf{LO eq. (%s stfs)}}'%order)
    axs[1].legend(handles,labels,fontsize = 25, loc=(0.00,0.55), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.01,wspace=0.15,top=0.97,right=0.99)

    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
 
    filename = 'plots/xiaochao_polarization_impact-%s'%order + ext 
    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

    #--save data
    DATA['X']      = X
    DATA['Q2']     = Q2
    DATA['LO eq u mean']  = u_mean
    DATA['LO eq u std']   = u_std
    DATA['LO eq d mean']  = d_mean
    DATA['LO eq d std']   = d_std

    DATA = pd.DataFrame(DATA)
  
    checkdir('%s/data'%cwd)
    filename = '%s/data/xiaochao_polarization-%s.csv'%(cwd,order)
    DATA.to_csv(filename,index=False)
    print('Saving polarization data to %s'%filename)

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
    axsL[0][0].text(-0.90, 0.60, r'\boldmath$\frac{\delta~({\rm +E12-06-110})}{\delta~({\rm JAM})}$',          transform=axsL[0][0].transAxes,fontsize=45)
  

    #handles, labels = [],[]
    #handles.append(hand[2])
    #handles.append(hand[1])
    #labels.append(r'\boldmath{$W^2 > 4$}  \rm \bf GeV\boldmath{$^2$}')
    #labels.append(r'\boldmath{$W^2 > 10$} \rm \bf GeV\boldmath{$^2$}')
    #axsL[0][0].legend(handles,labels,fontsize = 26, loc=(-1.10,0.72), frameon=0,handlelength=1,handletextpad=0.6)
  

    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.02,wspace=0.02,top=0.99,right=0.99,left=0.05)

    filename = 'plots/xiaochao_ppdfs_rel_uncert' + ext 

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


def plot_A1h(WDIR,Q2=5):

    nrows,ncols=1,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*7,nrows*5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    hand = {}

    wdir = WDIR[0]
    load_config(wdir + '/input.py')

    filename = wdir + '/data/stf-Q2=%3.5f.dat'%(Q2)
    if regen: stf.gen_stf(wdir,Q2=Q2) 
    try: STF=load(filename)
    except:
        stf.gen_stf(wdir,Q2=Q2)
        STF=load(filename)

    PSTF = {}
    for tar in ['p','n','h']:
        PSTF[tar] = {}
        for _stf in ['g1','g2']:
            filename = wdir + '/data/pstf-%s-%s-Q2=%3.5f.dat'   %(_stf,tar,Q2)
            if regen:
                pstf.gen_pstf(wdir,Q2=Q2,tar=tar,stf=_stf) 
            try: 
                PSTF[tar][_stf]=load(filename)
            except:
                pstf.gen_pstf(wdir,Q2=Q2,tar=tar,stf=_stf)
                PSTF[tar][_stf]=load(filename)

    X  = STF['X']
    Q2 = STF['Q2']
    F2p = STF['XF']['p']['F2']/X
    F2n = STF['XF']['n']['F2']/X
    F2h = STF['XF']['h']['F2']/X
    FLp = STF['XF']['p']['FL']/X
    FLn = STF['XF']['n']['FL']/X
    FLh = STF['XF']['h']['FL']/X

    M2 = 0.938**2
    rho2 = 1 + 4*M2*X**2/Q2
    F1p = (rho2 * F2p - FLp)/(2*X) 
    F1n = (rho2 * F2n - FLn)/(2*X)
    F1h = (rho2 * F2h - FLh)/(2*X)

    g1p = PSTF['p']['g1']['STF']
    g1n = PSTF['n']['g1']['STF']
    g1h = PSTF['h']['g1']['STF']
    g2p = PSTF['p']['g2']['STF']
    g2n = PSTF['n']['g2']['STF']
    g2h = PSTF['h']['g2']['STF']

    #--free neutron result
    A1n =  1/F1n * g1n - 1/F1n * (rho2 - 1) * g2n
    A2n =  np.sqrt(rho2-1)/F1n * (g1n + g2n)
  
    A1n_mean    = np.mean(A1n   ,axis=0) 
    A1n_std     = np.std (A1n   ,axis=0)
    A2n_mean    = np.mean(A2n   ,axis=0) 
    A2n_std     = np.std (A2n   ,axis=0)

    hand['n'] = ax11.fill_between(X,(A1n_mean-A1n_std),(A1n_mean+A1n_std   ),color='darkgreen',alpha=0.5,zorder=1.0)
    hand['n'] = ax12.fill_between(X,(A2n_mean-A2n_std),(A2n_mean+A2n_std   ),color='darkgreen',alpha=0.5,zorder=1.0)


    #--full smearing result


    A1h =  1/F1h * g1h - 1/F1h * (rho2 - 1) * g2h
    A2h =  np.sqrt(rho2-1)/F1h * (g1h + g2h)
 
    #--scaling factor
    factor = 1 + 2*F2p/F1n
    A1h *= 3*factor 
    A2h *= 3*factor 

    A1h_mean    = np.mean(A1h   ,axis=0) 
    A1h_std     = np.std (A1h   ,axis=0)
    A2h_mean    = np.mean(A2h   ,axis=0) 
    A2h_std     = np.std (A2h   ,axis=0)

    hand['full'] = ax11.fill_between(X,(A1h_mean-A1h_std),(A1h_mean+A1h_std),color='firebrick',alpha=0.5,zorder=1.0)
    hand['full'] = ax12.fill_between(X,(A2h_mean-A2h_std),(A2h_mean+A2h_std),color='firebrick',alpha=0.5,zorder=1.0)

    #--effective smearing result
    if 'pidis' not in conf:
        resman = RESMAN(nworkers=1,parallel=False,datasets=False)
        resman.setup_pidis()
    pol = conf['pidis'].pol['h']

    p, n = 2, 1
    
    F1h = (p * F1p + n * F1n)/(p+n)
    g1h = (pol['p'] * g1p + pol['n'] * g1n)/(p+n)

    A1h =  1/F1h * g1h - 1/F1h * (rho2 - 1) * g2h
    A2h =  np.sqrt(rho2-1)/F1h * (g1h + g2h)
    
    factor = 1 + 2*F2p/F1n
    A1h *= 3*factor 
    A2h *= 3*factor 
   
    A1h_mean    = np.mean(A1h   ,axis=0) 
    A1h_std     = np.std (A1h   ,axis=0)
    A2h_mean    = np.mean(A2h   ,axis=0) 
    A2h_std     = np.std (A2h   ,axis=0)

    hand['eff'] = ax11.fill_between(X,(A1h_mean-A1h_std),(A1h_mean+A1h_std),color='blue',alpha=0.5,zorder=1.0)
    hand['eff'] = ax12.fill_between(X,(A2h_mean-A2h_std),(A2h_mean+A2h_std),color='blue',alpha=0.5,zorder=1.0)


    for ax in [ax11,ax12]:
        ax.set_xlim(0.01,0.90)
        ax.set_xticks([0.2,0.4,0.6,0.8])
        ax.xaxis.set_minor_locator(MultipleLocator(0.10))


        ax.tick_params(axis='both', which='major', top=True, right=True, direction='in',labelsize=25,length=8)
        ax.tick_params(axis='both', which='minor', top=True, right=True, direction='in',labelsize=25,length=4)
        ax.axhline(0  ,color='k',alpha=0.5,zorder=6)

        ax.axhline(1  ,0.8,1.0,color='k',alpha=0.5,ls=':',zorder=6)

    ax11.set_ylim(-0.19,0.19)
    ax11.set_yticks([-0.15,-0.10,-0.05,0,0.05,0.10,0.15])
    ax11.set_yticklabels([r'$0.15$',r'$-0.10$',r'$-0.05$',r'$0$',r'$0.05$',r'$0.10$',r'$0.15$'])
    ax11.yaxis.set_minor_locator(MultipleLocator(0.01)) 

    ax12.set_ylim(-0.009,0.009)
    ax12.set_yticks([-0.008,-0.004,0,0.004,0.008])
    ax12.set_yticklabels([r'$-0.008$',r'$-0.004$',r'$0$',r'$0.004$',r'$0.008$'])
    ax12.yaxis.set_minor_locator(MultipleLocator(0.001)) 

    #ax12.tick_params(labelleft=False)

    ax11.set_xlabel(r'\boldmath$x$',size=45,loc = 'left')
    ax11.xaxis.set_label_coords(0.95,0.00)

    ax11.text(0.05,0.85,r'\boldmath$A_1$',                    transform=ax11.transAxes,size=50)
    ax12.text(0.05,0.85,r'\boldmath$A_2$',                    transform=ax12.transAxes,size=50)
    #ax12.text(0.35,0.85,r'\boldmath{$x (\Delta \bar{u} + \Delta \bar{d})$}',                    transform=ax12.transAxes,size=45)

    ax12.text(0.02,0.05,r'$Q^2 = %d$~'%Q2[0] + r'\textrm{GeV}'+r'$^2$', transform=ax12.transAxes,size=25)

    handles,labels = [],[]
    handles.append(hand['n'])
    handles.append(hand['eff'])
    handles.append(hand['full'])
    labels.append(r'\textrm{\textbf{free n}}')
    labels.append(r'\textrm{\textbf{eff pol}}')
    labels.append(r'\textrm{\textbf{smearing}}')
    ax11.legend(handles,labels,loc=(0.00,0.50),fontsize=22,frameon=0,handletextpad=0.3,handlelength=1.0)
 
    py.tight_layout()
    py.subplots_adjust(wspace=0.20,hspace=0.02,top=0.99,right=0.99)

    filename = 'plots/xiaochao_A1h'

    filename+=ext

    checkdir('%s/gallery'%wdir)
    py.savefig(filename)
    print ('Saving figure to %s'%filename)


def polarization_impact_pos(WDIR,Q2=10,order='NLO'):

    #M2 = 0.938**2 
    #regen = False 
    #W2min  = np.array([4,4])
    #xmax   = Q2/(Q2+W2min-M2)
 
    ncols = 2
    nrows = 1
    fig,axs = plt.subplots(nrows,ncols,figsize=(ncols*7,nrows*5))


    k = 0
    flavs = ['up','dp']
    hand = {}

    colors = ['red','cyan']
    alpha = [0.9,0.9]
    zorders = [1.1,1.0]

    DATA = {}

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

            if   flav=='up':     i = 0
            elif flav=='dp':     i = 1

            if flav=='up':
                result = (np.array(data['XF']['u']) + np.array(data['XF']['ub']))/(np.array(udata['XF']['u']) + np.array(udata['XF']['ub']))
            elif flav=='dp':
                result = (np.array(data['XF']['d']) + np.array(data['XF']['db']))/(np.array(udata['XF']['d']) + np.array(udata['XF']['db']))

            #result[:,X>xmax[k]] = 'nan'

            mean = np.mean(result,axis=0)
            std  = np.std (result,axis=0)

            if k==2:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color='none',edgecolor=colors[k],alpha=alpha[k],zorder=zorders[k],hatch='//')
            else:
                hand[k] = axs [i].fill_between(X,mean-std,mean+std,color=colors[k],alpha=alpha[k],zorder=zorders[k])

            if k==0:
                DATA['d%s/%s after mean'%(flav,flav)] = mean
                DATA['d%s/%s after std' %(flav,flav)] = std
            if k==1:
                DATA['d%s/%s before mean'%(flav,flav)] = mean
                DATA['d%s/%s before std' %(flav,flav)] = std
        k+=1 
   
    #--plot LO calculation
    wdir = WDIR[0]

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


    if order=='NLO': filename = wdir + '/data/stf-Q2=%3.5f.dat'%(Q2)
    if order=='LO':  filename = wdir + '/data/stf-Q2=%3.5f-LO.dat'%(Q2)
    if regen: stf.gen_stf(wdir,Q2=Q2,order=order) 
    try: STF=load(filename)
    except:
        stf.gen_stf(wdir,Q2=Q2,order=order)
        STF=load(filename)

    PSTF = {}
    for tar in ['p','n']:
        if order=='NLO': filename = wdir + '/data/pstf-%s-%s-Q2=%3.5f.dat'   %('g1',tar,Q2)
        if order=='LO':  filename = wdir + '/data/pstf-%s-%s-Q2=%3.5f-LO.dat'%('g1',tar,Q2)
        if regen:
            pstf.gen_pstf(wdir,Q2=Q2,tar=tar,stf='g1',order=order) 
        try: 
            PSTF[tar]=load(filename)
        except:
            pstf.gen_pstf(wdir,Q2=Q2,tar=tar,stf='g1',order=order)
            PSTF[tar]=load(filename)

    X  = STF['X']
    Q2 = STF['Q2']
    F2p = STF['XF']['p']['F2']/X
    F2n = STF['XF']['n']['F2']/X
    FLp = STF['XF']['p']['FL']/X
    FLn = STF['XF']['n']['FL']/X

    M2 = 0.938**2
    rho2 = 1 + 4*M2*X**2/Q2
    F1p = (rho2 * F2p - FLp)/(2*X) 
    F1n = (rho2 * F2n - FLn)/(2*X)

    g1p = PSTF['p']['STF']
    g1n = PSTF['n']['STF']

    u  = np.array(udata['XF']['u'])
    ub = np.array(udata['XF']['ub'])
    d  = np.array(udata['XF']['d'])
    db = np.array(udata['XF']['db'])
    rat = (d + db)/(u + ub) 

    test = np.std(F1p,axis=0)

    u_LO_eq = 4/15 * g1p/F1p * (4 + 1*rat) - 1/15 * g1n/F1n * (1 + 4*rat)
    d_LO_eq = 4/15 * g1n/F1n * (4 + 1/rat) - 1/15 * g1p/F1p * (1 + 4/rat)

    u_mean = np.mean(u_LO_eq,axis=0)
    u_std  = np.std (u_LO_eq,axis=0)
    d_mean = np.mean(d_LO_eq,axis=0)
    d_std  = np.std (d_LO_eq,axis=0)


    hand['LO eq'] = axs [0].fill_between(X,u_mean-u_std,u_mean+u_std,color='green',alpha=0.5,zorder=1.5)
    hand['LO eq'] = axs [1].fill_between(X,d_mean-d_std,d_mean+d_std,color='green',alpha=0.5,zorder=1.5)
 
    ls = 25 

    for i in range(len(axs)):
        axs[i].set_xlim(0.01,0.80)
        axs [i].axhline(0,lw = 1,color = 'k',alpha = 0.3)
        axs [i].tick_params(axis='both', which='major', top=False, direction='in',labelsize=ls,length=10)
        axs [i].tick_params(axis='both', which='minor', top=False, direction='in',labelsize=ls,length=5)
        axs [i].set_xticks([0.2,0.4,0.6])
        axs[i].set_xlabel(r'\boldmath{$x$}',fontsize = 35)
        axs[i].xaxis.set_minor_locator(MultipleLocator(0.10)) 

    axs[0].set_ylim(0.001,1.49)
    axs[0].yaxis.set_minor_locator(MultipleLocator(0.10)) 
    axs[0].set_yticks([0.2,0.4,0.6,0.8,1.0]) 
    axs[0].set_yticklabels([r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$',r'$1.0$']) 
   
    #axs[1].set_ylim(-1.10,1.00)
    axs[1].set_ylim(-1.99,1.99)
    axs[1].yaxis.set_minor_locator(MultipleLocator(0.10)) 
    axs[1].set_yticks([-1.0,-0.5,0,0.5,1.0]) 
    axs[1].set_yticklabels([r'$-1.0$',r'$-0.5$',r'$0$',r'$0.5$',r'1.0']) 
  

    axs[0].text(0.10, 0.05, r'\textrm{\textbf{PRELIMINARY}}',     transform=axs[0].transAxes,fontsize=30)
    axs[0].text(0.03, 0.80, r'\boldmath$\frac{\Delta u^+}{u^+}$', transform=axs[0].transAxes,fontsize=50)
    axs[1].text(0.03, 0.10, r'\boldmath$\frac{\Delta d^+}{d^+}$', transform=axs[1].transAxes,fontsize=50)
    
    axs[0].text(0.70, 0.05, r'$Q^2 = %d ~ {\rm GeV}^2$'%Q2[0],     transform=axs[0].transAxes,fontsize=20)


    handles,labels = [],[]
    #handles.append(hand[1])
    handles.append(hand[0])
    handles.append(hand['LO eq'])
    labels.append(r'\textrm{\textbf{JAM pos}}')
    #labels.append(r'\textrm{\textbf{+E12-06-110}}')
    labels.append(r'\textrm{\textbf{LO eq. (%s stfs)}}'%order)
    axs[1].legend(handles,labels,fontsize = 25, loc=(0.00,0.55), frameon=0,handlelength=1,handletextpad=0.6)
 
    plt.tight_layout() 
    plt.subplots_adjust(hspace=0.01,wspace=0.15,top=0.97,right=0.99)

    axs [0].set_rasterized(True)
    axs [1].set_rasterized(True)
 
    filename = 'plots/xiaochao_polarization_impact-%s'%order + ext 
    plt.savefig(filename)
    print('Saving figure to %s'%filename)
    plt.close()

    #--save data
    DATA['X']      = X
    DATA['Q2']     = Q2
    DATA['LO eq u mean']  = u_mean
    DATA['LO eq u std']   = u_std
    DATA['LO eq d mean']  = d_mean
    DATA['LO eq d std']   = d_std

    DATA = pd.DataFrame(DATA)
  
    checkdir('%s/data'%cwd)
    filename = '%s/data/xiaochao_polarization-%s.csv'%(cwd,order)
    DATA.to_csv(filename,index=False)
    print('Saving polarization data to %s'%filename)


if __name__ == "__main__":

    wdir1 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/jlab12_pdfpos'
    wdir2 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/jlab12_predict'
    
    WDIR = [wdir1,wdir2]
    #plot_JLab_E12_06_110(WDIR)

    wdir2 =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/HT4_pdfpos'
    WDIR = [wdir1,wdir2]


    polarization_impact(WDIR,Q2=10,order='NLO')
    polarization_impact(WDIR,Q2=10,order='LO')
    plot_ppdfs_rel_uncert(WDIR, Q2=10)
    ppdf_impact(WDIR, Q2=10)
    asym_impact(WDIR, Q2=10)


    #WDIR = [wdir2]
    #polarization_impact_pos(WDIR,Q2=10,order='LO' )
    #polarization_impact_pos(WDIR,Q2=10,order='NLO')


