#!/usr/bin/env python
import os,sys
from tools.config import load_config,conf
from fitlib.resman import RESMAN
import numpy as np

#--matplotlib
import matplotlib
import pylab  as py
from matplotlib.ticker import MultipleLocator

#--from tools
from tools.tools import load,lprint,save

#--from corelib
from analysis.corelib import core,classifier

from qcdlib import pdf as PDF
from obslib.idis.theory import OFFSHELL_MODEL

#--generate nuclear PDFs
def gen_nuclear_pdf(wdir,Q2):

    replicas = core.get_replicas(wdir)

    load_config('%s/input.py'%wdir)

    istep = core.get_istep()
    core.mod_conf(istep,replicas[0])
    resman=RESMAN(parallel=False,datasets=False)
    parman = resman.parman

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']
    names    = core.get_replicas_names(wdir)
    resman.setup_idis()
    
    if 'off pdf' not in conf['steps'][istep]['active distributions']: 
        if 'off pdf' not in conf['steps'][istep]['passive distributions']: 
            return

    idis=resman.idis_thy
    off = conf['off pdf']
    pdf = conf['pdf']
    offshell_model = OFFSHELL_MODEL()

    ##############################################
    #--generate offshell
    ##############################################
    X=10**np.linspace(-6,-1,200)
    X=np.append(X,np.linspace(0.101,0.99,200))
    q2 = Q2
    Q2 = Q2*np.ones(len(X))

    #nuclei = ['p','n']
    nuclei = ['p']
    NPDF = {_: {} for _ in ['u','d']}
    NPDF['u'] = {_:{} for _ in nuclei}
    NPDF['d'] = {_:{} for _ in nuclei}
    for _ in nuclei:
        NPDF['u'][_]['d'] = {_: [] for _ in ['onshell','offshell','total']}
        NPDF['u'][_]['h'] = {_: [] for _ in ['onshell','offshell','total']}
        NPDF['u'][_]['t'] = {_: [] for _ in ['onshell','offshell','total']}
        NPDF['d'][_]['d'] = {_: [] for _ in ['onshell','offshell','total']}
        NPDF['d'][_]['h'] = {_: [] for _ in ['onshell','offshell','total']}
        NPDF['d'][_]['t'] = {_: [] for _ in ['onshell','offshell','total']}
    NPDF['X'] = X

    grid = {}
    for tar in ['d','h','t']:
        idis.load_grid_custom(X,Q2,tar)
        idx = 'tar=%s,x=%s,Q2=%s'%(tar,X,Q2)
        grid[tar] = idis.grids[None][idx][tar]

    cnt = 0
    for par in replicas:
        lprint('Generating nuclear PDFs %s/%s'%(cnt+1,len(replicas)))
        parman.set_new_params(par)

        pdf.evolve(q2)
        off.evolve(q2)
        for nucleon in nuclei:
            for nucleus in ['d','h','t']:

                u_on  = pdf.storage[q2]['u']
                d_on  = pdf.storage[q2]['d']
                u_off = off.storage[q2]['u']
                d_off = off.storage[q2]['d']

                if nucleus in ['d']:
                    G_on  = grid[nucleus]['onshell'] ['f22']
                    G_off = grid[nucleus]['offshell']['f22']
                if nucleus in ['h','t']:
                    G_on  = grid[nucleus]['onshell'] ['f22%s'%nucleon]
                    G_off = grid[nucleus]['offshell']['f22%s'%nucleon]

                #--switch onshell pieces
                if nucleon=='n': d_on, u_on = u_on, d_on

                #--switch offshell pieces
                q = [0, u_off, d_off]
                q = offshell_model.get_model(q,nucleon,nucleus)
                u_off, d_off = q[1],q[2]

                Nu_on  = np.einsum('xi,i->x',G_on ,u_on )
                Nd_on  = np.einsum('xi,i->x',G_on ,d_on )
                Nu_off = np.einsum('xi,i->x',G_off,u_off)
                Nd_off = np.einsum('xi,i->x',G_off,d_off)

                #--mellin inversion
                #--multiply by x here as well
                phase = conf['mellin'].phase
                Nu_on  = X*np.imag(phase*Nu_on )/np.pi 
                Nd_on  = X*np.imag(phase*Nd_on )/np.pi 
                Nu_off = X*np.imag(phase*Nu_off)/np.pi 
                Nd_off = X*np.imag(phase*Nd_off)/np.pi 

                NPDF['u'][nucleon][nucleus]['onshell'] .append(Nu_on)
                NPDF['u'][nucleon][nucleus]['offshell'].append(Nu_off)
                NPDF['u'][nucleon][nucleus]['total']   .append(Nu_on+Nu_off)
                NPDF['d'][nucleon][nucleus]['onshell'] .append(Nd_on)
                NPDF['d'][nucleon][nucleus]['offshell'].append(Nd_off)
                NPDF['d'][nucleon][nucleus]['total']   .append(Nd_on+Nd_off)

        cnt +=1
     
    filename = '%s/data/nuclear-pdf-Q2=%3.5f.dat'%(wdir,Q2[0])
    save(NPDF,filename) 
    print()
    print('Saving nuclear PDF data to %s'%filename)

def plot_nuclear_pdf(wdir,Q2,mode=1):

    nrows,ncols=3,2
    fig = py.figure(figsize=(ncols*8,nrows*5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax21 = py.subplot(nrows,ncols,3)
    ax22 = py.subplot(nrows,ncols,4)
    ax31 = py.subplot(nrows,ncols,5)
    ax32 = py.subplot(nrows,ncols,6)

    hand = {}
    thy  = {}
    replicas = core.get_replicas(wdir)

    load_config('%s/input.py'%wdir)

    istep = core.get_istep()
    core.mod_conf(istep,replicas[0])

    replicas = core.get_replicas(wdir)
    names    = core.get_replicas_names(wdir)

    jar = load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas = jar['replicas']
   
    if 'off pdf' not in conf['steps'][istep]['active distributions']: 
        if 'off pdf' not in conf['steps'][istep]['passive distributions']: 
            return

    #--try to load data.  Generate it if it does not exist
    try:
        data = load('%s/data/nuclear-pdf-Q2=%3.5f.dat'%(wdir,Q2))
    except:
        gen_nuclear_pdf(wdir,Q2)
        data = load('%s/data/nuclear-pdf-Q2=%3.5f.dat'%(wdir,Q2))

    pdfs = load('%s/data/pdf-Q2=%3.5f.dat'%(wdir,Q2))

    ##############################################
    #--plot offshell
    ##############################################
    X   = np.array(data['X'])

    for pdf in data:
        if pdf=='X': continue
        for nucleon in data[pdf]:
            for nucleus in data[pdf][nucleon]:

                if pdf=='u':
                    if nucleon=='p': ax,axoff,axrat = ax11,ax21,ax31
                    else:            continue
                if pdf=='d':
                    if nucleon=='p': ax,axoff,axrat = ax12,ax22,ax32
                    else:            continue

                if nucleus=='d': color,alpha='darkblue' ,0.7
                if nucleus=='h': color,alpha='magenta'    ,0.7
                if nucleus=='t': color,alpha='goldenrod',0.4

                qon  = np.array(data[pdf][nucleon][nucleus]['onshell'])
                qoff = np.array(data[pdf][nucleon][nucleus]['offshell'])
                q    = qon + qoff
                qrat = qoff/q

                meanon  = np.mean(qon,axis=0)
                stdon   = np.std (qon,axis=0)               
                meanoff = np.mean(qoff,axis=0)
                stdoff  = np.std (qoff,axis=0)               
                mean    = np.mean(q,axis=0)
                std     = np.std (q,axis=0)               
                meanrat = np.mean(qrat,axis=0)
                stdrat  = np.std (qrat,axis=0)               

                zorder = 1
                style = '-'
                if mode == 0:
                    for i in range(len(q)):
                        hand[nucleus] ,= ax.plot(X,q[i]   ,color=color,alpha=alpha,zorder=zorder)
                        axoff              .plot(X,qoff[i],color=color,alpha=alpha,zorder=zorder)
                        axrat              .plot(X,qrat[i],color=color,alpha=alpha,zorder=zorder)
                if mode == 1:
                    hand[nucleus] = ax.fill_between(X,mean   -std   ,mean   +std   ,color=color,alpha=alpha,zorder=zorder,hatch=style)
                    axoff             .fill_between(X,meanoff-stdoff,meanoff+stdoff,color=color,alpha=alpha,zorder=zorder,hatch=style)
                    axrat             .fill_between(X,meanrat-stdrat,meanrat+stdrat,color=color,alpha=alpha,zorder=zorder,hatch=style)
 
    u = np.array(pdfs['XF']['u'])
    d = np.array(pdfs['XF']['d'])
    meanu = np.mean(u,axis=0)
    stdu  = np.std (u,axis=0)
    meand = np.mean(d,axis=0)
    stdd  = np.std (d,axis=0)
    color, alpha = 'darkgreen',0.5
    if mode==0:
        for i in range(len(u)):
            hand['free'] ,= ax11.plot(X,u[i],color=color,alpha=alpha,zorder=zorder)
            ax12                .plot(X,d[i],color=color,alpha=alpha,zorder=zorder)
    if mode==1:
        hand['free'] = ax11.fill_between(X,meanu-stdu,meanu+stdu,color=color,alpha=alpha,zorder=zorder,hatch=style)
        ax12               .fill_between(X,meand-stdd,meand+stdd,color=color,alpha=alpha,zorder=zorder,hatch=style)
    ax21.axhline(0,0,1,color=color,alpha=alpha)
    ax22.axhline(0,0,1,color=color,alpha=alpha)
    ax31.axhline(0,0,1,color=color,alpha=alpha)
    ax32.axhline(0,0,1,color=color,alpha=alpha)


 
    ##############################################

    ax12.text(0.60,0.30,r'$Q^2=%s{\rm~GeV^2}$'%Q2,size=30,transform=ax12.transAxes)

    ax11.text(0.05,0.05,r'\boldmath$x u^{p/A}$' ,transform=ax11.transAxes,size=40)
    ax12.text(0.05,0.05,r'\boldmath$x d^{p/A}$' ,transform=ax12.transAxes,size=40)
    ax21.text(0.60,0.80,r'\boldmath$x u^{p/A ({\rm off)}}$',transform=ax21.transAxes,size=40)
    ax22.text(0.60,0.80,r'\boldmath$x d^{p/A ({\rm off)}}$',transform=ax22.transAxes,size=40)
    ax31.text(0.05,0.80,r'\boldmath$  u^{p/A ({\rm off)}}/u^{p/A}$',transform=ax31.transAxes,size=40)
    ax32.text(0.05,0.80,r'\boldmath$  d^{p/A ({\rm off)}}/d^{p/A}$',transform=ax32.transAxes,size=40)
 
    for ax in [ax11,ax12,ax21,ax22,ax31,ax32]:
        ax.set_xlim(0.01,0.90)
        ax.set_xlabel(r'\boldmath$x$'         ,size=30)
        ax.xaxis.set_label_coords(0.98,0.00)
        minorLocator = MultipleLocator(0.05)
        majorLocator = MultipleLocator(0.2)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.set_xticks([0.2,0.4,0.6,0.8])

    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30,labelbottom=False)

    for ax in [ax21,ax22]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30,labelbottom=False)
        ax.axhline(0,0,1,color='black',ls=':',alpha=0.5)

    for ax in [ax31,ax32]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.axhline(0,0,1,color='black',ls=':',alpha=0.5)

    ax11.set_ylim(0,0.70)
    ax12.set_ylim(0,0.45)
    ax21.set_ylim(-0.015,0.022)
    ax22.set_ylim(-0.015,0.022)
    ax31.set_ylim(-0.04,0.12)
    ax32.set_ylim(-0.28,0.18)

    minorLocator = MultipleLocator(0.05)
    majorLocator = MultipleLocator(0.2)
    ax11.yaxis.set_minor_locator(minorLocator)
    ax11.yaxis.set_major_locator(majorLocator)
    ax11.set_yticks([0.2,0.4,0.6])

    minorLocator = MultipleLocator(0.05)
    majorLocator = MultipleLocator(0.2)
    ax12.yaxis.set_minor_locator(minorLocator)
    ax12.yaxis.set_major_locator(majorLocator)
    ax12.set_yticks([0.2,0.4])

    minorLocator = MultipleLocator(0.0025)
    majorLocator = MultipleLocator(0.01)
    ax21.yaxis.set_minor_locator(minorLocator)
    ax21.yaxis.set_major_locator(majorLocator)

    minorLocator = MultipleLocator(0.0025)
    majorLocator = MultipleLocator(0.01)
    ax22.yaxis.set_minor_locator(minorLocator)
    ax22.yaxis.set_major_locator(majorLocator)

    minorLocator = MultipleLocator(0.01)
    majorLocator = MultipleLocator(0.05)
    ax31.yaxis.set_minor_locator(minorLocator)
    ax31.yaxis.set_major_locator(majorLocator)

    minorLocator = MultipleLocator(0.02)
    majorLocator = MultipleLocator(0.10)
    ax32.yaxis.set_minor_locator(minorLocator)
    ax32.yaxis.set_major_locator(majorLocator)


    handles,labels = [],[]
    handles.append(hand['free'])
    handles.append(hand['d'])
    handles.append(hand['h'])
    handles.append(hand['t'])
    labels.append(r'\textrm{\textbf{Free}}')
    labels.append(r'\boldmath$d$')
    labels.append(r'\boldmath$^3 {\rm He}$')
    labels.append(r'\boldmath$^3 {\rm H}$')
    ax12.legend(handles,labels,frameon=False,loc='upper right',fontsize=28, handletextpad = 0.5, handlelength = 1.5, ncol = 1, columnspacing = 0.5)


    py.tight_layout()
    py.subplots_adjust(hspace=0)


    filename = '%s/gallery/nuclear-pdfs'%wdir
    if mode==1: filename+='-bands'
    filename += '.png'
    print('Saving figures to %s'%filename)
    py.savefig(filename)
    py.clf()

def plot_delta3(wdir,Q2,mode=1):

    nrows,ncols=1,2
    fig = py.figure(figsize=(ncols*8,nrows*5))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    hand = {}
    thy  = {}
    replicas = core.get_replicas(wdir)

    load_config('%s/input.py'%wdir)

    istep = core.get_istep()
    core.mod_conf(istep,replicas[0])

    replicas = core.get_replicas(wdir)
    names    = core.get_replicas_names(wdir)

    jar = load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas = jar['replicas']
    
    if 'off pdf' not in conf['steps'][istep]['active distributions']: 
        if 'off pdf' not in conf['steps'][istep]['passive distributions']: 
            return

    #--try to load data.  Generate it if it does not exist
    try:
        data = load('%s/data/nuclear-pdf-Q2=%3.5f.dat'%(wdir,Q2))
    except:
        gen_nuclear_pdf(wdir,q2)
        data = load('%s/data/nuclear-pdf-Q2=%3.5f.dat'%(wdir,Q2))


    ##############################################
    #--plot offshell
    ##############################################
    X   = np.array(data['X'])

    upT = np.array(data['u']['p']['t']['total'])
    upH = np.array(data['u']['p']['h']['total'])
    dpT = np.array(data['d']['p']['t']['total'])
    dpH = np.array(data['d']['p']['h']['total'])

    urat = (upT - upH)/(upT + upH) 
    drat = (dpT - dpH)/(dpT + dpH) 

    #new = []
    #for i in range(len(urat)):
    #    val = urat[i][-20] #--x = 0.81
    #    if val < 0: continue
    #    new.append(urat[i])

    #urat = np.array(new)

    meanu = np.mean(urat,axis=0)
    stdu  = np.std (urat,axis=0)
    meand = np.mean(drat,axis=0)
    stdd  = np.std (drat,axis=0)
    
    color = 'firebrick'
    alpha = 0.8

    zorder = 1

    if mode==0:
        for i in range(len(urat)): 
            ax11.plot(X,urat[i],color=color,alpha=alpha,zorder=zorder,ls=style)
        for i in range(len(drat)): 
            ax12.plot(X,drat[i],color=color,alpha=alpha,zorder=zorder,ls=style)
    if mode==1: 
        ax11.fill_between(X,meanu-stdu,meanu+stdu,color=color,alpha=alpha,zorder=zorder)
        ax12.fill_between(X,meand-stdd,meand+stdd,color=color,alpha=alpha,zorder=zorder)
  
 
    ##############################################

    ax12.text(0.05,0.05,r'$Q^2=%s{\rm~GeV^2}$'%Q2,size=30,transform=ax11.transAxes)

    #ax11.text(0.05,0.80,r'\boldmath$(u^{p/^3 {\rm H}}-u^{p/^3 {\rm He}})/(u^{p/^3 {\rm H}}+u^{p/^3 {\rm He}})$' ,transform=ax11.transAxes,size=40)
    #ax12.text(0.05,0.80,r'\boldmath$(d^{p/^3 {\rm H}}-d^{p/^3 {\rm He}})/(d^{p/^3 {\rm H}}+d^{p/^3 {\rm He}})$' ,transform=ax12.transAxes,size=40)
    ax11.text(0.02,0.82,r'\boldmath$\Delta_3^u$' ,transform=ax11.transAxes,size=50)
    ax12.text(0.02,0.82,r'\boldmath$\Delta_3^d$' ,transform=ax12.transAxes,size=50)
 
    for ax in [ax11,ax12]:
        ax.set_xlim(0.01,0.90)
        ax.set_xlabel(r'\boldmath$x$',size=40)
        ax.xaxis.set_label_coords(0.98,0.00)
        ax.xaxis.set_minor_locator(MultipleLocator(0.05))
        ax.xaxis.set_major_locator(MultipleLocator(0.2))
        ax.set_xticks([0.2,0.4,0.6,0.8])
        ax.set_ylim(-0.15,0.14)

    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.axhline(0,0,1,color='black',ls=':',alpha=0.5)

    ax12.tick_params(labelleft=False)

    #handles,labels = [],[]
    #handles.append(hand['d'])
    #handles.append(hand['h'])
    #handles.append(hand['t'])
    #labels.append(r'\boldmath$d$')
    #labels.append(r'\boldmath$^3 {\rm He}$')
    #labels.append(r'\boldmath$^3 {\rm H}$')
    #ax12.legend(handles,labels,frameon=False,loc='upper left',fontsize=28, handletextpad = 0.5, handlelength = 1.5, ncol = 1, columnspacing = 0.5)


    #py.tight_layout()
    py.subplots_adjust(left=0.08,top=0.99,right=0.99,wspace=0.01)

    filename = '%s/gallery/nuclear-pdfs-delta3'%wdir
    if mode==1: filename += '-bands'
    filename += '.png'
    print('Saving figures to %s'%filename)
    py.savefig(filename)
    py.clf()















