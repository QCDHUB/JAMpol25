#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT

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

#--from obslib
from obslib.pjet.reader import READER

import kmeanconf as kc

#ext = '.png'
ext = '.pdf'

cwd = os.getcwd()
checkdir('gallery')

def plot_pjets(wdir,kc):

    print('\ngenerating polarized jet plots from %s'%(wdir))
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))

    nrows,ncols=3,3
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax13 = py.subplot(nrows,ncols,3)
    ax21 = py.subplot(nrows,ncols,4)
    ax22 = py.subplot(nrows,ncols,5)
    ax23 = py.subplot(nrows,ncols,6)
    ax31 = py.subplot(nrows,ncols,7)
    ax32 = py.subplot(nrows,ncols,8)
    ax33 = py.subplot(nrows,ncols,9)

    filters = conf['datasets']['pjet']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['pjet']={}
    conf['datasets']['pjet']['filters']=filters
    conf['datasets']['pjet']['xlsx']={}
    conf['datasets']['pjet']['xlsx'][20001]='pjet/expdata/20001.xlsx' #STAR 2003
    conf['datasets']['pjet']['xlsx'][20002]='pjet/expdata/20002.xlsx' #STAR 2005
    conf['datasets']['pjet']['xlsx'][20003]='pjet/expdata/20003.xlsx' #STAR 2006
    conf['datasets']['pjet']['xlsx'][20004]='pjet/expdata/20004.xlsx' #STAR 2009
    conf['datasets']['pjet']['xlsx'][20005]='pjet/expdata/20005.xlsx' #PHENIX 2005
    conf['datasets']['pjet']['xlsx'][20006]='pjet/expdata/20006.xlsx' #STAR 2012
    conf['datasets']['pjet']['xlsx'][20007]='pjet/expdata/20007.xlsx' #STAR 2015
    conf['datasets']['pjet']['xlsx'][20008]='pjet/expdata/20008.xlsx' #STAR 2013
    conf['datasets']['pjet']['norm']={}
    conf['pjet tabs']=READER().load_data_sets('pjet')
    tabs = conf['pjet tabs']
   
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   
   
    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']

    data = predictions['reactions']['pjet']

    #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc)

    #--get theory by seperating solutions and taking mean
    for idx in tabs:
        predictions = copy.copy(data[idx]['prediction-rep'])
        predictions_ic = [predictions[i] for i in range(len(predictions))]
        data[idx]['thy']  = np.mean(predictions_ic,axis=0)
        data[idx]['dthy'] = np.std(predictions_ic,axis=0)

    #######################
    #--plot absolute values
    #######################

    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data
    for idx in tabs:
        pT = data[idx]['pT']
        values = data[idx]['value']
        alpha  = data[idx]['alpha']
        if idx==20001: ax = ax11
        if idx==20002: ax = ax12
        if idx==20003: ax = ax13
        if idx==20005: continue
        if idx==20006: ax = ax23
        if idx==20008: ax = ax31

        #--PHENIX in blue, STAR in black
        if idx in [20005]: color ='black' 
        else:              color ='black'

        ms = 5
        ts = 30 #text size for eta limits
      
        if idx in [20004]:
            idx1 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.0]
            idx2 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.5]
            hand[idx] = ax21.errorbar(pT[idx1],values[idx1],yerr=alpha[idx1],color=color,fmt='o',ms=ms,capsize=3.0)
            hand[idx] = ax22.errorbar(pT[idx2],values[idx2],yerr=alpha[idx2],color=color,fmt='o',ms=ms,capsize=3.0)
        elif idx in [20007]:
            idx1 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.0]
            idx2 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.5]
            hand[idx] = ax32.errorbar(pT[idx1],values[idx1],yerr=alpha[idx1],color=color,fmt='o',ms=ms,capsize=3.0)
            hand[idx] = ax33.errorbar(pT[idx2],values[idx2],yerr=alpha[idx2],color=color,fmt='o',ms=ms,capsize=3.0)
        else:
            hand[idx] = ax.errorbar(pT,values,yerr=alpha,color=color,fmt='o',ms=ms,capsize=3.0)


        xpos = 0.57
        ypos = 0.07
        if idx in [20004]:
            etamin = data[idx]['eta-abs-min'][idx1][0]
            etamax = data[idx]['eta-abs-max'][idx1][0]
            ax21.text(xpos, ypos, r'$|\eta_{\rm jet}| \in [%2.1f,%2.1f]$'%(etamin,etamax)  ,transform=ax21.transAxes,size=ts)
            etamin = data[idx]['eta-abs-min'][idx2][0]
            etamax = data[idx]['eta-abs-max'][idx2][0]
            ax22.text(xpos, ypos, r'$|\eta_{\rm jet}| \in [%2.1f,%2.1f]$'%(etamin,etamax)  ,transform=ax22.transAxes,size=ts)
        elif idx in [20007]:
            etamin = data[idx]['eta-abs-min'][idx1][0]
            etamax = data[idx]['eta-abs-max'][idx1][0]
            ax32.text(xpos, ypos, r'$|\eta_{\rm jet}| \in [%2.1f,%2.1f]$'%(etamin,etamax)  ,transform=ax32.transAxes,size=ts)
            etamin = data[idx]['eta-abs-min'][idx2][0]
            etamax = data[idx]['eta-abs-max'][idx2][0]
            ax33.text(xpos, ypos, r'$|\eta_{\rm jet}| \in [%2.1f,%2.1f]$'%(etamin,etamax)  ,transform=ax33.transAxes,size=ts)
        else:
            if 'eta-abs-min' in data[idx]: 
                etamin = data[idx]['eta-abs-min'][0]
                etamax = data[idx]['eta-abs-max'][0]
                ax.text(xpos, ypos, r'$|\eta_{\rm jet}| \in [%2.1f,%2.1f]$'%(etamin,etamax)  ,transform=ax.transAxes,size=ts)
            if 'eta-min'     in data[idx]: 
                etamin = data[idx]['eta-min'][0]
                etamax = data[idx]['eta-max'][0]
                ax.text(xpos, ypos, r'$\eta_{\rm jet}   \in [%2.1f,%2.1f]$'%(etamin,etamax)  ,transform=ax.transAxes,size=ts)




        thy = data[idx]['thy']
        std = data[idx]['dthy']
        down = thy - std
        up   = thy + std

        #--PHENIX in blue, STAR in red
        if idx in [20005]: color = 'firebrick'
        else:              color = 'firebrick'

        if idx in [20004]:
            idx1 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.0]
            idx2 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.5]
            thy_plot[idx] ,= ax21.plot(pT[idx1],thy[idx1],color=color)
            thy_plot[idx] ,= ax22.plot(pT[idx2],thy[idx2],color=color)
            thy_band[idx]  = ax21.fill_between(pT[idx1],down[idx1],up[idx1],color=color,alpha=0.5)
            thy_band[idx]  = ax22.fill_between(pT[idx2],down[idx2],up[idx2],color=color,alpha=0.5)
        elif idx in [20007]:
            idx1 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.0]
            idx2 = [i for i in range(len(pT)) if data[idx]['eta-abs-min'][i]==0.5]
            thy_plot[idx] ,= ax32.plot(pT[idx1],thy[idx1],color=color)
            thy_plot[idx] ,= ax33.plot(pT[idx2],thy[idx2],color=color)
            thy_band[idx]  = ax32.fill_between(pT[idx1],down[idx1],up[idx1],color=color,alpha=0.5)
            thy_band[idx]  = ax33.fill_between(pT[idx2],down[idx2],up[idx2],color=color,alpha=0.5)
        else:
            thy_plot[idx] ,= ax.plot(pT,thy,color=color)
            thy_band[idx]  = ax.fill_between(pT,down,up,color=color,alpha=0.5)

    for ax in [ax11,ax12,ax13,ax21,ax22,ax23]:
        ax.tick_params(axis='both',which='both',top=True,right=True,labelbottom=False,direction='in',labelsize=35)
        #ax.set_yticks([0.25,0.5,0.75])

    for ax in [ax31,ax32,ax33]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=35)
        #ax.set_yticks([0.25,0.5,0.75])

    for ax in [ax11,ax12,ax13,ax21,ax22,ax23,ax31,ax32,ax33]:
        ax.set_xlim(5,65)
        ax.set_xticks([10,20,30,40,50,60])
        minorLocator = MultipleLocator(2)
        majorLocator = MultipleLocator(10)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,color='black',alpha=0.2)

    #for ax in [ax31,ax32]:
    #    ax.set_xlabel(r'\boldmath$p_T^{\rm jet}$',size=35)
    #    #ax.xaxis.set_label_coords(0.90,0.00)
    for ax in [ax31,ax32,ax33]:
        ax.set_xlabel(r'\boldmath$p_T^{\rm jet}~[{\rm GeV}]$',size=35)
        # ax.xaxis.set_label_coords(0.85,-0.01)

    for ax in [ax12,ax13,ax22,ax23,ax32,ax33]:
        ax.tick_params(axis='both',which='both',labelleft=False)

    for ax in [ax11,ax12,ax13]:
        ax.set_ylim(-0.1,0.33)
        minorLocator = MultipleLocator(0.05)
        majorLocator = MultipleLocator(0.10)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_major_locator(majorLocator)
        ax.set_yticks([0,0.1,0.2,0.3])
        ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$'])

    for ax in [ax21,ax22,ax23]:
        ax.set_ylim(-0.05,0.08)
        minorLocator = MultipleLocator(0.01)
        majorLocator = MultipleLocator(0.02)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_major_locator(majorLocator)
        ax.set_yticks([-0.04,-0.02,0,0.02,0.04,0.06])
        ax.set_yticklabels([r'$-0.04$',r'$-0.02$',r'$0$',r'$0.02$',r'$0.04$',r'$0.06$'])

    for ax in [ax31,ax32,ax33]:
        ax.set_ylim(-0.03,0.08)
        minorLocator = MultipleLocator(0.01)
        majorLocator = MultipleLocator(0.02)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_major_locator(majorLocator)
        ax.set_yticks([-0.02,0,0.02,0.04,0.06])
        ax.set_yticklabels([r'$-0.02$',r'$0$',r'$0.02$',r'$0.04$',r'$0.06$'])




    ax11.text(0.05, 0.70, r'\boldmath$A_{LL}^{\rm jet}$'  ,transform=ax11.transAxes,size=60)

    ax11.text(0.78, 0.85, r'\textrm{\textbf{2003}}'  ,transform=ax11.transAxes,size=40)
    ax12.text(0.78, 0.85, r'\textrm{\textbf{2005}}'  ,transform=ax12.transAxes,size=40)
    ax13.text(0.78, 0.85, r'\textrm{\textbf{2006}}'  ,transform=ax13.transAxes,size=40)
    ax21.text(0.78, 0.85, r'\textrm{\textbf{2009}}'  ,transform=ax21.transAxes,size=40)
    ax22.text(0.78, 0.85, r'\textrm{\textbf{2009}}'  ,transform=ax22.transAxes,size=40)
    ax23.text(0.78, 0.85, r'\textrm{\textbf{2012}}'  ,transform=ax23.transAxes,size=40)
    ax31.text(0.78, 0.85, r'\textrm{\textbf{2013}}'  ,transform=ax31.transAxes,size=40)
    ax32.text(0.78, 0.85, r'\textrm{\textbf{2015}}'  ,transform=ax32.transAxes,size=40)
    ax33.text(0.78, 0.85, r'\textrm{\textbf{2015}}'  ,transform=ax33.transAxes,size=40)

    handles, labels = [],[]
    handles.append((thy_plot[20001],thy_band[20001],hand[20001]))

    labels.append(r'\textbf{\textrm{STAR}}')

    ax13.legend(handles,labels,frameon=False,fontsize=40,loc='upper left',handletextpad = 0.5, handlelength = 1.0, columnspacing = 0.5)
    

    #handles, labels = [],[]
    #handles.append(thy_band[20004])
    #handles.append(thy_band[20017])
    #labels.append(r'\textbf{\textrm{JAM (HERMES)}}')
    #labels.append(r'\textbf{\textrm{JAM (COMPASS)}}')

    #ax12.legend(handles,labels,frameon=False,fontsize=22,loc=(0.0,0.35),handletextpad = 0.5, handlelength = 1.0)

    py.tight_layout()
    py.subplots_adjust(hspace=0.02,wspace=0.02,right=0.99,top=0.99)

    checkdir('%s/gallery'%wdir)
    filename='%s/gallery/fig_pjets'%(cwd)
    filename += ext

    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    ax13.set_rasterized(True)
    ax21.set_rasterized(True)
    ax22.set_rasterized(True)
    ax23.set_rasterized(True)
    ax31.set_rasterized(True)
    ax32.set_rasterized(True)
    ax33.set_rasterized(True)


    py.savefig(filename)
    print()
    print('Saving polarized jet plot to %s'%filename)


if __name__ == "__main__":

    wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g'

    plot_pjets(wdir, kc)




