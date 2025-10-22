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
from obslib.psidis.reader import READER

import kmeanconf as kc

#ext = '.png'
ext = '.pdf'

cwd = os.getcwd()

def plot_pion(wdir,kc):

    print('\ngenerating pion PSIDIS plots from %s'%(wdir))
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))

    nrows,ncols=2,2
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax21 = py.subplot(nrows,ncols,3)
    ax22 = py.subplot(nrows,ncols,4)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][20004]='psidis/expdata/20004.xlsx'
    conf['datasets']['psidis']['xlsx'][20005]='psidis/expdata/20005.xlsx'
    conf['datasets']['psidis']['xlsx'][20017]='psidis/expdata/20017.xlsx'
    conf['datasets']['psidis']['xlsx'][20018]='psidis/expdata/20018.xlsx'
    conf['datasets']['psidis']['xlsx'][20008]='psidis/expdata/20008.xlsx'
    conf['datasets']['psidis']['xlsx'][20009]='psidis/expdata/20009.xlsx'
    conf['datasets']['psidis']['xlsx'][20021]='psidis/expdata/20021.xlsx'
    conf['datasets']['psidis']['xlsx'][20022]='psidis/expdata/20022.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   
   
    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']

    data = predictions['reactions']['psidis']

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
        X = tabs[idx]['X']
        values = tabs[idx]['value']
        alpha  = data[idx]['alpha']
        if idx==20004: ax,color,fmt,ms = ax11,'firebrick','o',5
        if idx==20005: ax,color,fmt,ms = ax12,'firebrick','o',5 
        if idx==20017: ax,color,fmt,ms = ax11,'firebrick','^',5
        if idx==20018: ax,color,fmt,ms = ax12,'firebrick','^',5
        if idx==20008: ax,color,fmt,ms = ax21,'firebrick','o',5
        if idx==20009: ax,color,fmt,ms = ax22,'firebrick','o',5 
        if idx==20021: ax,color,fmt,ms = ax21,'firebrick','^',5
        if idx==20022: ax,color,fmt,ms = ax22,'firebrick','^',5

        hand[idx] = ax.errorbar(X,values,yerr=alpha,color=color,fmt=fmt,ms=ms,capsize=3.0)

        #if idx in [20004,20005,20008,20009]: continue

        thy = data[idx]['thy']
        std = data[idx]['dthy']
        down = thy - std
        up   = thy + std
        thy_plot[idx] ,= ax.plot(X,thy,color=color)
        thy_band[idx]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,labelbottom=False,direction='in',labelsize=30)
        ax.set_ylim(-0.05,1.0)
        ax.set_yticks([0,0.2,0.4,0.6,0.8])
        ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
        ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    for ax in [ax21,ax22]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_ylim(-0.05,0.29)
        ax.set_yticks([0,0.1,0.2])
        ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$'])
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    for ax in [ax11,ax12,ax21,ax22]:
        ax.semilogx()
        ax.set_xlim(0.01,0.6)
        ax.set_xticks([0.02,0.1,0.4])
        ax.set_xticklabels([r'$0.02$',r'$0.1$',r'$0.4$'])
        #minorLocator = MultipleLocator(0.1)
        #majorLocator = MultipleLocator(1.0)
        #ax.xaxis.set_minor_locator(minorLocator)
        #ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,color='black',alpha=0.2)

    for ax in [ax21,ax22]:
        # ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=40)
        ax.set_xlabel(r'\boldmath$x$',size=40)
        # ax.xaxis.set_label_coords(0.94,0.00)

    for ax in [ax12,ax22]:
        ax.tick_params(axis='both',which='both',labelleft=False)



    ax11.text(0.02, 0.75, r'\boldmath$A^{\pi^+}_{1}$',transform=ax11.transAxes,size=55)
    ax12.text(0.02, 0.75, r'\boldmath$A^{\pi^-}_{1}$',transform=ax12.transAxes,size=55)
    ax21.text(0.02, 0.75, r'\boldmath$A^{\pi^+}_{1}$',transform=ax21.transAxes,size=55)
    ax22.text(0.02, 0.75, r'\boldmath$A^{\pi^-}_{1}$',transform=ax22.transAxes,size=55)

    ax11.text(0.02,0.55,r'\rm \bf proton'  ,transform=ax11.transAxes,size=30)
    ax12.text(0.02,0.55,r'\rm \bf proton'  ,transform=ax12.transAxes,size=30)
    ax21.text(0.02,0.55,r'\rm \bf deuteron',transform=ax21.transAxes,size=30)
    ax22.text(0.02,0.55,r'\rm \bf deuteron',transform=ax22.transAxes,size=30)

    handles, labels = [],[]
    #handles.append((thy_plot[20004],thy_band[20004],hand[20004]))
    handles.append((thy_plot[20017],thy_band[20017],hand[20017]))

    #labels.append(r'\textbf{\textrm{HERMES}}')
    labels.append(r'\textbf{\textrm{COMPASS}}')

    #ax12.legend(handles,labels,frameon=False,fontsize=30,loc=(0.33,0.53),handletextpad = 0.5, handlelength = 1.0)
    ax12.legend(handles,labels,frameon=False,fontsize=30,loc=(0.50,0.75),handletextpad = 0.5, handlelength = 1.0)
    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%wdir)
    filename='%s/gallery/fig_psidis_pion'%(cwd)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving PSIDIS pion plot to %s'%filename)

def plot_kaon(wdir,kc):

    print('\ngenerating kaon PSIDIS plots from %s'%(wdir))
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))

    nrows,ncols=2,2
    fig = py.figure(figsize=(ncols*7,nrows*4))

    ax11 = py.subplot(nrows,ncols,1) # p, K+
    ax12 = py.subplot(nrows,ncols,2) # p, K-
    ax21 = py.subplot(nrows,ncols,3) # d, K+
    ax22 = py.subplot(nrows,ncols,4) # d, K-

    filters = conf['datasets']['psidis']['filters']
    
    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][20012]='psidis/expdata/20012.xlsx'
    conf['datasets']['psidis']['xlsx'][20013]='psidis/expdata/20013.xlsx'
    conf['datasets']['psidis']['xlsx'][20014]='psidis/expdata/20014.xlsx'
    conf['datasets']['psidis']['xlsx'][20019]='psidis/expdata/20019.xlsx'
    conf['datasets']['psidis']['xlsx'][20020]='psidis/expdata/20020.xlsx'
    conf['datasets']['psidis']['xlsx'][20025]='psidis/expdata/20025.xlsx'
    conf['datasets']['psidis']['xlsx'][20026]='psidis/expdata/20026.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   

    names = core.get_replicas_names(wdir)  

    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']

    data = predictions['reactions']['psidis']

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
        X = tabs[idx]['X']
        values = tabs[idx]['value']
        alpha  = data[idx]['alpha']
        if idx==20012: ax,color,fmt,ms = ax21,'firebrick','o',5
        if idx==20013: ax,color,fmt,ms = ax22,'firebrick','o',5 
        if idx==20014: ax,color,fmt,ms = ax23,'firebrick','o',5
        if idx==20019: ax,color,fmt,ms = ax11,'firebrick','^',5
        if idx==20020: ax,color,fmt,ms = ax12,'firebrick','^',5
        if idx==20025: ax,color,fmt,ms = ax21,'firebrick','^',5 
        if idx==20026: ax,color,fmt,ms = ax22,'firebrick','^',5

        hand[idx] = ax.errorbar(X,values,yerr=alpha,color=color,fmt=fmt,ms=ms,capsize=3.0)

        thy = data[idx]['thy']
        std = data[idx]['dthy']
        down = thy - std
        up   = thy + std
        thy_plot[idx] ,= ax.plot(X,thy,color=color)
        thy_band[idx]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,labelbottom=False,direction='in',labelsize=30)
        ax.set_ylim(-0.05,1.0)
        ax.set_yticks([0,0.2,0.4,0.6,0.8])
        ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])
        ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    for ax in [ax21,ax22]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_ylim(-0.19,0.3)
        ax.set_yticks([-0.1,0,0.1,0.2])
        ax.set_yticklabels([r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))

    for ax in [ax11,ax12,ax21,ax22]:
        ax.semilogx()
        ax.set_xlim(0.01,0.6)
        ax.set_xticks([0.02,0.1,0.4])
        ax.set_xticklabels([r'$0.02$',r'$0.1$',r'$0.4$'])
        #minorLocator = MultipleLocator(0.1)
        #majorLocator = MultipleLocator(1.0)
        #ax.xaxis.set_minor_locator(minorLocator)
        #ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,color='black',alpha=0.2)

    for ax in [ax21,ax22]:
        # ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=40)
        ax.set_xlabel(r'\boldmath$x$',size=40)
        # ax.xaxis.set_label_coords(0.94,0.00)

    for ax in [ax12,ax22]:
        ax.tick_params(axis='both',which='both',labelleft=False)


    ax11.text(0.02, 0.75, r'\boldmath$A^{K^+}_{1}$',transform=ax11.transAxes,size=55)
    ax12.text(0.02, 0.75, r'\boldmath$A^{K^-}_{1}$',transform=ax12.transAxes,size=55)
    ax21.text(0.02, 0.75, r'\boldmath$A^{K^+}_{1}$',transform=ax21.transAxes,size=55)
    ax22.text(0.02, 0.75, r'\boldmath$A^{K^-}_{1}$',transform=ax22.transAxes,size=55)

    ax11.text(0.02,0.55,r'\rm \bf proton'  ,transform=ax11.transAxes,size=30)
    ax12.text(0.02,0.55,r'\rm \bf proton'  ,transform=ax12.transAxes,size=30)
    ax21.text(0.02,0.55,r'\rm \bf deuteron',transform=ax21.transAxes,size=30)
    ax22.text(0.02,0.55,r'\rm \bf deuteron',transform=ax22.transAxes,size=30)

    handles, labels = [],[]
    #handles.append((thy_plot[20012],thy_band[20012],hand[20012]))
    handles.append((thy_plot[20019],thy_band[20019],hand[20019]))

    #labels.append(r'\textbf{\textrm{HERMES}}')
    labels.append(r'\textbf{\textrm{COMPASS}}')

    ax12.legend(handles,labels,frameon=False,fontsize=30,loc=(0.40,0.75),handletextpad = 0.5, handlelength = 1.0)
    
    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%wdir)
    filename='%s/gallery/fig_psidis_kaon'%(cwd)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving PSIDIS kaon plot to %s'%filename)

def plot_hadron(wdir,kc):

    print('\ngenerating hadron PSIDIS plots from %s'%(wdir))
    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))

    nrows,ncols=2,2
    fig = py.figure(figsize=(ncols*7,nrows*4))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)
    ax21 = py.subplot(nrows,ncols,3)
    ax22 = py.subplot(nrows,ncols,4)

    filters = conf['datasets']['psidis']['filters']

    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['psidis']={}
    conf['datasets']['psidis']['filters']=filters
    conf['datasets']['psidis']['xlsx']={}
    conf['datasets']['psidis']['xlsx'][20000]='psidis/expdata/20000.xlsx'
    conf['datasets']['psidis']['xlsx'][20001]='psidis/expdata/20001.xlsx'
    conf['datasets']['psidis']['xlsx'][20002]='psidis/expdata/20002.xlsx'
    conf['datasets']['psidis']['xlsx'][20003]='psidis/expdata/20003.xlsx'
    conf['datasets']['psidis']['xlsx'][20006]='psidis/expdata/20006.xlsx'
    conf['datasets']['psidis']['xlsx'][20007]='psidis/expdata/20007.xlsx'
    conf['datasets']['psidis']['xlsx'][20010]='psidis/expdata/20010.xlsx'
    conf['datasets']['psidis']['xlsx'][20011]='psidis/expdata/20011.xlsx'
    conf['datasets']['psidis']['xlsx'][20015]='psidis/expdata/20015.xlsx'
    conf['datasets']['psidis']['xlsx'][20016]='psidis/expdata/20016.xlsx'
    conf['datasets']['psidis']['xlsx'][20023]='psidis/expdata/20023.xlsx'
    conf['datasets']['psidis']['xlsx'][20024]='psidis/expdata/20024.xlsx'
    conf['datasets']['psidis']['norm']={}
    conf['psidis tabs']=READER().load_data_sets('psidis')
    tabs = conf['psidis tabs']
    
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   
   
    jar=load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas=jar['replicas']

    data = predictions['reactions']['psidis']

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
        X = tabs[idx]['X']
        values = tabs[idx]['value']
        alpha  = data[idx]['alpha']
        if idx==20000: ax,color,fmt,ms = ax11,'blue'     , '*', 5
        if idx==20001: ax,color,fmt,ms = ax12,'blue'     , '*', 5 
        if idx==20002: ax,color,fmt,ms = ax21,'blue'     , '*', 5
        if idx==20003: ax,color,fmt,ms = ax22,'blue'     , '*', 5
        if idx==20006: ax,color,fmt,ms = ax11,'firebrick', 'o', 5
        if idx==20007: ax,color,fmt,ms = ax12,'firebrick', 'o', 5 
        if idx==20010: ax,color,fmt,ms = ax21,'firebrick', 'o', 5
        if idx==20011: ax,color,fmt,ms = ax22,'firebrick', 'o', 5
        if idx==20015: ax,color,fmt,ms = ax31,'firebrick', 'o', 5
        if idx==20016: ax,color,fmt,ms = ax32,'firebrick', 'o', 5
        if idx==20023: ax,color,fmt,ms = ax21,'firebrick', '^', 5
        if idx==20024: ax,color,fmt,ms = ax22,'firebrick', '^', 5

        hand[idx] = ax.errorbar(X,values,yerr=alpha,color=color,fmt=fmt,ms=ms,capsize=3.0)

        thy = data[idx]['thy']
        std = data[idx]['dthy']
        down = thy - std
        up   = thy + std
        thy_plot[idx] ,= ax.plot(X,thy,color=color)
        thy_band[idx]  = ax.fill_between(X,down,up,color=color,alpha=0.5)

    for ax in [ax11,ax12]:
        ax.tick_params(axis='both',which='both',top=True,right=True,labelbottom=False,direction='in',labelsize=30)
        ax.set_ylim(-0.10,0.60)
        ax.set_yticks([0,0.2,0.4])
        ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$'])
        ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    for ax in [ax21,ax22]:
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_ylim(-0.10,0.70)
        ax.set_yticks([0,0.2,0.4,0.6])
        ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$'])
        ax.yaxis.set_minor_locator(MultipleLocator(0.10))

    for ax in [ax11,ax12,ax21,ax22]:
        ax.semilogx()
        ax.set_xlim(0.004,0.6)
        ax.set_xticks([0.01,0.1,0.4])
        ax.set_xticklabels([r'$0.01$',r'$0.1$',r'$0.4$'])
        #minorLocator = MultipleLocator(0.1)
        #majorLocator = MultipleLocator(1.0)
        #ax.xaxis.set_minor_locator(minorLocator)
        #ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,color='black',alpha=0.2)

    for ax in [ax21,ax22]:
        # ax.set_xlabel(r'\boldmath$x_{\rm bj}$',size=40)
        ax.set_xlabel(r'\boldmath$x$',size=40)
        # ax.xaxis.set_label_coords(0.94,0.00)

    for ax in [ax12,ax22]:
        ax.tick_params(axis='both',which='both',labelleft=False)



    ax11.text(0.02, 0.75, r'\boldmath$A^{h^+}_{1}$',transform=ax11.transAxes,size=55)
    ax12.text(0.02, 0.75, r'\boldmath$A^{h^-}_{1}$',transform=ax12.transAxes,size=55)
    ax21.text(0.02, 0.75, r'\boldmath$A^{h^+}_{1}$',transform=ax21.transAxes,size=55)
    ax22.text(0.02, 0.75, r'\boldmath$A^{h^-}_{1}$',transform=ax22.transAxes,size=55)

    ax11.text(0.02,0.55,r'\rm \bf proton'  ,transform=ax11.transAxes,size=30)
    ax12.text(0.02,0.55,r'\rm \bf proton'  ,transform=ax12.transAxes,size=30)
    ax21.text(0.02,0.55,r'\rm \bf deuteron',transform=ax21.transAxes,size=30)
    ax22.text(0.02,0.55,r'\rm \bf deuteron',transform=ax22.transAxes,size=30)

    handles, labels = [],[]
    #handles.append((thy_plot[20006],thy_band[20006],hand[20006]))
    handles.append((thy_plot[20023],thy_band[20023],hand[20023]))
    handles.append((thy_plot[20000],thy_band[20000],hand[20000]))

    #labels.append(r'\textbf{\textrm{HERMES}}')
    labels.append(r'\textbf{\textrm{COMPASS}}')
    labels.append(r'\textbf{\textrm{SMC}}')

    ax12.legend(handles,labels,frameon=False,fontsize=30,loc=(0.50,0.60),handletextpad = 0.5, handlelength = 1.0)
    


    py.tight_layout()
    py.subplots_adjust(hspace=0.01,wspace=0.01,top=0.99,right=0.99)

    checkdir('%s/gallery'%wdir)
    filename='%s/gallery/fig_psidis_hadron'%(cwd)
    filename += ext

    py.savefig(filename)
    print()
    print('Saving PSIDIS hadron plot to %s'%filename)


if __name__ == "__main__":

    wdir =  '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g'

    plot_pion  (wdir, kc)
    plot_kaon  (wdir, kc)
    plot_hadron(wdir, kc)




