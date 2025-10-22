#!/usr/bin/env python
import sys,os
import matplotlib
import matplotlib as plt
matplotlib.use('Agg')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT


#--matplotlib
from matplotlib.ticker import MultipleLocator
import pylab as py
import matplotlib.gridspec as gridspec
# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#-- from qcdlib
from qcdlib import aux

from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
# from analysis.corelib import classifier
import kmeanconf as kc

#--from obslib
from obslib.wzrv.theory import WZRV
from obslib.wzrv.reader import READER

#ext = '.png'
ext = '.pdf'


cwd = os.getcwd()
checkdir('gallery')

def plot_Wlep_backup(wdir):
    nrows,ncols=1,2
    fig = py.figure(figsize=(ncols*7,nrows*6))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    conf['path2wzrvtab'] = '%s/grids/grids-wzrv'%os.environ['FITPACK']
    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['wzrv']={}
    conf['datasets']['wzrv']['xlsx']={}
    conf['datasets']['wzrv']['xlsx'][1000]='wzrv/expdata/1000.xlsx'
    conf['datasets']['wzrv']['xlsx'][1020]='wzrv/expdata/1020.xlsx'
    conf['datasets']['wzrv']['xlsx'][1021]='wzrv/expdata/1021.xlsx'
    conf['datasets']['wzrv']['norm']={}
    conf['datasets']['wzrv']['filters']=[]
    conf['wzrv tabs']=READER().load_data_sets('wzrv')
    tabs = conf['wzrv tabs']
    #--plot data
    hand = {}
    for idx in tabs:
        if   idx==1000: ax,color = ax11,'black'
        elif idx==1020: ax,color = ax12,'blue'
        elif idx==1021: ax,color = ax12,'purple'
        values = tabs[idx]['value']
        alpha = np.sqrt(np.array(tabs[idx]['stat_u'])**2 + np.array(tabs[idx]['syst_u'])**2)
        if 'eta' in tabs[idx]:
            eta = tabs[idx]['eta']
        else: 
            eta = (tabs[idx]['eta_min'] + tabs[idx]['eta_max'])/2.0
            eta_min = tabs[idx]['eta_min']
            eta_max = tabs[idx]['eta_max']
        n = int(len(eta)/2)
        if   idx==1000: X = 0
        elif idx==1020: X =  0.08
        elif idx==1021: X =  0.08
        if idx in [1000]:
            hand[idx] = ax.errorbar(eta[:n]+X,values[:n],yerr=alpha[:n],color=color,fmt='o',ms=5.0,capsize=3.0)
            ax            .errorbar(eta[n:]-X,values[n:],yerr=alpha[n:],color=color,fmt='o',ms=5.0,capsize=3.0)
        elif idx in [1020,1021]:
            hand[idx] = {}
            hand[idx]['+'] = ax.errorbar(eta[:n]+X,values[:n],yerr=alpha[:n],color=color,fmt='o',ms=5.0,capsize=3.0)
            hand[idx]['-'] = ax.errorbar(eta[n:]-X,values[n:],yerr=alpha[n:],color=color,fmt='s',ms=6.0,capsize=3.0)

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   
   
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))
    data = predictions['reactions']['wzrv']
 
    #--get theory by seperating solutions and taking mean
    #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc)
    for idx in tabs:
        if   idx==1000: ax,color = ax11,'black'
        elif idx==1020: ax,color = ax12,'blue'
        elif idx==1021: ax,color = ax12,'purple'
        predictions = copy.copy(data[idx]['prediction-rep'])
        predictions_ic = [predictions[i] for i in range(len(predictions))]
        data[idx]['thy']  = np.mean(predictions_ic,axis=0)
        data[idx]['dthy'] = np.std(predictions_ic,axis=0)

        if 'eta' in tabs[idx]:
            eta = tabs[idx]['eta']
        else: 
            eta = (tabs[idx]['eta_min'] + tabs[idx]['eta_max'])/2.0
            eta_min = tabs[idx]['eta_min']
            eta_max = tabs[idx]['eta_max']
        n = int(len(eta)/2)

        #--compute cross-section for all replicas        
        thy  = data[idx]['thy']
        dthy = data[idx]['dthy']
        
        up   = thy + dthy 
        down = thy - dthy 

        if idx in [1000]:
            thy_plot0 ,= ax.plot(eta[:n],thy[:n],color='red',alpha=1.0,lw=2)
            ax             .plot(eta[n:],thy[n:],color='red',alpha=1.0,lw=2)
            thy_band0  = ax.fill_between(eta[:n],down[:n],up[:n],color='red',alpha=0.2)
            ax             .fill_between(eta[n:],down[n:],up[n:],color='red',alpha=0.2)
        elif idx in [1020,1021]:
            if idx==1020:   X =  0
            elif idx==1021: X =  0.12
            thy_point0 = ax.errorbar(eta[:n]+X,thy[:n],yerr=dthy[:n],color='red',fmt='o',ms=5.0,capsize=3.0)
            ax             .errorbar(eta[n:]-X,thy[n:],yerr=dthy[n:],color='red',fmt='s',ms=6.0,capsize=3.0)


    for ax in [ax11,ax12]:
        ax.set_xlim(-1.9,1.9)
        minorLocator = MultipleLocator(0.25)
        majorLocator = MultipleLocator(1.0)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        minorLocator = MultipleLocator(0.1)
        majorLocator = MultipleLocator(0.5)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,ls='-',color='black',alpha=0.2,lw=1.0)
        ax.set_xlabel(r'\boldmath$\eta$',size=40)
        # ax.xaxis.set_label_coords(0.95,-0.01)
        ax.set_xticks([-1.0,0.0,1.0])
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_ylim(-1.00,1.00)
        ax.set_yticks([-0.5,0,0.5])
        ax.set_yticklabels([r'$-0.5$',r'$0$',r'$0.5$'])

    ax12.tick_params(labelleft=False)

    ax11.text(0.03,0.07,r'\boldmath$A_L^{W}$'                                  ,transform=ax11.transAxes,size=50)
    ax11.text(0.20,0.08,r'\textrm{\textbf{(STAR)}}',transform=ax11.transAxes,size=30)
    ax11.text(0.70,0.43,r'$\sqrt{s}=510~$'+r'\textrm{GeV}'                 ,transform=ax11.transAxes,size=20)
    ax11.text(0.72,0.36,r'$p_T > 25~$'+r'\hspace{0.00cm}\textrm{GeV}' ,transform=ax11.transAxes,size=20)
    ax11.text(0.04,0.61,r'\boldmath$W^-$'                                  ,transform=ax11.transAxes,size=30)
    ax11.text(0.04,0.30,r'\boldmath$W^+$'                                  ,transform=ax11.transAxes,size=30)


    ax12.text(0.34,0.82,r'\boldmath$A_L^{W/Z}$'      ,transform=ax12.transAxes,size=50)
    ax12.text(0.62,0.84,r'\textrm{\textbf{(PHENIX)}}',transform=ax12.transAxes,size=30)
    #ax12.text(0.04,0.82,r'$W^-$'                                  ,transform=ax12.transAxes,size=25)
    #ax12.text(0.10,0.60,r'$W^+$'                                  ,transform=ax12.transAxes,size=25)
    #ax12.text(0.45,0.66,r'$W^-$'                                  ,transform=ax12.transAxes,size=25)
    #ax12.text(0.45,0.22,r'$W^+$'                                  ,transform=ax12.transAxes,size=25)
    #ax12.text(0.82,0.60,r'$W^-$'                                  ,transform=ax12.transAxes,size=25)
    #ax12.text(0.85,0.10,r'$W^+$'                                  ,transform=ax12.transAxes,size=25)

    handles, labels = [],[]
    handles.append((thy_band0,thy_plot0))
    #handles.append(hand[1000])

    labels.append(r'\textrm{\textbf{JAM}}')
    #labels.append(r'\textrm{\textbf{STAR}}')
    ax11.legend(handles,labels,frameon=False,fontsize=25,loc='upper left',handletextpad=0.3,handlelength=1.2,ncol=1,columnspacing=2.3)


    handles, labels = [],[]
    handles.append(hand[1020]['+'])
    handles.append(hand[1021]['+'])
    handles.append(hand[1020]['-'])
    handles.append(hand[1021]['-'])

    #labels.append(r'\textrm{\textbf{PHENIX (\boldmath$p_T^{\ell} > 30$)}}')
    #labels.append(r'\textrm{\textbf{PHENIX (\boldmath$p_T^{\ell} > 16$)}}')
    labels.append(r'\boldmath$W^+, p_T^{} > 16$')
    labels.append(r'\boldmath$W^+, p_T^{} > 30$')
    labels.append(r'\boldmath$W^-, p_T^{} > 16$')
    labels.append(r'\boldmath$W^-, p_T^{} > 30$')
    ax12.legend(handles,labels,frameon=False,fontsize=20,loc='lower left',handletextpad=0.3,handlelength=1.2,ncol=2,columnspacing=2.3,handleheight=1.5)


    #--custom legend
    ax12.errorbar(-1.68,-0.45,yerr=0.06,color='red',fmt='o',ms=5.0,capsize=3.0)
    ax12.errorbar(-1.61,-0.45,yerr=0.06,color='red',fmt='s',ms=6.0,capsize=3.0)
    ax12.text(0.10,0.25,r'\textrm{\textbf{JAM}}'      ,transform=ax12.transAxes,size=25)

  
    py.tight_layout()
    py.subplots_adjust(hspace=0,wspace=0.02,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/fig_polWlep'%(cwd)
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)

    py.savefig(filename)
    print('Saving polarized W-lepton plot to %s'%filename)

def plot_Wlep(wdir):
    nrows,ncols=1,2
    fig = py.figure(figsize=(ncols*7,nrows*6))
    ax11 = py.subplot(nrows,ncols,1)
    ax12 = py.subplot(nrows,ncols,2)

    conf['path2wzrvtab'] = '%s/grids/grids-wzrv'%os.environ['FITPACK']
    conf['aux']=aux.AUX()
    conf['datasets'] = {}
    conf['datasets']['wzrv']={}
    conf['datasets']['wzrv']['xlsx']={}
    conf['datasets']['wzrv']['xlsx'][1000]='wzrv/expdata/1000.xlsx'
    conf['datasets']['wzrv']['xlsx'][1020]='wzrv/expdata/1020.xlsx'
    conf['datasets']['wzrv']['xlsx'][1021]='wzrv/expdata/1021.xlsx'
    conf['datasets']['wzrv']['norm']={}
    conf['datasets']['wzrv']['filters']=[]
    conf['wzrv tabs']=READER().load_data_sets('wzrv')
    tabs = conf['wzrv tabs']
    #--plot data
    hand = {}
    for idx in tabs:
        if   idx==1000: ax,color = ax11,'firebrick'
        elif idx==1020: ax,color = ax12,'blue'
        elif idx==1021: ax,color = ax12,'darkgreen'
        values = tabs[idx]['value']
        alpha = np.sqrt(np.array(tabs[idx]['stat_u'])**2 + np.array(tabs[idx]['syst_u'])**2)
        if 'eta' in tabs[idx]:
            eta = tabs[idx]['eta']
        else: 
            eta = (tabs[idx]['eta_min'] + tabs[idx]['eta_max'])/2.0
            eta_min = tabs[idx]['eta_min']
            eta_max = tabs[idx]['eta_max']
        n = int(len(eta)/2)
        if   idx==1000: X = 0
        elif idx==1020: X =  0.08
        elif idx==1021: X =  0.08
        if idx in [1000]:
            hand[idx] = ax.errorbar(eta[:n]+X,values[:n],yerr=alpha[:n],color=color,fmt='o',ms=5.0,capsize=3.0)
            ax            .errorbar(eta[n:]-X,values[n:],yerr=alpha[n:],color=color,fmt='o',ms=5.0,capsize=3.0)
        elif idx in [1020,1021]:
            hand[idx] = {}
            hand[idx]['+'] = ax.errorbar(eta[:n]+X,values[:n],yerr=alpha[:n],color=color,fmt='o',ms=5.0,capsize=3.0)
            hand[idx]['-'] = ax.errorbar(eta[n:]-X,values[n:],yerr=alpha[n:],color=color,fmt='s',ms=6.0,capsize=3.0)

    load_config('%s/input.py'%wdir)
    istep=core.get_istep()
    replicas=core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0]) #--set conf as specified in istep   
   
    predictions = load('%s/data/predictions-%d.dat'%(wdir,istep))
    data = predictions['reactions']['wzrv']

    thy_plot, thy_band = {},{} 
    #--get theory by seperating solutions and taking mean
    #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc)
    for idx in tabs:
        if   idx==1000: ax,color = ax11,'firebrick'
        elif idx==1020: ax,color = ax12,'blue'
        elif idx==1021: ax,color = ax12,'darkgreen'
        predictions = copy.copy(data[idx]['prediction-rep'])
        predictions_ic = [predictions[i] for i in range(len(predictions))]
        data[idx]['thy']  = np.mean(predictions_ic,axis=0)
        data[idx]['dthy'] = np.std(predictions_ic,axis=0)

        if 'eta' in tabs[idx]:
            eta = tabs[idx]['eta']
        else: 
            eta = (tabs[idx]['eta_min'] + tabs[idx]['eta_max'])/2.0
            eta_min = tabs[idx]['eta_min']
            eta_max = tabs[idx]['eta_max']
        n = int(len(eta)/2)

        #--compute cross-section for all replicas        
        thy  = data[idx]['thy']
        dthy = data[idx]['dthy']
      
        boson = data[idx]['boson'] 
        up   = thy + dthy 
        down = thy - dthy 

        if   idx==1000: X = 0
        elif idx==1020: X =  0.08
        elif idx==1021: X =  0.08
        if idx in [1000]:
            thy_plot[idx] ,= ax.plot(eta[:n],thy[:n],color=color,alpha=1.0,lw=2)
            ax                 .plot(eta[n:],thy[n:],color=color,alpha=1.0,lw=2)
            thy_band[idx]  = ax.fill_between(eta[:n],down[:n],up[:n],color=color,alpha=0.2)
            ax                 .fill_between(eta[n:],down[n:],up[n:],color=color,alpha=0.2)
        elif idx in [1020,1021]:
            for i in range(len(eta)):
                if '-' in boson[i]: _eta = [eta[i]-0.06-X,eta[i]+0.06-X]
                if '+' in boson[i]: _eta = [eta[i]-0.06+X,eta[i]+0.06+X]
                _up,_down,_thy = [up[i],up[i]],[down[i],down[i]],[thy[i],thy[i]]
                thy_plot[idx] ,= ax.plot(_eta,_thy,color=color,alpha=1.0,lw=2)
                thy_band[idx]  = ax.fill_between(_eta,_down,_up,color=color,alpha=0.2)


    for ax in [ax11,ax12]:
        ax.set_xlim(-1.9,1.9)
        minorLocator = MultipleLocator(0.25)
        majorLocator = MultipleLocator(1.0)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_tick_params(which='major',length=6)
        ax.xaxis.set_tick_params(which='minor',length=3)
        minorLocator = MultipleLocator(0.1)
        majorLocator = MultipleLocator(0.5)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_major_locator(majorLocator)
        ax.yaxis.set_tick_params(which='major',length=6)
        ax.yaxis.set_tick_params(which='minor',length=3)
        ax.axhline(0,0,1,ls='-',color='black',alpha=0.2,lw=1.0)
        ax.set_xlabel(r'\boldmath$\eta$',size=40)
        # ax.xaxis.set_label_coords(0.95,-0.01)
        ax.set_xticks([-1.0,0.0,1.0])
        ax.tick_params(axis='both',which='both',top=True,right=True,direction='in',labelsize=30)
        ax.set_ylim(-1.00,1.00)
        #ax.set_yticks([-0.5,0,0.5])
        #ax.set_yticklabels([r'$-0.5$',r'$0$',r'$0.5$'])
        ax.set_yticks([-0.8,-0.4,0,0.4,0.8])
        ax.set_yticklabels([r'$-0.8$',r'$-0.4$',r'$0$',r'$0.4$',r'$0.8$'])

    ax12.tick_params(labelleft=False)

    ax11.text(0.03,0.85,r'\boldmath$A_L^{W}$'                                  ,transform=ax11.transAxes,size=50)
    ax11.text(0.71,0.93,r'$\sqrt{s}=510~$'+r'\textrm{GeV}'                 ,transform=ax11.transAxes,size=20)
    ax11.text(0.73,0.86,r'$p_T > 25~$'+r'\hspace{0.00cm}\textrm{GeV}' ,transform=ax11.transAxes,size=20)
    ax11.text(0.04,0.61,r'\boldmath$W^-$'                                  ,transform=ax11.transAxes,size=30)
    ax11.text(0.04,0.30,r'\boldmath$W^+$'                                  ,transform=ax11.transAxes,size=30)


    ax12.text(0.03,0.08,r'\boldmath$A_L^{W/Z}$'      ,transform=ax12.transAxes,size=50)
    ax12.text(0.20,0.90,r'\textrm{\textbf{PHENIX \boldmath$W^+$}}',transform=ax12.transAxes,size=30)
    ax12.text(0.72,0.90,r'\boldmath$W^-$',transform=ax12.transAxes,size=30)

    handles, labels = [],[]
    handles.append((thy_band[1000],thy_plot[1000],hand[1000]))
    #handles.append(hand[1000])

    labels.append(r'\textrm{\textbf{STAR}}')
    ax11.legend(handles,labels,frameon=False,fontsize=30,loc=(0.00,0.00),handletextpad=0.3,handlelength=1.2,ncol=1,columnspacing=2.3)


    handles, labels = [],[]
    handles.append((thy_band[1020],thy_plot[1020],hand[1020]['+']))
    handles.append((thy_band[1021],thy_plot[1021],hand[1021]['+']))
    handles.append((thy_band[1020],thy_plot[1020],hand[1020]['-']))
    handles.append((thy_band[1021],thy_plot[1021],hand[1021]['-']))

    labels.append(r'\boldmath$p_T > 16$')
    labels.append(r'\boldmath$p_T > 30$')
    labels.append(r'\boldmath$p_T > 16$')
    labels.append(r'\boldmath$p_T > 30$')
    ax12.legend(handles,labels,frameon=False,fontsize=20,loc=(0.18,0.70),handletextpad=0.3,handlelength=1.2,ncol=2,columnspacing=6.8,labelspacing=0.7)


  
    py.tight_layout()
    py.subplots_adjust(hspace=0,wspace=0.02,top=0.99,right=0.99)

    checkdir('%s/gallery'%cwd)
    filename='%s/gallery/fig_polWlep'%(cwd)
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)

    py.savefig(filename)
    print('Saving polarized W-lepton plot to %s'%filename)

if __name__=="__main__":

    wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g'

    plot_Wlep(wdir)






