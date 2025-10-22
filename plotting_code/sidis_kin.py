#!/usr/bin/env python

import sys,os
import numpy as np
import copy

path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

#--matplotlib
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import pylab as py


#--from tools
from tools.tools     import load,save,checkdir,lprint
from tools.config    import conf,load_config

#--from fitpack
from obslib.idis.reader    import READER as idisREAD
from obslib.pidis.reader   import READER as pidisREAD
from obslib.dy.reader      import READER as dyREAD
from obslib.wasym.reader   import READER as wasymREAD
from obslib.zrap.reader    import READER as zrapREAD
from obslib.wzrv.reader    import READER as wzrvREAD
from obslib.jet.reader    import READER as jetREAD
from obslib.pjet.reader   import READER as pjetREAD
from obslib.sidis.reader   import READER as sidisREAD
from obslib.psidis.reader  import READER as psidisREAD
from obslib.sia.reader     import READER as siaREAD
from qcdlib import aux

conf['aux']=aux.AUX()

ext = '.png'
#ext = '.pdf'

cwd = ''
#--make kinematic plot

def get_kin_X(exp,data):
    #--get X, Q2
    if exp == 'idis' or exp=='pidis' or exp=='sidis' or exp=='psidis':
        X  = data['X']
        Q2 = data['Q2']
    elif exp == 'dy' or exp == 'wasym' or exp == 'zrap':
        Y = data['Y']
        if exp == 'dy':
            Q2  = data['Q2']
            tau = data['tau']
        if exp == 'wasym':
            Q2  = 80.398**2*np.ones(len(Y))
            S   = data['S']
            tau = Q2/S
        if exp == 'zrap':
            Q2  = 91.1876**2*np.ones(len(Y))
            S   = data['S']
            tau = Q2/S
        X = np.sqrt(tau)*np.cosh(Y)
    elif exp == 'wzrv':
        if 'eta' in data:
            eta = data['eta']
        else:
            eta = (data['eta_min'] + data['eta_max'])/2.0
        Q2  = 80.398**2*np.ones(len(eta))
        S   = data['cms']**2
        tau = Q2/S
        X = np.sqrt(tau)*np.cosh(eta)
    elif exp == 'pjet' or exp=='jet':
        pT   = (data['pt-max'] + data['pt-min'])/2.0
        Q2   = pT**2
        S    = data['RS']*(10**3)
        tau  = Q2/S
        if 'eta-max' in data:
            etamin,etamax = data['eta-max'],data['eta-min']
            Xmin = 2*np.sqrt(tau)*np.cosh(etamin)
            Xmax = 2*np.sqrt(tau)*np.cosh(etamax)
            X = (Xmin+Xmax)/2.0
        else:
            eta = (data['eta-abs-max']+data['eta-abs-min'])/2.0
            X = 2*np.sqrt(tau)*np.cosh(eta)

    return X,Q2

def get_kin_Z(exp,data):
    #--get Z, Q2
    if  exp=='sidis' or exp=='psidis':
        Z  = data['Z']
        Q2 = data['Q2']

    elif exp=='sia':
        Z  = data['z']
        Q2 = data['RS']**2

    return Z,Q2

def load_data(Q2cut = 3.5, W2SIDIScut=20, zmincut = 0.2, zmaxcut = 0.8):

    conf['datasets'] = {}
    data = {}

    ##--SIDIS 
    conf['datasets']['sidis']={}
    conf['datasets']['sidis']['filters']=[]
    if Q2cut     != None: conf['datasets']['sidis']['filters'].append("Q2>%f"%Q2cut) 
    if W2SIDIScut!= None: conf['datasets']['sidis']['filters'].append("W2SIDIS>%f"%W2SIDIScut) 
    if zmincut   != None: conf['datasets']['sidis']['filters'].append('Z>%f'%(zmincut))
    if zmaxcut   != None: conf['datasets']['sidis']['filters'].append('Z<%f'%(zmaxcut))
    conf['datasets']['sidis']['xlsx']={}
    conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx' # deuteron , mult , pi+ , COMPASS
    conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx' # deuteron , mult , pi- , COMPASS
    #conf['datasets']['sidis']['xlsx'][2005]='sidis/expdata/2005.xlsx' # deuteron , mult , K+  , COMPASS
    #conf['datasets']['sidis']['xlsx'][2006]='sidis/expdata/2006.xlsx' # deuteron , mult , K-  , COMPASS
    #conf['datasets']['sidis']['xlsx'][3000]='sidis/expdata/3000.xlsx' # deuteron , mult , h+  , COMPASS
    #conf['datasets']['sidis']['xlsx'][3001]='sidis/expdata/3001.xlsx' # deuteron , mult , h-  , COMPASS
    conf['datasets']['sidis']['norm']={}
    data['sidis'] = sidisREAD().load_data_sets('sidis')

    return data

def plot_kin_sidis():

    data = load_data(Q2cut = None, W2SIDIScut = None, zmincut = None, zmaxcut = None)

    data_JAM = load_data(Q2cut = 1.3**2, W2SIDIScut = 20,   zmincut = 0.2, zmaxcut = 0.8)
    data_MAP = load_data(Q2cut = 4,      W2SIDIScut = None, zmincut = 0.2, zmaxcut = 0.8)

    nrows,ncols=2,1
    fig = py.figure(figsize=(ncols*14,nrows*8))
    ax11=py.subplot(nrows,ncols,1)
    ax12=py.subplot(nrows,ncols,2)

    divider = make_axes_locatable(ax11)
    ax11L = divider.append_axes("right",size=6,pad=0,sharey=ax11)
    ax11L.set_xlim(0.1,0.85)
    ax11L.spines['left'].set_visible(False)
    ax11L.yaxis.set_ticks_position('right')
    py.setp(ax11L.get_xticklabels(),visible=True)

    ax11.spines['right'].set_visible(False)

    hand = {}

    exp = 'sidis'
    #--plot x and Q2 (unpolarized)
    hand[exp] = {}
    for idx in data[exp]:
        #--plot all data
        X,Q2 = get_kin_X(exp,data[exp][idx])
        marker,color,label,s = 'o','black','all', 15
        zorder,edgecolors,linewidths = 1.5,'face',1.5
        ax11               .scatter(X,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)
        hand[label] = ax11L.scatter(X,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)

        #--plot z and Q2
        Z,Q2 = get_kin_Z(exp,data[exp][idx])
        label = None
        hand[label] = ax12.scatter(Z,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)

        #--plot data with JAM cuts
        X,Q2 = get_kin_X(exp,data_JAM[exp][idx])
        marker,color,label,s = 'o','red','JAM', 250
        zorder,edgecolors,linewidths = 1.3,'face',1.5
        ax11               .scatter(X,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)
        hand[label] = ax11L.scatter(X,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)

        #--plot z and Q2
        Z,Q2 = get_kin_Z(exp,data_JAM[exp][idx])
        label = None
        hand[label] = ax12.scatter(Z,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)

        #--plot data with MAP cuts
        X,Q2 = get_kin_X(exp,data_MAP[exp][idx])
        marker,color,label,s = 'o','green','MAP', 100
        zorder,edgecolors,linewidths = 1.4,'face',1.5
        ax11               .scatter(X,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)
        hand[label] = ax11L.scatter(X,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)

        #--plot z and Q2
        Z,Q2 = get_kin_Z(exp,data_MAP[exp][idx])
        label = None
        hand[label] = ax12.scatter(Z,Q2,c=color,s=s,marker=marker,zorder=zorder,edgecolors=edgecolors,linewidths=linewidths)

    #--Plot cuts
    x = np.linspace(1e-3,0.999,100)
    #W2cut3  = (3.5 - 0.938**2)*x/(1-x)
    #W2cut10 = (10.0 - 0.938**2)*x/(1-x)

    #hand['W2=3.5']  ,= ax11L.plot(x,W2cut3 ,'k-' ,zorder=2)
    #hand['W2=10']   ,= ax11L.plot(x,W2cut10,'k--',zorder=2)

    ax11 .plot(x,1.3**2*np.ones(len(x)),ls=':',color='red')
    ax11L.plot(x,1.3**2*np.ones(len(x)),ls=':',color='red')

    ax12 .plot(x,1.3**2*np.ones(len(x)),ls=':',color='red')

    ax11 .plot(x,4*     np.ones(len(x)),ls=':',color='green')
    ax11L.plot(x,4*     np.ones(len(x)),ls=':',color='green')

    ax12 .plot(x,4*     np.ones(len(x)),ls=':',color='green')

    #MB = 0.138
    #z = 0.3

    for ax in [ax11, ax11L, ax12]:
        ax.set_yscale('log')
        ax.set_ylim(0.8,39)
 

    ax11 .tick_params(axis='both',which='both',top=True,right=False,direction='in',labelsize=30)
    ax11L.tick_params(axis='both',which='both',top=True,right=True,labelright=False,direction='in',labelsize=30)

    ax11 .tick_params(axis='x',which='major',pad=8)
    ax11L.tick_params(axis='x',which='major',pad=8)

    ax11 .tick_params(axis='both',which='minor',size=4)
    ax11 .tick_params(axis='both',which='major',size=8)
    ax11L.tick_params(axis='both',which='minor',size=4)
    ax11L.tick_params(axis='both',which='major',size=8)

    ax11.set_xscale('log')

    ax11.set_xlim(5e-3,0.1)
    ax11. set_xticks([1e-2])
    ax11. set_xticklabels([r'$10^{-2}$'])
    ax11L.set_xlim(0.1,0.29)
    ax11L.set_xticks([0.1,0.2])
    ax11L.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax11.set_ylabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)
    ax12.set_ylabel(r'\boldmath$Q^2$' + '  ' + r'\textbf{\textrm{(GeV}}' + r'\boldmath$^2)$', size=40)


    ax12 .tick_params(axis='both',which='both',top=True,right=False,direction='in',labelsize=30)

    ax12 .tick_params(axis='x',which='major',pad=8)

    ax12 .tick_params(axis='both',which='minor',size=4)
    ax12 .tick_params(axis='both',which='major',size=8)

    ax12.set_xlim(0.1,1.0)
    ax12.set_xticks([0.2,0.4,0.6,0.8])
    ax12.xaxis.set_minor_locator(MultipleLocator(0.1))

    ax11L.set_xlabel(r'\boldmath$x$',size=50)
    ax12. set_xlabel(r'\boldmath$z$',size=50)
    
    ax11L.xaxis.set_label_coords(0.95,0.00)
    ax12. xaxis.set_label_coords(0.95,0.00)

    ax11L.text(-1.00,0.90,r'\textrm{\textbf{COMPASS kinematics}}',transform = ax11L.transAxes, size=50)

    ax11L.text(0.70,0.13,r'$Q^2 = 1.69$',transform = ax11L.transAxes, size=30,color='red')
    ax11L.text(0.75,0.34,r'$Q^2 = 4$'   ,transform = ax11L.transAxes, size=30,color='green')
    
    ax12 .text(0.90,0.34,r'$Q^2 = 4$'   ,transform = ax12 .transAxes, size=30,color='green')
    ax12 .text(0.85,0.13,r'$Q^2 = 1.69$',transform = ax12 .transAxes, size=30,color='red')

    hand['blank'] ,= ax11.plot(0,0,alpha=0)

    fs = 40

    handles,labels = [], []
    handles.append(hand['all'])
    handles.append(hand['JAM'])
    handles.append(hand['MAP'])
    labels.append(r'\textbf{\textrm{no cuts}}')
    labels.append(r'\textbf{\textrm{JAM}}')
    labels.append(r'\textbf{\textrm{MAP}}')
    ax11.legend(handles,labels,loc=(0.00,0.50),fontsize=fs,frameon=False, handlelength = 0.95, handletextpad = 0.1, ncol = 1, columnspacing = 1.0)

    py.tight_layout()
    py.subplots_adjust(top=0.99,right=0.99)
    filename='gallery/kin_sidis'
    filename+='.png'

    py.savefig(filename)
    print ('Saving figure to %s'%filename)


if __name__ == "__main__":


    plot_kin_sidis()





