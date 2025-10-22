#!/usr/bin/env python
import sys, os
import numpy as np
import copy
import pandas as pd
import scipy as sp
from scipy.interpolate import griddata

# sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))
path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

## matplotlib
import matplotlib
import matplotlib as plt
matplotlib.use('Agg')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
import pylab as py
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
# from analysis.corelib import classifier

import kmeanconf as kc

ext = '.png'
#ext = '.pdf'

legend_fs = 45
text_fs   = 30
label_fs  = 40
xaxis_fs  = 60

cwd = os.getcwd()
checkdir('gallery')

def get_Q2bins(data,kind):
    query = {}
    #--proton A1/Apa
    if kind == 'p A1':
        query[6]  = data.query('Q2 > 60.0')                
        query[5]  = data.query('Q2 > 40.0 and Q2 <= 60.0') 
        query[4]  = data.query('Q2 > 20.0 and Q2 <= 40.0') 
        query[3]  = data.query('Q2 > 10.0 and Q2 <= 20.0') 
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 10.0') 
        query[1]  = data.query('Q2 > 2.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 2.00')
    if kind == 'p A1 SANE':
        query[3]  = data.query('Q2 > 4.00') 
        query[2]  = data.query('Q2 > 3.00 and Q2 <= 4.00') 
        query[1]  = data.query('Q2 > 2.50 and Q2 <= 3.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 2.50') 


    if kind == 'p Apa':
        query[8]  = data.query('Q2 > 20.0') 
        query[7]  = data.query('Q2 > 15.0 and Q2 <= 20.0') 
        query[6]  = data.query('Q2 > 10.0 and Q2 <= 15.0') 
        query[5]  = data.query('Q2 > 7.00 and Q2 <= 10.0')
        query[4]  = data.query('Q2 > 5.00 and Q2 <= 7.00')
        query[3]  = data.query('Q2 > 4.00 and Q2 <= 5.00') 
        query[2]  = data.query('Q2 > 3.00 and Q2 <= 4.00') 
        query[1]  = data.query('Q2 > 2.00 and Q2 <= 3.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 2.00') 
    if kind == 'p Apa DVCS':
        query[6]  = data.query('Q2 > 4.20') 
        query[5]  = data.query('Q2 > 3.50 and Q2 <= 4.20') 
        query[4]  = data.query('Q2 > 2.90 and Q2 <= 3.50') 
        query[3]  = data.query('Q2 > 2.50 and Q2 <= 2.90') 
        query[2]  = data.query('Q2 > 2.40 and Q2 <= 2.50') 
        query[1]  = data.query('Q2 > 2.00 and Q2 <= 2.40') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 2.00')
    if kind == 'p Apa eg1b E4.2':
        query[3]  = data.query('Q2 > 2.40') 
        query[2]  = data.query('Q2 > 2.10 and Q2 <= 2.40') 
        query[1]  = data.query('Q2 > 1.80 and Q2 <= 2.10')
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 1.80')
    if kind == 'p Apa eg1b E5.7':
        query[6]  = data.query('Q2 > 4.10') 
        query[5]  = data.query('Q2 > 3.40 and Q2 <= 4.10') 
        query[4]  = data.query('Q2 > 2.90 and Q2 <= 3.40') 
        query[3]  = data.query('Q2 > 2.60 and Q2 <= 2.90') 
        query[2]  = data.query('Q2 > 2.20 and Q2 <= 2.60') 
        query[1]  = data.query('Q2 > 1.80 and Q2 <= 2.20')
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 1.80')


    #--proton A2/Ape/Atpe
    if kind == 'p A2':
        query[3]  = data.query('Q2 > 5.00') 
        query[2]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[1]  = data.query('Q2 > 2.20 and Q2 <= 3.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 2.20')


    if kind == 'p Ape':
        query[4]  = data.query('Q2 > 15.0') 
        # iter = query[4][((query[4].X == 0.750) &( query[4].Q2 == 23.72))].index
        # query[4] = query[4].drop(iter)
        iter = query[4][((query[4].X == 0.845) &( query[4].Q2 == 24.84))].index
        query[4] = query[4].drop(iter)
        query[3]  = data.query('Q2 > 7.00 and Q2 <= 15.0')
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 7.00')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 
    if kind == 'p Atpe':
        query[4]  = data.query('Q2 > 15.0') 
        iter = query[4][((query[4].X == 0.844) &( query[4].Q2 == 20.45))].index
        query[4] = query[4].drop(iter)
        iter = query[4][((query[4].X == 0.845) &( query[4].Q2 == 17.10))].index
        query[4] = query[4].drop(iter)
        query[3]  = data.query('Q2 > 7.00 and Q2 <= 15.0')
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 7.00')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 





    #--deuteron A1/Apa
    if kind == 'd A1':
        query[4]  = data.query('Q2 > 40.0') 
        query[3]  = data.query('Q2 > 20.0 and Q2 <= 40.0') 
        query[2]  = data.query('Q2 > 10.0 and Q2 <= 20.0') 
        query[1]  = data.query('Q2 > 5.00 and Q2 <= 10.0') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 5.00')


    if kind == 'd Apa':
        query[8]  = data.query('Q2 > 20.0') 
        query[7]  = data.query('Q2 > 15.0 and Q2 <= 20.0') 
        query[6]  = data.query('Q2 > 10.0 and Q2 <= 15.0') 
        query[5]  = data.query('Q2 > 7.00 and Q2 <= 10.0')
        query[4]  = data.query('Q2 > 5.00 and Q2 <= 7.00')
        query[3]  = data.query('Q2 > 4.00 and Q2 <= 5.00') 
        query[2]  = data.query('Q2 > 3.00 and Q2 <= 4.00') 
        query[1]  = data.query('Q2 > 2.00 and Q2 <= 3.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 2.00') 
    if kind == 'd Apa JLab':
        query[10] = data.query('Q2 > 4.50') 
        query[9]  = data.query('Q2 > 4.00 and Q2 <= 4.50') 
        query[8]  = data.query('Q2 > 3.50 and Q2 <= 4.00') 
        query[7]  = data.query('Q2 > 3.00 and Q2 <= 3.50') 
        query[6]  = data.query('Q2 > 2.50 and Q2 <= 3.00') 
        query[5]  = data.query('Q2 > 2.40 and Q2 <= 2.50') 
        query[4]  = data.query('Q2 > 2.04 and Q2 <= 2.40') 
        query[3]  = data.query('Q2 > 2.00 and Q2 <= 2.04') 
        query[2]  = data.query('Q2 > 1.71 and Q2 <= 2.00') 
        query[1]  = data.query('Q2 > 1.70 and Q2 <= 1.71')
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 1.70')


    #--deuteron Ape/Atpe
    if kind == 'd Ape':
        query[4]  = data.query('Q2 > 15.0') 
        iter = query[4][((query[4].X == 0.845) &( query[4].Q2 == 24.94))].index
        query[4] = query[4].drop(iter)
        query[3]  = data.query('Q2 > 7.00 and Q2 <= 15.0')
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 7.00')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 
    if kind == 'd Atpe':
        query[4]  = data.query('Q2 > 15.0') 
        iter = query[4][((query[4].X == 0.844) &( query[4].Q2 == 20.45))].index
        query[4] = query[4].drop(iter)
        iter = query[4][((query[4].X == 0.845) &( query[4].Q2 == 17.10))].index
        query[4] = query[4].drop(iter)
        query[3]  = data.query('Q2 > 7.00 and Q2 <= 15.0')
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 7.00')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 

    #--helium 
    if kind == 'h A1':
        query[0]  = data.query('Q2 > 1.60') 
    if kind == 'h Apa':
        query[3]  = data.query('Q2 > 10.00') 
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 10.0')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 
    if kind == 'h A2':
        query[0]  = data.query('Q2 > 1.60') 
    if kind == 'h Ape':
        query[3]  = data.query('Q2 > 10.0') 
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 10.0')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 

    if kind == 'JLabHAE06014':
        query[3]  = data.query('Q2 > 10.0') 
        query[2]  = data.query('Q2 > 5.00 and Q2 <= 10.0')
        query[1]  = data.query('Q2 > 3.00 and Q2 <= 5.00') 
        query[0]  = data.query('Q2 > 1.60 and Q2 <= 3.00') 
 
    return query

def get_plot(query,cluster_i=0):

    #--generate dictionary with everything needed for plot

    plot = {_:{} for _ in ['theory','X','Q2','value','alpha','std']}
    for key in query:
        theory = query[key]['thy-%d' % cluster_i]
        std    = query[key]['dthy-%d' % cluster_i]
        X      = query[key]['X']
        Q2     = query[key]['Q2']
        value  = query[key]['value']
        alpha  = query[key]['alpha']
        #--sort by ascending Q2
        zx = sorted(zip(Q2,X))
        zt = sorted(zip(Q2,theory))
        zv = sorted(zip(Q2,value))
        za = sorted(zip(Q2,alpha))
        zs = sorted(zip(Q2,std))
        plot['X'][key]      = np.array([zx[i][1] for i in range(len(zx))])
        plot['theory'][key] = np.array([zt[i][1] for i in range(len(zt))])
        plot['value'][key]  = np.array([zv[i][1] for i in range(len(zv))])
        plot['alpha'][key]  = np.array([za[i][1] for i in range(len(za))])
        plot['std'][key]    = np.array([zs[i][1] for i in range(len(zs))])
        plot['Q2'][key]     = np.array(sorted(Q2))

    return plot

def get_theory(PLOT,nbins,loop=True,funcQ2=True,funcX=False):

    #--interpolate theoretical values across Q2 or X

    theory = {}

    if funcQ2: svar = 'Q2'
    if funcX:  svar = 'X'

    theory = {_:{} for _ in [svar,'value','std']}
    for key in range(nbins):
        var = []
        thy = []
        std = []
        
        #--if plotting for multiple experiments, loop over and combine
        if loop:
            for exp in PLOT:
                var. extend(PLOT[exp][svar][key])
                thy.extend(PLOT[exp]['theory'][key])
                std.extend(PLOT[exp]['std'][key])
        else:
            var.extend(PLOT[svar][key])
            thy.extend(PLOT['theory'][key])
            std.extend(PLOT['std'][key])

        #--if nothing in bin, skip
        if len(var) == 0: continue

        vmin = np.min(var)
        vmax = np.max(var)
        theory[svar][key]  = np.geomspace(vmin,vmax,100)

        #--if more than one value, interpolate between them
        if len(var) > 1:
            theory['value'][key] = griddata(np.array(var),np.array(thy),theory[svar][key],method='linear')
            theory['std'][key]   = griddata(np.array(var),np.array(std),theory[svar][key],method='linear')
        else:
            theory['value'][key] = np.ones(100)*thy 
            theory['std'][key]   = np.ones(100)*std


    return theory

def get_details(exp):

    #--get details for plotting
    combo1 = ('firebrick', '*', 10)
    combo2 = ('darkgreen', '^', 8)
    combo3 = ('blue'     , 'o', 8)
    combo4 = ('purple'   , 'v', 8)
    combo5 = ('black'    , 's', 7)
    combo6 = ('orange'   , 'D', 7)

    if exp=='COMPASS':     color, marker, ms = combo1
    if exp=='COMPASS1':    color, marker, ms = combo1
    if exp=='COMPASS2':    color, marker, ms = combo1
    if exp=='EMC' :        color, marker, ms = combo2
    if exp=='SMC1':        color, marker, ms = combo3
    if exp=='SMC2':        color, marker, ms = combo3
    if exp=='SANE1':       color, marker, ms = combo4
    if exp=='SANE2':       color, marker, ms = combo5

    if exp=='HERMES':            color, marker, ms = combo1
    if exp=='SLACE143':          color, marker, ms = combo2
    if exp=='SLACE155':          color, marker, ms = combo3
    if exp=='SLACE80E130':       color, marker, ms = combo4
    if exp=='JLabHBEG1DVCSE4.8': color, marker, ms = combo1
    if exp=='JLabHBEG1DVCSE6.0': color, marker, ms = combo2
    if exp=='JLabHBEG1b4.2':     color, marker, ms = combo3
    if exp=='JLabHBEG1b5.7':     color, marker, ms = combo4
    if exp=='JLabHBEG1DVCS': color, marker, ms = combo2


    if exp=='SLACE142':      color, marker, ms = combo2
    if exp=='SLACE154':      color, marker, ms = combo2
    if exp=='JLabHAE06014':  color, marker, ms = combo3
    if exp=='JLabHAE99117':  color, marker, ms = combo4
    
    if exp=='SLACE155x_E29_t2.75': color, marker, ms = combo1
    if exp=='SLACE155x_E29_t5.50': color, marker, ms = combo2
    if exp=='SLACE155x_E29_t10.5': color, marker, ms = combo3
    if exp=='SLACE155x_E32_t2.75': color, marker, ms = combo4
    if exp=='SLACE155x_E32_t5.50': color, marker, ms = combo5
    if exp=='SLACE155x_E32_t10.5': color, marker, ms = combo6

    if exp=='JLab E12 Apa':     color, marker, ms = combo1
    if exp=='JLab E12 Ape':     color, marker, ms = combo2

    return color,marker,ms


#--proton plots
def plot_A1_proton(data):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #######################################
    #--Plot A1 proton for COMPASS, EMC, SMC 
    #######################################

    ax = ax11
    DATA = {}
    DATA['COMPASS1']  = pd.DataFrame(data[10002]) #--COMPASS A1
    DATA['COMPASS2']  = pd.DataFrame(data[10003]) #--COMPASS A1
    DATA['EMC']       = pd.DataFrame(data[10004]) #--EMC A1
    DATA['SMC1']      = pd.DataFrame(data[10035]) #--SMC A1
    DATA['SMC2']      = pd.DataFrame(data[10036]) #--SMC A1

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p A1')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot, thy_band = {}, {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='EMC': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy   = theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            color,marker,ms = get_details(exp)
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)

            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)
            

    #--plot labels
    ax.text(0.49, 0.89,  r'$Q^2 > 60~{\rm GeV^2} \, (i=6)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.49, 0.74,  r'$40 < Q^2 < 60$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.59,  r'$20 < Q^2 < 40$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.24, 0.44,  r'$10 < Q^2 < 20$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.68, 0.32,  r'$5 < Q^2 < 10$'                   , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.57, 0.17,  r'$2 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.27, 0.03,  r'$m_c^2 < Q^2 < 2\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.16, 0.01
    ax.text(x,   y    ,r'\boldmath$A_1$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',             transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.005, 0.7)
    ax.set_ylim(-0.2, 7.2)
    ax.set_xticks([0.01,0.1,0.5])
    ax.set_xticklabels([r'$0.01$',r'$0.1$',r'$0.5$'])

    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    legend_fs = 45
    handles,labels = [],[]
    handles.append((thy_band['COMPASS1'],thy_plot['COMPASS1'],hand['COMPASS1']))
    handles.append((thy_band['EMC'],thy_plot['EMC'],hand['EMC']))
    handles.append((thy_band['SMC1'],thy_plot['SMC1'],hand['SMC1']))
    labels.append(r'\textbf{\textrm{COMPASS}}')
    labels.append(r'\textbf{\textrm{EMC}}')
    labels.append(r'\textbf{\textrm{SMC}}')
    ax.legend(handles,labels,loc=(0.01,0.60), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    #######################################
    #--Plot A1 proton for SANE 
    #######################################

    ax = ax12
    DATA = {}
    DATA['SANE1']     = pd.DataFrame(data[20001]) #--SANE A1
    DATA['SANE2']     = pd.DataFrame(data[20002]) #--SANE A1

    N = 1.0
    hand = {}
    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p A1 SANE')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot, thy_band = {}, {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SANE1': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy   = theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            color,marker,ms = get_details(exp)
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)

            if len(X)==1:
                x = [X_thy[0]-0.005,X_thy[0]+0.005]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)

    #--plot labels
    ax.text(0.44, 0.86,  r'$Q^2 > 4~{\rm GeV^2} \, (i=3)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.20, 0.60,  r'$3 < Q^2 < 4$'                   , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.70, 0.40,  r'$2.5 < Q^2 < 3$'                 , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.60, 0.17,  r'$m_c^2 < Q^2 < 2.5\, (i=0)$'     , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)
    #ax.set_xticks([0.01,0.1])
    #ax.set_xticklabels([r'$0.01$',r'$0.1$'])


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.22, 0.58)
    ax.set_ylim(-0.2, 4.2)

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))

    legend_fs = 45
    ax.text(0.02, 0.95,  r'\textrm{\textbf{SANE}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SANE1'],thy_plot['SANE1'],hand['SANE1']))
    handles.append((thy_band['SANE2'],thy_plot['SANE2'],hand['SANE2']))
    labels.append(r'\textbf{\textrm{E4.7}}')
    labels.append(r'\textbf{\textrm{E5.9}}')
    ax.legend(handles,labels,loc=(0.01,0.78), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.title(r'\rm \bf proton target',fontsize = 60,y = 1.02, x = 0)

    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_A1_proton'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_Apa_proton(data):


    nrows, ncols = 2, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)
    ax21 = py.subplot(nrows, ncols, 3)
    ax22 = py.subplot(nrows, ncols, 4)

    #######################################
    #--Plot Apa proton for HERMES and SLAC 
    #######################################

    ax = ax11
    DATA = {}
    DATA['HERMES']        = pd.DataFrame(data[10007])  #--HERMES Apa
    DATA['SLACE143']      = pd.DataFrame(data[10022])  #--SLAC E143 Apa
    DATA['SLACE155']      = pd.DataFrame(data[10029])  #--SLAC E155 Apa
    DATA['SLACE80E130']   = pd.DataFrame(data[10032])  #--SLAC E80E130 Apa

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p Apa')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='HERMES': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.44, 0.91,  r'$Q^2 > 20~{\rm GeV^2} \, (i=8)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.49, 0.79,  r'$15 < Q^2 < 20$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.40, 0.69,  r'$10 < Q^2 < 15$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.57,  r'$7 < Q^2 < 10$'                   , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.35, 0.45,  r'$5 < Q^2 < 7$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.27, 0.34,  r'$4 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.32, 0.24,  r'$3 < Q^2 < 4$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.18, 0.13,  r'$2 < Q^2 < 3$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.02, 0.04,  r'$m_c^2 < Q^2 < 2\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)
    ax.set_xticks([0.01,0.1])
    ax.set_xticklabels([r'$0.01$',r'$0.1$'])

    x, y   = 0.03, 0.89
    dx, dy = 0.16, 0.01
    ax.text(x,   y    ,r'\boldmath$A_{\parallel}$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                       transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.005, 0.9)
    ax.set_ylim(-0.2, 9.2)
    ax.set_xticks([0.01,0.1])
    ax.set_xticklabels([r'$0.01$',r'$0.1$'])

    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    legend_fs = 45
    handles,labels = [],[]
    handles.append((thy_band['HERMES'],thy_plot['HERMES'],hand['HERMES']))
    labels.append(r'\textbf{\textrm{HERMES}}')

    legend = ax.legend(handles,labels,loc=(0.01,0.75), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)
    ax.add_artist(legend)

    ax.text(0.02, 0.70,  r'\textrm{\textbf{SLAC}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SLACE143'],thy_plot['SLACE143'],hand['SLACE143']))
    handles.append((thy_band['SLACE155'],thy_plot['SLACE155'],hand['SLACE155']))
    handles.append((thy_band['SLACE80E130'],thy_plot['SLACE80E130'],hand['SLACE80E130']))
    labels.append(r'\textbf{\textrm{E143}}')
    labels.append(r'\textbf{\textrm{E155}}')
    labels.append(r'\textbf{\textrm{E80/E130}}')
    ax.legend(handles,labels,loc=(0.01,0.45), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #######################################
    #--Plot Apa proton for JLab 
    #######################################

    ax = ax12
    DATA = {}
    DATA['JLabHBEG1DVCSE4.8'] = pd.DataFrame(data[10017]).query('Elab==4.8')  #--JLabHBEG1DVCS
    DATA['JLabHBEG1DVCSE6.0'] = pd.DataFrame(data[10017]).query('Elab==6.0')  #--JLabHBEG1DVCS

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p Apa DVCS')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    for exp in PLOT:
        for key in range(nbins):
            if exp=='JLabHBEG1DVCSE4.8': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.40, 0.90,  r'$Q^2 > 4.2~{\rm GeV^2} \, (i=6)$', transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.33, 0.76,  r'$3.5 < Q^2 < 4.2$'                , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.17, 0.62,  r'$2.9 < Q^2 < 3.5$'                , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.67, 0.48,  r'$2.5 < Q^2 < 2.9$'                , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.58, 0.34,  r'$2.4 < Q^2 < 2.5$'                , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.50, 0.20,  r'$2.0 < Q^2 < 2.4$'                , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.06,  r'$m_c^2 < Q^2 < 2.0\, (i=0)$'      , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.15, 0.65)
    ax.set_ylim(-0.2, 7.2)
    ax.set_xticks([0.2,0.3,0.4,0.5,0.6])

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))

    ax.text(0.02, 0.95,  r'\textrm{\textbf{JLab Hall B eg1 DVCS}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['JLabHBEG1DVCSE4.8'],thy_plot['JLabHBEG1DVCSE4.8'],hand['JLabHBEG1DVCSE4.8']))
    handles.append((thy_band['JLabHBEG1DVCSE6.0'],thy_plot['JLabHBEG1DVCSE6.0'],hand['JLabHBEG1DVCSE6.0']))
    labels.append(r'\textbf{\textrm{E4.8}}')
    labels.append(r'\textbf{\textrm{E6.0}}')
    ax.legend(handles,labels,loc=(0.01,0.78), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    ax.set_title(r'\rm \bf proton target',fontsize = 60,y = 1.02, x = -0.06)

    #######################################
    #--Plot Apa proton for JLab 
    #######################################
    ax = ax21
    DATA = {}
    DATA['JLabHBEG1b4.2'] = pd.DataFrame(data[10043])  #--JLabHBEG1b, Elab = 4.2

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p Apa eg1b E4.2')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    for exp in PLOT:
        for key in range(nbins):
            if exp=='JLabHBEG1b4.2': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.65, 0.57,  r'$Q^2 > 2.4~{\rm GeV^2} \, (i=3)$', transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.15, 0.40,  r'$2.1 < Q^2 < 2.4$'               , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.06, 0.24,  r'$1.8 < Q^2 < 2.1$'               , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.43, 0.07,  r'$m_c^2 < Q^2 < 1.8\, (i=0)$'     , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.15, 0.65)
    ax.set_ylim(-0.2, 6.2)
    ax.set_xticks([0.2,0.3,0.4,0.5,0.6])

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))

    ax.text(0.02, 0.95,  r'\textrm{\textbf{JLab Hall B}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['JLabHBEG1b4.2'],thy_plot['JLabHBEG1b4.2'],hand['JLabHBEG1b4.2']))
    labels.append(r'\textbf{\textrm{eg1b (E4.2)}}')
    ax.legend(handles,labels,loc=(0.01,0.84), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)



    #######################################
    #--Plot Apa proton for JLab 
    #######################################
    ax = ax22
    DATA = {}
    DATA['JLabHBEG1b5.7'] = pd.DataFrame(data[10044])  #--JLabHBEG1b, Elab = 5.7

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p Apa eg1b E5.7')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    for exp in PLOT:
        for key in range(nbins):
            if exp=='JLabHBEG1b5.7': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    t1=ax.text(0.45, 0.80,  r'$Q^2 > 4.1~{\rm GeV^2} \, (i=6)$' , transform = ax.transAxes, fontsize = text_fs)
    t2=ax.text(0.29, 0.66,  r'$3.4 < Q^2 < 4.1$'                , transform = ax.transAxes, fontsize = text_fs)
    t3=ax.text(0.20, 0.55,  r'$2.9 < Q^2 < 3.4$'                , transform = ax.transAxes, fontsize = text_fs)
    t4=ax.text(0.68, 0.42,  r'$2.6 < Q^2 < 2.9$'                , transform = ax.transAxes, fontsize = text_fs)
    t5=ax.text(0.59, 0.29,  r'$2.2 < Q^2 < 2.6$'                , transform = ax.transAxes, fontsize = text_fs)
    t6=ax.text(0.50, 0.17,  r'$1.8 < Q^2 < 2.2$'                , transform = ax.transAxes, fontsize = text_fs)
    t7=ax.text(0.42, 0.04,  r'$m_c^2 < Q^2 < 1.8\, (i=0)$'      , transform = ax.transAxes, fontsize = text_fs)

    t1.set_bbox(dict(facecolor='white',alpha=0.8,edgecolor='white'))
    t4.set_bbox(dict(facecolor='white',alpha=0.8,edgecolor='white'))

    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.15, 0.65)
    ax.set_ylim(-0.2, 8.2)
    ax.set_xticks([0.2,0.3,0.4,0.5,0.6])

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))

    ax.text(0.02, 0.95,  r'\textrm{\textbf{JLab Hall B}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['JLabHBEG1b5.7'],thy_plot['JLabHBEG1b5.7'],hand['JLabHBEG1b5.7']))
    labels.append(r'\textbf{\textrm{eg1b (E5.7)}}')
    ax.legend(handles,labels,loc=(0.01,0.84), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_Apa_proton'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    ax21.set_rasterized(True)
    ax22.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_A2_proton(data):


    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)

    #######################################
    #--Plot A2 proton for HERMES and SANE
    #######################################

    ax = ax11
    DATA = {}
    DATA['HERMES']  = pd.DataFrame(data[10008]) #--HERMES A2
    DATA['SANE1']   = pd.DataFrame(data[20003]) #--SANE A2
    DATA['SANE2']   = pd.DataFrame(data[20004]) #--SANE A2

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p A2')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot, thy_band = {}, {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='HERMES': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]-0.005,X_thy[0]+0.005]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)

    #--plot labels
    t1=ax.text(0.02, 0.76,  r'$Q^2 > 5~{\rm GeV^2} \, (i=3)$'  , transform = ax.transAxes, fontsize = text_fs)
    t2=ax.text(0.02, 0.52,  r'$3 < Q^2 < 5$'                   , transform = ax.transAxes, fontsize = text_fs)
    t3=ax.text(0.66, 0.28,  r'$2.2 < Q^2 < 3$'                 , transform = ax.transAxes, fontsize = text_fs)
    t4=ax.text(0.50, 0.06,  r'$m_c^2 < Q^2 < 2.2\, (i=0)$'     , transform = ax.transAxes, fontsize = text_fs)

    t4.set_bbox(dict(facecolor='white',alpha=0.8,edgecolor='white'))

    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.17, 0.01
    ax.text(x,   y    ,r'\boldmath$A_2$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',             transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    #ax.semilogx()
    ax.set_xlim(0.03, 0.7)
    ax.set_ylim(-0.2, 4.2)
    ax.set_xticks([0.2,0.4,0.6])
    #ax.set_xticklabels([r'$0.01$',r'$0.1$'])
    minorLocator = MultipleLocator(0.05)
    ax.xaxis.set_minor_locator(minorLocator)
    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    #handles,labels = [],[]
    #handles.append((thy_band['HERMES'],thy_plot['HERMES'],hand['HERMES']))
    #handles.append((thy_band['SANE1'],thy_plot['SANE1'],hand['SANE1']))
    #labels.append(r'\textbf{\textrm{HERMES}}')
    #labels.append(r'\textbf{\textrm{SANE}}')
    #ax.legend(handles,labels,loc=(0.64,0.85), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    ax.text(0.42, 0.95,  r'\textrm{\textbf{SANE}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SANE1'],thy_plot['SANE1'],hand['SANE1']))
    handles.append((thy_band['SANE2'],thy_plot['SANE2'],hand['SANE2']))
    handles.append((thy_band['HERMES'],thy_plot['HERMES'],hand['HERMES']))
    labels.append(r'\textbf{\textrm{E4.7}}')
    labels.append(r'\textbf{\textrm{E5.9}}')
    labels.append(r'\textbf{\textrm{HERMES}}')
    ax.legend(handles,labels,loc=(0.40,0.78), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0, ncol = 2,columnspacing=0.8)


    py.title(r'\rm \bf proton target',fontsize = 60)

    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_A2_proton'%cwd
    filename += ext
    ax11.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_Ape_proton(data):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #######################################
    #--Plot Ape proton for SLAC
    #######################################

    ax = ax11
    DATA = {}
    DATA['SLACE143']      = pd.DataFrame(data[10023])  #--SLAC E143 Ape
    DATA['SLACE155']      = pd.DataFrame(data[10028])  #--SLAC E155 Ape

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p Ape')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE155': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.33, 0.71,  r'$Q^2 > 15~{\rm GeV^2} \, (i=4)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.16, 0.54,  r'$7 < Q^2 < 15$'                   , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.14, 0.37,  r'$5 < Q^2 < 7$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.02, 0.16,  r'$3 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.59, 0.03,  r'$m_c^2 < Q^2 < 3\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)
    ax.set_xticks([0.01,0.1])
    ax.set_xticklabels([r'$0.01$',r'$0.1$'])

    x, y   = 0.03, 0.89
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$A_{\perp}$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                       transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.04, 0.9)
    ax.set_ylim(-0.2, 5.6)
    ax.set_xticks([0.05,0.1,0.3,0.7])
    ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.3$',r'$0.7$'])

    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    ax.text(0.02, 0.80,  r'\textrm{\textbf{SLAC}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SLACE143'],thy_plot['SLACE143'],hand['SLACE143']))
    handles.append((thy_band['SLACE155'],thy_plot['SLACE155'],hand['SLACE155']))
    labels.append(r'\textbf{\textrm{E143}}')
    labels.append(r'\textbf{\textrm{E155}}')
    ax.legend(handles,labels,loc=(0.01,0.62), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #######################################
    #--Plot Atpe proton for SLAC
    #######################################

    ax = ax12
    DATA = {}
    DATA['SLACE155x_E29_t2.75'] = pd.DataFrame(data[10031]).query('Elab==29.1 and theta==2.75')
    DATA['SLACE155x_E29_t5.50'] = pd.DataFrame(data[10031]).query('Elab==29.1 and theta==5.5')
    DATA['SLACE155x_E29_t10.5'] = pd.DataFrame(data[10031]).query('Elab==29.1 and theta==10.5')
    DATA['SLACE155x_E32_t2.75'] = pd.DataFrame(data[10031]).query('Elab==32.3 and theta==2.75')
    DATA['SLACE155x_E32_t5.50'] = pd.DataFrame(data[10031]).query('Elab==32.3 and theta==5.5')
    DATA['SLACE155x_E32_t10.5'] = pd.DataFrame(data[10031]).query('Elab==32.3 and theta==10.5')

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'p Atpe')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE155': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.44, 0.72,  r'$Q^2 > 15~{\rm GeV^2} \, (i=4)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.23, 0.54,  r'$7 < Q^2 < 15$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.18, 0.37,  r'$5 < Q^2 < 7$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.03, 0.19,  r'$3 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.68, 0.05,  r'$m_c^2 < Q^2 < 3\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)

    x, y   = 0.03, 0.86
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$\tilde{A}_{\perp}$' , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                       transform = ax.transAxes, size = 60)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    ax.semilogx()
    ax.set_xlim(0.04, 0.9)
    ax.set_ylim(-0.2, 5.6)
    ax.set_xticks([0.05,0.1,0.3,0.7])
    ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.3$',r'$0.7$'])

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))
    ax.text(0.42, 0.95,  r'\textrm{\textbf{SLAC E29}}'      , transform = ax.transAxes, fontsize = legend_fs)
    ax.text(0.77, 0.95,  r'\textrm{\textbf{E32}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SLACE155x_E29_t2.75'],thy_plot['SLACE155x_E29_t2.75'],hand['SLACE155x_E29_t2.75']))
    handles.append((thy_band['SLACE155x_E29_t5.50'],thy_plot['SLACE155x_E29_t5.50'],hand['SLACE155x_E29_t5.50']))
    handles.append((thy_band['SLACE155x_E29_t10.5'],thy_plot['SLACE155x_E29_t10.5'],hand['SLACE155x_E29_t10.5']))
    handles.append((thy_band['SLACE155x_E32_t2.75'],thy_plot['SLACE155x_E32_t2.75'],hand['SLACE155x_E32_t2.75']))
    handles.append((thy_band['SLACE155x_E32_t5.50'],thy_plot['SLACE155x_E32_t5.50'],hand['SLACE155x_E32_t5.50']))
    handles.append((thy_band['SLACE155x_E32_t10.5'],thy_plot['SLACE155x_E32_t10.5'],hand['SLACE155x_E32_t10.5']))
    labels.append(r'\boldmath$\theta=2.75$')
    labels.append(r'\boldmath$\theta=5.50$')
    labels.append(r'\boldmath$\theta=10.5$')
    labels.append(r'\boldmath$\theta=2.75$')
    labels.append(r'\boldmath$\theta=5.50$')
    labels.append(r'\boldmath$\theta=10.5$')
    ax.legend(handles,labels,loc=(0.40,0.78), fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0, ncol = 2, columnspacing = 4.5)

    py.title(r'\rm \bf proton target',fontsize = 60,y = 1.02, x = 0)


    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_Ape_proton'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()


#--deuteron plots
def plot_A1_deuteron(data):


    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)

    #######################################
    #--Plot A1 deuteron for COMPASS and SMC 
    #######################################

    ax = ax11
    DATA = {}
    DATA['COMPASS']   = pd.DataFrame(data[10001]) #--COMPASS A1
    DATA['SMC1']      = pd.DataFrame(data[10033]) #--SMC A1
    DATA['SMC2']      = pd.DataFrame(data[10034]) #--SMC A1

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'d A1')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot, thy_band = {},{}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='COMPASS': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)

    #--plot labels
    ax.text(0.51, 0.81,  r'$Q^2 > 40~{\rm GeV^2} \, (i=4)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.59,  r'$20 < Q^2 < 40$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.32, 0.41,  r'$10 < Q^2 < 20$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.58, 0.23,  r'$5 < Q^2 < 10$'                   , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.43, 0.03,  r'$m_c^2 < Q^2 < 5\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$A_1$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',             transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.005, 0.7)
    ax.set_ylim(-0.2, 5.2)
    ax.set_xticks([0.01,0.1,0.5])
    ax.set_xticklabels([r'$0.01$',r'$0.1$',r'$0.5$'])

    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    handles,labels = [],[]
    handles.append((thy_band['COMPASS'],thy_plot['COMPASS'],hand['COMPASS']))
    handles.append((thy_band['SMC1'],thy_plot['SMC1'],hand['SMC1']))
    labels.append(r'\textbf{\textrm{COMPASS}}')
    labels.append(r'\textbf{\textrm{SMC}}')
    ax.legend(handles,labels,loc=(0.01,0.70), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.title(r'\rm \bf deuteron target',fontsize = 60)
    
    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_A1_deuteron'%cwd
    filename += ext
    ax11.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_Apa_deuteron(data):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #######################################
    #--Plot Apa proton for HERMES and SLAC 
    #######################################

    ax = ax11
    DATA = {}
    DATA['HERMES']        = pd.DataFrame(data[10006])  #--HERMES Apa
    DATA['SLACE143']      = pd.DataFrame(data[10021])  #--SLAC E143 Apa
    DATA['SLACE155']      = pd.DataFrame(data[10027])  #--SLAC E155 Apa

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'d Apa')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='HERMES': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.44, 0.90,  r'$Q^2 > 20~{\rm GeV^2} \, (i=8)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.49, 0.78,  r'$15 < Q^2 < 20$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.40, 0.67,  r'$10 < Q^2 < 15$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.55,  r'$7 < Q^2 < 10$'                   , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.35, 0.45,  r'$5 < Q^2 < 7$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.28, 0.34,  r'$4 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.27, 0.24,  r'$3 < Q^2 < 4$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.18, 0.13,  r'$2 < Q^2 < 3$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.02, 0.04,  r'$m_c^2 < Q^2 < 2\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$A_{\parallel}$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                       transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.005, 0.9)
    ax.set_ylim(-0.2, 9.2)
    ax.set_xticks([0.01,0.1,0.5])
    ax.set_xticklabels([r'$0.01$',r'$0.1$',r'$0.5$'])

    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    handles,labels = [],[]
    handles.append((thy_band['HERMES'],thy_plot['HERMES'],hand['HERMES']))
    labels.append(r'\textbf{\textrm{HERMES}}')

    legend = ax.legend(handles,labels,loc=(0.01,0.75), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)
    ax.add_artist(legend)

    ax.text(0.02, 0.72,  r'\textrm{\textbf{SLAC}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SLACE143'],thy_plot['SLACE143'],hand['SLACE143']))
    handles.append((thy_band['SLACE155'],thy_plot['SLACE155'],hand['SLACE155']))
    labels.append(r'\textbf{\textrm{E143}}')
    labels.append(r'\textbf{\textrm{E155}}')
    ax.legend(handles,labels,loc=(0.01,0.54), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #######################################
    #--Plot Apa proton for JLab 
    #######################################

    ax = ax12
    DATA = {}
    DATA['JLabHBEG1DVCS'] = pd.DataFrame(data[10016])  #--JLabHBEG1DVCS
    DATA['JLabHBEG1b4.2'] = pd.DataFrame(data[10039])  #--JLabHBEG1b, Elab = 4.2
    DATA['JLabHBEG1b5.7'] = pd.DataFrame(data[10040])  #--JLabHBEG1b, Elab = 5.7

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'d Apa JLab')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    for exp in PLOT:
        for key in range(nbins):
            if exp=='JLabHBEG1DVCS': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X    = theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.46, 0.83,  r'$Q^2 > 4.5~{\rm GeV^2} \, (i=10)$', transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.48, 0.74,  r'$4 < Q^2 < 4.5$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.38, 0.66,  r'$3.5 < Q^2 < 4$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.24, 0.57,  r'$3 < Q^2 < 3.5$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.68, 0.49,  r'$2.5 < Q^2 < 3$'                  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.60, 0.40,  r'$2.4 < Q^2 < 2.5$'                , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.53, 0.32,  r'$2.04 < Q^2 < 2.4$'               , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.50, 0.25,  r'$2 < Q^2 < 2.04$'                 , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.47, 0.18,  r'$1.71 < Q^2 < 2$'                 , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.10,  r'$1.7 < Q^2 < 1.71$'               , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.42, 0.02,  r'$m_c^2 < Q^2 < 1.7\, (i=0)$'      , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.15, 0.65)
    ax.set_ylim(-0.2, 12.5)
    ax.set_xticks([0.2,0.3,0.4,0.5,0.6])

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))

    ax.text(0.02, 0.95,  r'\textrm{\textbf{JLab Hall B}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['JLabHBEG1DVCS'],thy_plot['JLabHBEG1DVCS'],hand['JLabHBEG1DVCS']))
    handles.append((thy_band['JLabHBEG1b4.2'],thy_plot['JLabHBEG1b4.2'],hand['JLabHBEG1b4.2']))
    handles.append((thy_band['JLabHBEG1b5.7'],thy_plot['JLabHBEG1b5.7'],hand['JLabHBEG1b5.7']))
    labels.append(r'\textbf{\textrm{eg1 DVCS}}')
    labels.append(r'\textbf{\textrm{eg1b (E4.2)}}')
    labels.append(r'\textbf{\textrm{eg1b (E5.7)}}')
    ax.legend(handles,labels,loc=(0.01,0.68), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.title(r'\rm \bf deuteron target',fontsize = 60,y = 1.02, x = 0)

    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_Apa_deuteron'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_Ape_deuteron(data):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #######################################
    #--Plot Ape deuteron for SLAC
    #######################################

    ax = ax11
    DATA = {}
    DATA['SLACE143']      = pd.DataFrame(data[10020])  #--SLAC E143 Ape
    DATA['SLACE155']      = pd.DataFrame(data[10026])  #--SLAC E155 Ape

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'d Ape')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE155': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.33, 0.71,  r'$Q^2 > 15~{\rm GeV^2} \, (i=4)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.16, 0.54,  r'$7 < Q^2 < 15$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.14, 0.37,  r'$5 < Q^2 < 7$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.02, 0.22,  r'$3 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.59, 0.03,  r'$m_c^2 < Q^2 < 3\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$A_{\perp}$'     , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                       transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.04, 0.9)
    ax.set_ylim(-0.2, 5.6)
    ax.set_xticks([0.05,0.1,0.3,0.7])
    ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.3$',r'$0.7$'])

    majorLocator = MultipleLocator(1.00)
    minorLocator = MultipleLocator(0.50)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)

    ax.text(0.02, 0.80,  r'\textrm{\textbf{SLAC}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SLACE143'],thy_plot['SLACE143'],hand['SLACE143']))
    handles.append((thy_band['SLACE155'],thy_plot['SLACE155'],hand['SLACE155']))
    labels.append(r'\textbf{\textrm{E143}}')
    labels.append(r'\textbf{\textrm{E155}}')
    ax.legend(handles,labels,loc=(0.01,0.62), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    #######################################
    #--Plot Atpe deuteron for SLAC
    #######################################

    ax = ax12
    DATA = {}
    DATA['SLACE155x_E29_t2.75'] = pd.DataFrame(data[10031]).query('Elab==29.1 and theta==2.75')
    DATA['SLACE155x_E29_t5.50'] = pd.DataFrame(data[10031]).query('Elab==29.1 and theta==5.5')
    DATA['SLACE155x_E29_t10.5'] = pd.DataFrame(data[10031]).query('Elab==29.1 and theta==10.5')
    DATA['SLACE155x_E32_t2.75'] = pd.DataFrame(data[10031]).query('Elab==32.3 and theta==2.75')
    DATA['SLACE155x_E32_t5.50'] = pd.DataFrame(data[10031]).query('Elab==32.3 and theta==5.5')
    DATA['SLACE155x_E32_t10.5'] = pd.DataFrame(data[10031]).query('Elab==32.3 and theta==10.5')

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'d Atpe')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE155': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    #--plot labels
    ax.text(0.45, 0.72,  r'$Q^2 > 15~{\rm GeV^2} \, (i=4)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.24, 0.54,  r'$7 < Q^2 < 15$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.18, 0.36,  r'$5 < Q^2 < 7$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.02, 0.19,  r'$3 < Q^2 < 5$'                    , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.69, 0.03,  r'$m_c^2 < Q^2 < 3\, (i=0)$'        , transform = ax.transAxes, fontsize = text_fs)

    x, y   = 0.03, 0.86
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$\tilde{A}_{\perp}$' , transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                       transform = ax.transAxes, size = 60)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    ax.semilogx()
    ax.set_xlim(0.04, 0.9)
    ax.set_ylim(-0.2, 5.6)
    ax.set_xticks([0.05,0.1,0.3,0.7])
    ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.3$',r'$0.7$'])

    ax.yaxis.set_major_locator(MultipleLocator(1.00))
    ax.yaxis.set_minor_locator(MultipleLocator(0.50))

    ax.text(0.42, 0.95,  r'\textrm{\textbf{SLAC E29}}'      , transform = ax.transAxes, fontsize = legend_fs)
    ax.text(0.77, 0.95,  r'\textrm{\textbf{E32}}'      , transform = ax.transAxes, fontsize = legend_fs)
    handles,labels = [],[]
    handles.append((thy_band['SLACE155x_E29_t2.75'],thy_plot['SLACE155x_E29_t2.75'],hand['SLACE155x_E29_t2.75']))
    handles.append((thy_band['SLACE155x_E29_t5.50'],thy_plot['SLACE155x_E29_t5.50'],hand['SLACE155x_E29_t5.50']))
    handles.append((thy_band['SLACE155x_E29_t10.5'],thy_plot['SLACE155x_E29_t10.5'],hand['SLACE155x_E29_t10.5']))
    handles.append((thy_band['SLACE155x_E32_t2.75'],thy_plot['SLACE155x_E32_t2.75'],hand['SLACE155x_E32_t2.75']))
    handles.append((thy_band['SLACE155x_E32_t5.50'],thy_plot['SLACE155x_E32_t5.50'],hand['SLACE155x_E32_t5.50']))
    handles.append((thy_band['SLACE155x_E32_t10.5'],thy_plot['SLACE155x_E32_t10.5'],hand['SLACE155x_E32_t10.5']))
    labels.append(r'\boldmath$\theta=2.75$')
    labels.append(r'\boldmath$\theta=5.50$')
    labels.append(r'\boldmath$\theta=10.5$')
    labels.append(r'\boldmath$\theta=2.75$')
    labels.append(r'\boldmath$\theta=5.50$')
    labels.append(r'\boldmath$\theta=10.5$')
    ax.legend(handles,labels,loc=(0.40,0.78), fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0, ncol = 2, columnspacing = 4.5)

    py.title(r'\rm \bf deuteron target',fontsize = 60,y = 1.02, x = 0)


    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_Ape_deuteron'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()


#--neutron/helium plots
def plot_A1_Apa_helium(data):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #######################################
    #--Plot A1 neutron/helium for HERMES and SLAC
    #######################################

    ax = ax11
    DATA = {}
    DATA['HERMES']    = pd.DataFrame(data[10005]) #--HERMES A1 n
    DATA['SLACE142']  = pd.DataFrame(data[10018]) #--SLAC A1 h

    theory = {}
    PLOT = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'h A1')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='HERMES': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.16, 0.01
    ax.text(x,   y    ,r'\boldmath$A_1$'     , transform = ax.transAxes, size = 100)
    #ax.text(x+dx,y+dy ,r'$(+\, i)$',           transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.05, 0.7)
    ax.set_ylim(-0.25, 0.25)
    ax.set_xticks([0.05,0.1,0.2,0.5])
    ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.2$',r'$0.5$'])

    majorLocator = MultipleLocator(0.10)
    minorLocator = MultipleLocator(0.05)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks([-0.2,-0.1,0,0.1,0.2])
    ax.set_yticklabels([r'$-0.2$',r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])

    handles,labels = [],[]
    handles.append((thy_band['HERMES'],thy_plot['HERMES'],hand['HERMES']))
    handles.append((thy_band['SLACE142'],thy_plot['SLACE142'],hand['SLACE142']))
    labels.append(r'\textbf{\textrm{HERMES \boldmath$n$}}')
    labels.append(r'\textbf{\textrm{SLAC (E142) \boldmath$^3$He}}')
    ax.legend(handles,labels,loc=(0.01,0.68), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    ax.set_title(r'\boldmath$^3 {\rm He}/n$ \rm \bf target',fontsize = 60)


    #######################################
    #--Plot Apa helium for SLAC and JLab
    #######################################

    ax = ax12
    DATA = {}
    DATA['JLabHAE06014']     = pd.DataFrame(data[10010]) #--JLab Apa h
    DATA['JLabHAE99117']     = pd.DataFrame(data[10014]) #--JLab Apa h
    DATA['SLACE154']         = pd.DataFrame(data[10025]) #--SLAC Apa h

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'h Apa')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 0.1
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE154': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)

    #--plot labels
    ax.text(0.02, 0.70,  r'$Q^2 > 10~{\rm GeV^2} \, (i=0.3)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.02, 0.50,  r'$5 < Q^2 < 10$'                     , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.80, 0.36,  r'$3 < Q^2 < 5$'                      , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.66, 0.07,  r'$m_c^2 < Q^2 < 3\, (i=0)$'          , transform = ax.transAxes, fontsize = text_fs)

    x, y   = 0.03, 0.89
    dx, dy = 0.17, 0.01
    ax.text(x,   y    ,r'\boldmath$A_{\parallel}$', transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                              transform = ax.transAxes, size = 60)

    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)
    #ax.set_xticks([0.01,0.1])
    #ax.set_xticklabels([r'$0.01$',r'$0.1$'])


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.02, 0.65)
    ax.set_ylim(-0.05, 0.45)

    ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6])
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))

    ax.yaxis.set_major_locator(MultipleLocator(0.10))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_yticks([0,0.1,0.2,0.3,0.4])
    ax.set_yticklabels([r'$0$',r'$0.1$',r'$0.2$',r'$0.3$',r'$0.4$'])

    handles,labels = [],[]
    handles.append((thy_band['SLACE154'],thy_plot['SLACE154'],hand['SLACE154']))
    handles.append((thy_band['JLabHAE06014'],thy_plot['JLabHAE06014'],hand['JLabHAE06014']))
    handles.append((thy_band['JLabHAE99117'],thy_plot['JLabHAE99117'],hand['JLabHAE99117']))
    labels.append(r'\textbf{\textrm{SLAC (E154)}}')
    labels.append(r'\textbf{\textrm{JLab (E06-014)}}')
    labels.append(r'\textbf{\textrm{JLab (E99-117)}}')
    ax.legend(handles,labels,loc=(0.45,0.73), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.title(r'\boldmath{$^3 {\rm He}$} \rm \bf target',fontsize = 60)


    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_A1_Apa_helium'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_A2_Ape_helium(data):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)
    ax12 = py.subplot(nrows, ncols, 2)

    #######################################
    #--Plot A1 neutron/helium for HERMES and SLAC
    #######################################

    ax = ax11
    DATA = {}
    DATA['SLACE142']  = pd.DataFrame(data[10019]) #--SLAC A2 h

    theory = {}
    PLOT = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'h A2')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE142': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.89
    dx, dy = 0.16, 0.01
    ax.text(x,   y    ,r'\boldmath$A_2$'     , transform = ax.transAxes, size = 100)
    #ax.text(x+dx,y+dy ,r'$(+\, i)$',           transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.semilogx()
    ax.set_xlim(0.05, 0.7)
    ax.set_ylim(-0.25, 0.25)
    ax.set_xticks([0.05,0.1,0.2,0.5])
    ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.2$',r'$0.5$'])

    majorLocator = MultipleLocator(0.10)
    minorLocator = MultipleLocator(0.05)
    ax.yaxis.set_major_locator(majorLocator)
    ax.yaxis.set_minor_locator(minorLocator)
    ax.set_yticks([-0.2,-0.1,0,0.1,0.2])
    ax.set_yticklabels([r'$-0.2$',r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])

    handles,labels = [],[]
    handles.append((thy_band['SLACE142'],thy_plot['SLACE142'],hand['SLACE142']))
    labels.append(r'\textbf{\textrm{SLAC (E142)}}')
    ax.legend(handles,labels,loc=(0.01,0.75), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    #######################################
    #--Plot Apa helium for SLAC and JLab
    #######################################

    ax = ax12
    DATA = {}
    DATA['JLabHAE06014']     = pd.DataFrame(data[10011]) #--JLab Ape h
    DATA['JLabHAE99117']     = pd.DataFrame(data[10015]) #--JLab Ape h
    DATA['SLACE154']         = pd.DataFrame(data[10024]) #--SLAC Ape h

    PLOT = {}
    theory = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'h Ape')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 0.2
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='SLACE154': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            X     = PLOT[exp]['X'][key]
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.97,X_thy[0]*1.03]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)

    #--plot labels
    ax.text(0.02, 0.65,  r'$Q^2 > 10~{\rm GeV^2} \, (i=0.6)$'  , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.68, 0.48,  r'$5 < Q^2 < 10$'                     , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.80, 0.30,  r'$3 < Q^2 < 5$'                      , transform = ax.transAxes, fontsize = text_fs)
    ax.text(0.66, 0.05,  r'$m_c^2 < Q^2 < 3\, (i=0)$'          , transform = ax.transAxes, fontsize = text_fs)

    x, y   = 0.03, 0.89
    dx, dy = 0.18, 0.01
    ax.text(x,   y    ,r'\boldmath$A_{\perp}$', transform = ax.transAxes, size = 100)
    ax.text(x+dx,y+dy ,r'$(+\, i)$',                          transform = ax.transAxes, size = 60)

    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)
    #ax.set_xticks([0.01,0.1])
    #ax.set_xticklabels([r'$0.01$',r'$0.1$'])


    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')

    #ax.semilogx()
    ax.set_xlim(0.02, 0.65)
    ax.set_ylim(-0.05, 0.88)

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6])

    ax.yaxis.set_major_locator(MultipleLocator(0.20))
    ax.yaxis.set_minor_locator(MultipleLocator(0.10))
    ax.set_yticks([0,0.2,0.4,0.6,0.8])
    ax.set_yticklabels([r'$0$',r'$0.2$',r'$0.4$',r'$0.6$',r'$0.8$'])

    handles,labels = [],[]
    handles.append((thy_band['SLACE154'],thy_plot['SLACE154'],hand['SLACE154']))
    handles.append((thy_band['JLabHAE06014'],thy_plot['JLabHAE06014'],hand['JLabHAE06014']))
    handles.append((thy_band['JLabHAE99117'],thy_plot['JLabHAE99117'],hand['JLabHAE99117']))
    labels.append(r'\textbf{\textrm{SLAC (E154)}}')
    labels.append(r'\textbf{\textrm{JLab (E06-014)}}')
    labels.append(r'\textbf{\textrm{JLab (E99-117)}}')
    ax.legend(handles,labels,loc=(0.45,0.73), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    py.title(r'\boldmath{$^3 {\rm He}$} \rm \bf target',fontsize = 60,y = 1.02, x = 0)

    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_A2_Ape_helium'%cwd
    filename += ext
    ax11.set_rasterized(True)
    ax12.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()


def plot_JLab_E06_014_Ape(data,data_noHT4):


    nrows, ncols = 1, 2
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    axs = {}
    axs[0] = py.subplot(nrows, ncols, 1)
    axs[1] = py.subplot(nrows, ncols, 2)

    DATA = {}
    DATA['JLabHAE06014 E4.74']     = pd.DataFrame(data[10011]).query("Elab==4.74") #--JLab Ape h
    DATA['JLabHAE06014 E5.89']     = pd.DataFrame(data[10011]).query("Elab==5.89") #--JLab Ape h

    DATA_noHT4 = {}
    DATA_noHT4['JLabHAE06014 E4.74']     = pd.DataFrame(data_noHT4[10011]).query("Elab==4.74") #--JLab Ape h
    DATA_noHT4['JLabHAE06014 E5.89']     = pd.DataFrame(data_noHT4[10011]).query("Elab==5.89") #--JLab Ape h

    j = 0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in DATA:
        X     = DATA[exp]['X']
        val   = DATA[exp]['value']
        alpha = DATA[exp]['alpha']
        #color,marker,ms = get_details(exp)
        color,marker,ms = 'black','o', 8
        hand[exp] = axs[j].errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=4,linestyle='none')

        mean = DATA[exp]['thy-0']
        std  = DATA[exp]['dthy-0']
        down = mean - std
        up   = mean + std
        color = 'firebrick'
        hand['JAM mean'] ,= axs[j].plot(X,mean,linestyle='solid',color=color)
        hand['JAM']  = axs[j].fill_between(X,down,up,color=color,alpha=0.3)

        mean_noHT4 = DATA_noHT4[exp]['thy-0']
        std_noHT4  = DATA_noHT4[exp]['dthy-0']
        down_noHT4 = mean_noHT4 - std_noHT4
        up_noHT4   = mean_noHT4 + std_noHT4
        color = 'gold'
        hand['noHT mean'] ,= axs[j].plot(X,mean_noHT4,linestyle='solid',color=color)
        hand['noHT']  = axs[j].fill_between(X,down_noHT4,up_noHT4,color=color,alpha=0.3)

        j+=1

    for i in range(len(axs)):
        axs[i].tick_params(axis = 'both', which = 'both', labelleft = True, labelright = False, direction='in')
        axs[i].tick_params(axis = 'both', labelsize = label_fs)

        axs[i].yaxis.set_tick_params(which = 'major', length = 10)
        axs[i].yaxis.set_tick_params(which = 'minor', length = 5)

        axs[i].xaxis.set_tick_params(which = 'major', length = 10)
        axs[i].xaxis.set_tick_params(which = 'minor', length = 5)
        axs[i].set_xlim(0.25,0.65)
        axs[i].set_xticks([0.3,0.4,0.5,0.6])
        axs[i].xaxis.set_minor_locator(MultipleLocator(0.05))


        axs[i].set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
        # ax.xaxis.set_label_coords(0.95,0.00)
        axs[i].axhline(0,0,1,color='black',alpha=0.5)

        axs[i].set_ylim(-0.016,0.02)
        axs[i].set_yticks([-0.01,0,0.01])
        axs[i].set_yticklabels([r'$-0.01$',r'$0$',r'$0.01$'])
        axs[i].yaxis.set_minor_locator(MultipleLocator(0.005))

    axs[1].tick_params(labelleft=False)
    #--plot labels
    axs[0].text(0.60, 0.90,  r'\boldmath$E=4.7~{\rm GeV}$'  , transform = axs[0].transAxes, fontsize = 40)
    axs[1].text(0.60, 0.90,  r'\boldmath$E=5.9~{\rm GeV}$'  , transform = axs[1].transAxes, fontsize = 40)

    axs[0].text(0.05,   0.85    ,r'\boldmath$A_{\perp}^{^3{\rm He}}$', transform = axs[0].transAxes, size = 100)


    handles,labels = [],[]
    handles.append(hand['JLabHAE06014 E4.74'])
    handles.append((hand['JAM mean'],hand['JAM'] ))
    handles.append((hand['noHT mean'],hand['noHT']))
    labels.append(r'\textbf{\textrm{JLab (E06-014)}}')
    labels.append(r'\textbf{\textrm{JAM}}')
    labels.append(r'\textbf{\textrm{no HT}}')
    axs[1].legend(handles,labels,loc=(0.00,0.00), fontsize = legend_fs, frameon = 0, handletextpad = 0.3, handlelength = 1.0)


    #py.title(r'\boldmath{$^3 {\rm He}$} \rm \bf target',fontsize = 60,y = 1.02, x = 0)

    py.tight_layout()
    py.subplots_adjust(right=0.99,top=0.99)
    filename = '%s/gallery/fig_JLab_E06_014_Ape'%cwd
    filename += ext
    axs[0].set_rasterized(True)
    axs[1].set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()

def plot_JLab_E12_06_110(data):


    nrows, ncols = 1, 1
    py.figure(figsize = (ncols * 12.0, nrows * 14.0))
    ax11 = py.subplot(nrows, ncols, 1)

    #######################################
    #--Plot A1 neutron/helium for HERMES and SLAC
    #######################################

    ax = ax11
    DATA = {}
    DATA['JLab E12 Apa']  = pd.DataFrame(data[20005]) #--
    DATA['JLab E12 Ape']  = pd.DataFrame(data[20006]) #--

    theory = {}
    PLOT = {}
    for cluster_i in range(kc.nc[istep]):
        if cluster_i != 0: continue

        for exp in DATA:
            query = get_Q2bins(DATA[exp],'h A1')
            PLOT[exp] = get_plot(query)

            nbins = len(query)
            theory[exp] = get_theory(PLOT[exp],nbins,funcX=True,funcQ2=False,loop=False)

    N = 1.0
    hand = {}
    thy_plot = {}
    thy_band = {}
    #--plot data points
    for exp in PLOT:
        for key in range(nbins):
            if exp=='JLab E12 Apa': ax.axhline(N*key,0,1,color='black',alpha=0.1)
            if exp=='JLab E12 Apa': shift = 0.0 
            if exp=='JLab E12 Ape': shift = 0.005
            X     = PLOT[exp]['X'][key] + shift
            val   = PLOT[exp]['value'][key] + N*float(key)
            alpha = PLOT[exp]['alpha'][key]
            color,marker,ms = get_details(exp)
            hand[exp] = ax.errorbar(X,val,alpha,marker=marker,color=color,ms=ms,capsize=3,linestyle='none')

            #--plot theory interpolated between data points
            if key not in theory[exp]['value']: continue
            X_thy= theory[exp]['X'][key]
            mean = theory[exp]['value'][key] + N*float(key)
            std  = theory[exp]['std'][key]
            down = mean - std
            up   = mean + std
            thy_plot[exp] ,= ax.plot(X_thy,mean,linestyle='solid',color=color)
            thy_band[exp]  = ax.fill_between(X_thy,down,up,color=color,alpha=0.3)
            if len(X)==1:
                x = [X_thy[0]*0.95,X_thy[0]*1.05]
                down, up, mean = [down[0], down[0]], [up[0],up[0]], [mean[0],mean[0]]
                thy_band[exp] = ax.fill_between(x,down,up,color=color,alpha=0.3)
                thy_plot[exp] ,= ax.plot(x,mean,linestyle='solid',color=color)


    ax.tick_params(axis = 'both', labelsize = label_fs)

    ax.yaxis.set_tick_params(which = 'major', length = 10)
    ax.yaxis.set_tick_params(which = 'minor', length = 5)

    ax.xaxis.set_tick_params(which = 'major', length = 10)
    ax.xaxis.set_tick_params(which = 'minor', length = 5)

    x, y   = 0.03, 0.91
    dx, dy = 0.16, 0.01
    ax.text(x,   y       ,r'\textrm{\textbf{JLab E12-06-110}}' , transform = ax.transAxes, size = 50)
    ax.text(0.03,0.82    ,r'\textrm{\textbf{PRELIMINARY}}'     , transform = ax.transAxes, size = 50)
    #ax.text(x+dx,y+dy ,r'$(+\, i)$',           transform = ax.transAxes, size = 60)

    # ax.set_xlabel(r'\boldmath$x_{\rm bj}$', size=xaxis_fs)
    ax.set_xlabel(r'\boldmath$x$', size=xaxis_fs)
    # ax.xaxis.set_label_coords(0.95,0.00)
    ax.tick_params(axis='x',pad=8)

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.tick_params(axis = 'y', which = 'both', labelleft = True, labelright = False, direction='in')
    ax.tick_params(axis = 'x', which = 'both', labeltop = False, labelbottom = True, direction='in')


    ax.set_xlim(0.41, 0.79)
    #ax.set_xticks([0.05,0.1,0.2,0.5])
    #ax.set_xticklabels([r'$0.05$',r'$0.1$',r'$0.2$',r'$0.5$'])

    ax.set_ylim(-0.08, 0.08)
    #ax.yaxis.set_major_locator(MultipleLocator(0.10))
    #ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    #ax.set_yticks([-0.2,-0.1,0,0.1,0.2])
    #ax.set_yticklabels([r'$-0.2$',r'$-0.1$',r'$0$',r'$0.1$',r'$0.2$'])

    handles,labels = [],[]
    handles.append((thy_band['JLab E12 Apa'],thy_plot['JLab E12 Apa'],hand['JLab E12 Apa']))
    handles.append((thy_band['JLab E12 Ape'],thy_plot['JLab E12 Ape'],hand['JLab E12 Ape']))
    labels.append(r'\boldmath$A_{\parallel}$')
    labels.append(r'\boldmath$A_{\perp}$')
    ax.legend(handles,labels,loc=(0.00,0.00), fontsize = 50, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.title(r'\boldmath{$^3 {\rm He}$} \rm \bf target',fontsize = 60)


    py.tight_layout()
    py.subplots_adjust(right=0.99)
    filename = '%s/gallery/fig_pidis_JLab_E12_06_110'%cwd
    filename += ext
    ax11.set_rasterized(True)
    py.savefig(filename)
    print('saving figure to %s'%filename)
    py.close()


if __name__ == "__main__":

    #wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4/pos_g/'
    #wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/xz_predict'
    wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/jlab12/pos_g'

    print('\nplotting pidis data from %s' % (wdir))

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
        for ic in range(kc.nc[istep]):
            #predictions_ic = [predictions[i] for i in range(len(predictions)) if cluster[i] == ic]
            predictions_ic = predictions
            data[idx]['thy-%d' % ic] = np.mean(predictions_ic, axis = 0)
            data[idx]['dthy-%d' % ic] = np.std(predictions_ic, axis = 0)
            if 'X' in data[idx]: data[idx]['x'] = data[idx]['X']
            data[idx]['rQ2'] = np.around(data[idx]['Q2'], decimals = 0)
            data[idx]['rx'] = np.around(data[idx]['x'], decimals = 2)

    wdir_noHT4 = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/noHT4/pos_g/'

    load_config('%s/input.py' % wdir_noHT4)
    istep = core.get_istep()
    replicas = core.get_replicas(wdir_noHT4)
    core.mod_conf(istep, replicas[0]) #--set conf as specified in istep

    predictions = load('%s/data/predictions-%d.dat' % (wdir_noHT4, istep))

    data_noHT4 = predictions['reactions']['pidis']

    for idx in data_noHT4:
        predictions = copy.copy(data_noHT4[idx]['prediction-rep'])
        del data_noHT4[idx]['prediction-rep']
        del data_noHT4[idx]['residuals-rep']
        del data_noHT4[idx]['shift-rep']
        del data_noHT4[idx]['rres-rep']
        del data_noHT4[idx]['r-residuals']
        del data_noHT4[idx]['n-residuals']
        for ic in range(kc.nc[istep]):
            #predictions_ic = [predictions[i] for i in range(len(predictions)) if cluster[i] == ic]
            predictions_ic = predictions
            data_noHT4[idx]['thy-%d' % ic] = np.mean(predictions_ic, axis = 0)
            data_noHT4[idx]['dthy-%d' % ic] = np.std(predictions_ic, axis = 0)
            if 'X' in data_noHT4[idx]: data_noHT4[idx]['x'] = data_noHT4[idx]['X']
            data_noHT4[idx]['rQ2'] = np.around(data_noHT4[idx]['Q2'], decimals = 0)
            data_noHT4[idx]['rx'] = np.around(data_noHT4[idx]['x'], decimals = 2)

    plot_JLab_E12_06_110(data)
    sys.exit()

    plot_A1_proton (data)
    plot_Apa_proton(data)

    plot_A2_proton (data) 
    plot_Ape_proton(data) 


    plot_A1_deuteron (data)
    plot_Apa_deuteron(data)

    plot_Ape_deuteron(data) 

    plot_A1_Apa_helium(data)
    plot_A2_Ape_helium(data)


