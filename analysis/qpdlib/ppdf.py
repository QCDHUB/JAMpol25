import sys, os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MultipleLocator
from matplotlib import cm
import pylab as py


#--from scipy stack 
from scipy.integrate import quad
#from scipy.integrate import cumulative_trapezoid as cumtrapz
from scipy.integrate import cumulative_trapezoid

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core
from analysis.corelib import classifier

import kmeanconf as kc

#--calculate other polarized PDFs
#from qcdlib.qpdcalc import QPDCALC
#try: 
#    from analysis.qpdlib.sets.DSSV import dssvlib
#    from analysis.qpdlib.sets.DSSV.DSSVcalc import DSSV
#except:
#    pass

flavors = []
flavors.append('up')
flavors.append('dp')
flavors.append('u')
flavors.append('d')
flavors.append('sp')
flavors.append('sm')
flavors.append('s')
flavors.append('sb')
flavors.append('g')
flavors.append('uv')
flavors.append('dv')
flavors.append('ub')
flavors.append('db')
flavors.append('ub+db')
flavors.append('ub-db')

cmap = matplotlib.cm.get_cmap('plasma')

def gen_xf(wdir,Q2 = None):
 
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    if Q2==None: Q2 = conf['Q20']

    _replicas = core.get_replicas(wdir)
    names    = core.get_replicas_names(wdir)
    core.mod_conf(istep,_replicas[0])

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        if 'ppdf' not in conf['steps'][istep]['passive distributions']:
            print('ppdf-proton not an active or passive distribution')
            return

    if 'ppdf' in conf['steps'][istep]['active distributions']:
        passive = False
    else:
        passive = True

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = _replicas[0]['order'][istep]

    jar = load('%s/data/jar-%d.dat'%(wdir,istep))
    replicas = jar['replicas']

    ppdf = conf['ppdf']

    ## setup kinematics
    X=10**np.linspace(-6,-1,200)
    X=np.append(X,np.linspace(0.101,0.99,200))
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating polarized pdf-proton from %s at Q2 = %f' % (wdir, Q2))

    ## compute XF for all replicas
    XF = {}
    n_replicas = len(replicas)

    cnt = 0
    for par in replicas:
        lprint('%d/%d' % (cnt, n_replicas))

        if passive: core.mod_conf(istep, _replicas[cnt])
        parman.set_new_params(par, initial = True)

        for flavor in flavors:
            if flavor not in XF: XF[flavor] = []
            if flavor == 'up':
                func = lambda x: ppdf.get_xF(x, Q2, 'u') + ppdf.get_xF(x, Q2, 'ub')
            elif flavor == 'dp':
                func = lambda x: ppdf.get_xF(x, Q2, 'd') + ppdf.get_xF(x, Q2, 'db')
            elif flavor == 'u':
                func = lambda x: ppdf.get_xF(x, Q2, 'u')
            elif flavor == 'd':
                func = lambda x: ppdf.get_xF(x, Q2, 'd')
            elif flavor == 'sp':
                func = lambda x: ppdf.get_xF(x, Q2, 's') + ppdf.get_xF(x, Q2, 'sb')
            elif flavor == 'ub+db':
                func = lambda x: ppdf.get_xF(x, Q2, 'ub') + ppdf.get_xF(x, Q2, 'db')
            elif flavor == 'ub-db':
                func = lambda x: ppdf.get_xF(x, Q2, 'ub') - ppdf.get_xF(x, Q2, 'db')
            else:
                func = lambda x: ppdf.get_xF(x, Q2, flavor)

            XF[flavor].append(np.array([func(x) for x in X]))

        cnt +=1

    print()
    checkdir('%s/data' % wdir)
    filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
    save({'X': X, 'Q2': Q2, 'XF': XF}, filename)
    print('Saving data to %s'%filename)
        
def plot_xf_main(wdir,Q2=None,mode=0,SETS={}):
  #--mode 0: plot each replica
  #--mode 1: plot average and standard deviation of replicas 


  nrows,ncols=3,2
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*7,nrows*4))
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
  if Q2==None: Q2 = conf['Q20']
  istep=core.get_istep()
  if 'ppdf' not in conf['steps'][istep]['active distributions']:
      if 'ppdf' not in conf['steps'][istep]['passive distributions']:
          print('ppdf-proton not an active or passive distribution')
          return

  filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
  try:
      data=load(filename)
  except:
      gen_xf(wdir,Q2)
      data=load(filename)

  replicas=core.get_replicas(wdir)
  #cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
  #best_cluster=cluster_order[0]

  #scale = classifier.get_scale(wdir)

  X=data['X']

  SAVE = {}
  SAVE['X'] = X
  SAVE['Q2'] = Q2*np.ones(len(X))

  u_m_d   = np.array(data['XF']['u'] ) - np.array(data['XF']['d'] )
  ub_m_db = np.array(data['XF']['ub']) - np.array(data['XF']['db'])

  u_m_d_mean   = np.mean(u_m_d  ,axis=0)
  u_m_d_std    = np.std (u_m_d  ,axis=0)
  ub_m_db_mean = np.mean(ub_m_db,axis=0)
  ub_m_db_std  = np.std (ub_m_db,axis=0)

  SAVE['u minus d mean']   = u_m_d_mean
  SAVE['u minus d std']    = u_m_d_std
  SAVE['ub minus db mean'] = ub_m_db_mean
  SAVE['ub minus db std']  = ub_m_db_std

  SAVE = pd.DataFrame(SAVE)

  #SAVE.to_excel('Jinchen.xlsx')
  
  for flav in data['XF']:

      mean = np.mean(data['XF'][flav],axis=0)
      std  = np.std (data['XF'][flav],axis=0)

    
      y = 5 
      perc_min = np.percentile(data['XF'][flav],y,    axis=0) 
      perc_max = np.percentile(data['XF'][flav],100-y,axis=0) 

      if flav=='up'      : ax,axL = axs[1],axLs[1]
      elif flav=='dp'    : ax,axL = axs[2],axLs[2]
      elif flav=='g'     : ax,axL = axs[3],axLs[3]
      elif flav=='sp'    : ax,axL = axs[4],axLs[4]
      elif flav=='ub'    : ax,axL = axs[5],axLs[5]
      elif flav=='db'    : ax,axL = axs[6],axLs[6]
      else: continue


      #--plot each replica
      if mode==0:
          for i in range(len(data['XF'][flav])):
              thy ,= ax.plot(X,np.array(data['XF'][flav][i]),color='red',alpha=0.3,zorder=2)
              axL.      plot(X,np.array(data['XF'][flav][i]),color='red',alpha=0.3,zorder=2)
    
      #--plot average and standard deviation
      if mode==1:
          thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)
          axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)

      #--plot percentiles
      if mode==2:
          thy = ax.fill_between(X,perc_min,perc_max,color='red',alpha=0.6,zorder=2)
          axL.     fill_between(X,perc_min,perc_max,color='red',alpha=0.6,zorder=2)

      #--plot other PPDF sets
      for SET in SETS:
          _SET, _color, _alpha = SETS[SET][0], SETS[SET][1], SETS[SET][2]

          if SET=='DSSV':
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
              hand[SET] = ax.fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=1)

          if SET=='NNPDF' or SET=='JAM17':
              ppdf = _SET.get_xpdf(flav,X,Q2)
              hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=1)



  for i in range(N):
        axs[i+1].set_xlim(8e-3,0.1)
        axs[i+1].semilogx()

        axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=10)
        axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=5)
        axs[i+1].set_xticks([0.01,0.1])
        axs[i+1].set_xticklabels([r'$0.01$',r'$0.1$'])

        axLs[i+1].set_xlim(0.1,1.0)

        axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=10)
        axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=5)
        axLs[i+1].set_xticks([0.3,0.5,0.7])
        axLs[i+1].set_xticklabels([r'$0.3$',r'$0.5$',r'$0.7$'])

  for i in [1,2,3,4]:
        axs[i] .tick_params(labelbottom=False)
        axLs[i].tick_params(labelbottom=False)

  axs[1].set_ylim(-0.10,0.50)   
  axs[2].set_ylim(-0.20,0.10)   
  axs[3].set_ylim(-0.30,0.40)  
  axs[4].set_ylim(-0.07,0.15)
  axs[5].set_ylim(-0.07,0.15)
  axs[6].set_ylim(-0.07,0.15)

  axs[1].set_yticks([0,0.2,0.4])
  axs[2].set_yticks([-0.15,-0.10,-0.05,0,0.05])
  axs[3].set_yticks([-0.2,0.0,0.2])
  axs[4].set_yticks([-0.04,0.00,0.04,0.08,0.12])
  axs[5].set_yticks([-0.04,0.00,0.04,0.08,0.12])
  axs[6].set_yticks([-0.04,0.00,0.04,0.08,0.12])

  for i in range(N):
      axs [i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axLs[i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
      axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

  axLs[5].set_xlabel(r'\boldmath$x$',size=30)
  axLs[6].set_xlabel(r'\boldmath$x$',size=30)   
  axLs[5].xaxis.set_label_coords(0.95,0.00)
  axLs[6].xaxis.set_label_coords(0.95,0.00)

  axs[1].text(0.10,0.80,r'\boldmath{$x \Delta u^+$}',                        transform=axs[1].transAxes,size=30)
  axs[2].text(0.10,0.80,r'\boldmath{$x \Delta d^+$}',                        transform=axs[2].transAxes,size=30)
  axs[3].text(0.10,0.85,r'\boldmath{$x \Delta g$}'  ,                        transform=axs[3].transAxes,size=30)
  axs[4].text(0.10,0.80,r'\boldmath{$x \Delta s^+$}',                        transform=axs[4].transAxes,size=30)
  axs[5].text(0.10,0.80,r'\boldmath{$x \Delta \bar{u}$}',                    transform=axs[5].transAxes,size=30)
  axs[6].text(0.10,0.80,r'\boldmath{$x \Delta \bar{d}$}',                    transform=axs[6].transAxes,size=30)

  if Q2 == 1.27**2: axs[1].text(0.10,0.65,r'$Q^2 = m_c^2$',                            transform=axs[1].transAxes,size=30)
  else:             axs[1].text(0.10,0.65,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axs[1].transAxes,size=30)

  if len(SETS) > 0:
      handles,labels = [],[]
      for _ in SETS:
          handles.append(hand[_])
          labels.append(SETS[_][3])

      #axLs[2].legend(handles,labels,loc=(-0.30,0.80),fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
      axLs[2].legend(handles,labels,loc='upper right',fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
    
  if mode==0:
      sm   = py.cm.ScalarMappable(cmap=cmap)
      sm.set_array([])
      cax = fig.add_axes([0.77,0.95,0.20,0.02])
      cax.tick_params(axis='both',which='both',labelsize=20,direction='in')
      cax.xaxis.set_label_coords(0.65,-1.4)
      cbar = py.colorbar(sm,cax=cax,orientation='horizontal',ticks=[0.2,0.4,0.6,0.8])
      cbar.set_label(r'\boldmath${\rm scaled}~\chi^2_{\rm red}$',size=30)
 
  py.tight_layout()
  py.subplots_adjust(wspace=0.2,hspace=0)

  filename = '%s/gallery/ppdfs-Q2=%3.5f'%(wdir,Q2)
  if mode==1: filename+='-bands'
  if mode==2: filename+='-percentile'

  filename+='.png'

  checkdir('%s/gallery'%wdir)
  py.savefig(filename)
  print ('Saving figure to %s'%filename)

def plot_lowx(wdir,Q2=None,mode=0,SETS={}):
  #--mode 0: plot each replica
  #--mode 1: plot average and standard deviation of replicas 


  nrows,ncols=3,2
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*7,nrows*4))
  axs,axLs = {},{}
  for i in range(N):
      axs[i+1] = py.subplot(nrows,ncols,i+1)
      divider = make_axes_locatable(axs[i+1])
      axLs[i+1] = divider.append_axes("right",size=0.00,pad=0,sharey=axs[i+1])
      axLs[i+1].set_xlim(0.1,0.9)
      axLs[i+1].spines['left'].set_visible(False)
      axLs[i+1].yaxis.set_ticks_position('right')
      py.setp(axLs[i+1].get_xticklabels(),visible=True)

      axs[i+1].spines['right'].set_visible(False)


  thy  = {}
  hand = {}
  load_config('%s/input.py'%wdir)
  if Q2==None: Q2 = conf['Q20']
  istep=core.get_istep()
  if 'ppdf' not in conf['steps'][istep]['active distributions']:
      if 'ppdf' not in conf['steps'][istep]['passive distributions']:
          print('ppdf-proton not an active or passive distribution')
          return

  filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
  try:
      data=load(filename)
  except:
      gen_xf(wdir,Q2)
      data=load(filename)

  replicas=core.get_replicas(wdir)
  cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
  best_cluster=cluster_order[0]

  scale = classifier.get_scale(wdir)

  X=data['X']

  for flav in data['XF']:
      mean = np.mean(data['XF'][flav]/X,axis=0)
      std  = np.std (data['XF'][flav]/X,axis=0)

      if flav=='up'      : ax,axL = axs[1],axLs[1]
      elif flav=='dp'    : ax,axL = axs[2],axLs[2]
      elif flav=='g'     : ax,axL = axs[3],axLs[3]
      elif flav=='sp'    : ax,axL = axs[4],axLs[4]
      elif flav=='ub'    : ax,axL = axs[5],axLs[5]
      elif flav=='db'    : ax,axL = axs[6],axLs[6]
      else: continue


      #--plot each replica
      if mode==0:
          for i in range(len(data['XF'][flav])):
              #--if plotting one step, use clusters
              thy ,= ax.plot(X,np.array(data['XF'][flav][i])/X,color='red',alpha=0.3,zorder=2)
              axL.      plot(X,np.array(data['XF'][flav][i])/X,color='red',alpha=0.3,zorder=2)
    
      #--plot average and standard deviation
      if mode==1:
          thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)
          axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)

      #--plot other PPDF sets
      for SET in SETS:
          _SET, _color, _alpha = SETS[SET][0], SETS[SET][1], SETS[SET][2]

          if SET=='DSSV':
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
              hand[SET] = ax.fill_between(_X,(mean-std)/_X,(mean+std)/_X,color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(_X,(mean-std)/_X,(mean+std)/_X,color=_color,alpha=_alpha,zorder=1)

          if SET=='NNPDF' or SET=='JAM17':
              ppdf = _SET.get_xpdf(flav,X,Q2)
              hand[SET] = ax.fill_between(X,ppdf['xfmin']/X,ppdf['xfmax']/X,color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(X,ppdf['xfmin']/X,ppdf['xfmax']/X,color=_color,alpha=_alpha,zorder=1)



  for i in range(N):
        axs[i+1].set_xlim(8e-5,1.0)
        axs[i+1].semilogx()

        axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,pad=10,length=10)
        axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,pad=10,length=5)
        axs[i+1].set_xticks([0.0001,0.001,0.01,0.1,1.0])
        axs[i+1].set_xticklabels([r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$10^{-1}$',r'$1$'])

        axLs[i+1].set_xlim(0.99,1.0)

        axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=10)
        axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=5)
        axLs[i+1].set_xticks([])

  for i in [1,2,3,4]:
        axs[i] .tick_params(labelbottom=False)
        axLs[i].tick_params(labelbottom=False)

    
  axs[1].set_ylim(-35,35)   
  axs[2].set_ylim(-35,35)   
  axs[3].set_ylim(-35,35)  
  axs[4].set_ylim(-35,35)
  axs[5].set_ylim(-35,35)
  axs[6].set_ylim(-35,35)

  #axs[1].set_yticks([0,0.2,0.4])
  #axs[2].set_yticks([-0.15,-0.10,-0.05,0,0.05])
  #axs[3].set_yticks([-0.2,0.0,0.2])
  #axs[4].set_yticks([-0.40,-0.20,0,0.20])
  #axs[5].set_yticks([-0.04,0.00,0.04,0.08,0.12])
  #axs[6].set_yticks([-0.04,0.00,0.04,0.08,0.12])

  for i in range(N):
      axs [i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axLs[i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      #axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
      #axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

  axs[5].set_xlabel(r'\boldmath$x$',size=40)
  axs[6].set_xlabel(r'\boldmath$x$',size=40)   
  axs[5].xaxis.set_label_coords(0.95,-0.02)
  axs[6].xaxis.set_label_coords(0.95,-0.02)

  axs[1].text(0.80,0.10,r'\boldmath{$\Delta u^+$}',                        transform=axs[1].transAxes,size=30)
  axs[2].text(0.80,0.10,r'\boldmath{$\Delta d^+$}',                        transform=axs[2].transAxes,size=30)
  axs[3].text(0.80,0.15,r'\boldmath{$\Delta g$}'  ,                        transform=axs[3].transAxes,size=30)
  axs[4].text(0.80,0.10,r'\boldmath{$\Delta s^+$}',                        transform=axs[4].transAxes,size=30)
  axs[5].text(0.80,0.10,r'\boldmath{$\Delta \bar{u}$}',                    transform=axs[5].transAxes,size=30)
  axs[6].text(0.80,0.10,r'\boldmath{$\Delta \bar{d}$}',                    transform=axs[6].transAxes,size=30)

  if Q2 == 1.27**2: axs[5].text(0.50,0.80,r'$Q^2 = m_c^2$',                            transform=axs[5].transAxes,size=30)
  else:             axs[5].text(0.50,0.80,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axs[5].transAxes,size=30)

  if len(SETS) > 0:
      handles,labels = [],[]
      for _ in SETS:
          handles.append(hand[_])
          labels.append(SETS[_][3])

      #axLs[2].legend(handles,labels,loc=(-0.30,0.80),fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
      axLs[2].legend(handles,labels,loc='upper right',fontsize=28,frameon=0,handletextpad=0.3,handlelength=1.0)
    
  if mode==0:
      sm   = py.cm.ScalarMappable(cmap=cmap)
      sm.set_array([])
      cax = fig.add_axes([0.77,0.95,0.20,0.02])
      cax.tick_params(axis='both',which='both',labelsize=20,direction='in')
      cax.xaxis.set_label_coords(0.65,-1.4)
      cbar = py.colorbar(sm,cax=cax,orientation='horizontal',ticks=[0.2,0.4,0.6,0.8])
      cbar.set_label(r'\boldmath${\rm scaled}~\chi^2_{\rm red}$',size=30)
 
  py.tight_layout()
  py.subplots_adjust(wspace=0.2,hspace=0)

  filename = '%s/gallery/ppdfs_lowx-Q2=%3.5f'%(wdir,Q2)
  if mode==1: filename+='-bands'
  filename+='.png'

  checkdir('%s/gallery'%wdir)
  py.savefig(filename)
  print ('Saving figure to %s'%filename)
       
def plot_asymmetry(wdir,Q2=None,mode=0,SETS={},reweight=None):
  #--mode 0: plot each replica
  #--mode 1: plot average and standard deviation of replicas 

  nrows,ncols=1,1
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*9,nrows*5))
  axs,axLs = {},{}
  for i in range(N):
      axs[i+1] = py.subplot(nrows,ncols,i+1)
      divider = make_axes_locatable(axs[i+1])
      axLs[i+1] = divider.append_axes("right",size=3.5,pad=0,sharey=axs[i+1])
      axLs[i+1].set_xlim(0.1,0.9)
      axLs[i+1].spines['left'].set_visible(False)
      axLs[i+1].yaxis.set_ticks_position('right')
      py.setp(axLs[i+1].get_xticklabels(),visible=True)

      axs[i+1].spines['right'].set_visible(False)


  if reweight!=None:
      filename+='-reweight-%s'%reweight
      try: w    = np.load('%s/data/weights-%s.npy'%(wdir,reweight))
      except:
          print('weights not found')
          w = np.ones(len(data['XF'][flav]))/len(data['XF'][flav])

  thy  = {}
  hand = {}
  load_config('%s/input.py'%wdir)
  istep=core.get_istep()
  if 'ppdf' not in conf['steps'][istep]['active distributions']:
      if 'ppdf' not in conf['steps'][istep]['passive distributions']:
          print('ppdf-proton not an active or passive distribution')
          return

  filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
  try:
      data=load(filename)
  except:
      gen_xf(wdir,Q2)
      data=load(filename)

  replicas=core.get_replicas(wdir)
  cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
  best_cluster=cluster_order[0]

  scale = classifier.get_scale(wdir)
  
  X=data['X']

  for flav in data['XF']:
      if reweight==None:
          mean = np.mean(data['XF'][flav],axis=0)
          std  = np.std (data['XF'][flav],axis=0)
      else:
          mean = np.einsum('i,ij->j',w,data['XF'][flav])
          std  = np.einsum('i,ij->j',w,(data['XF'][flav]-mean)**2)**0.5

      if flav=='ub-db' : ax,axL = axs[1],axLs[1]
      else: continue


      #--plot each replica
      if mode==0:
          for i in range(len(data['XF'][flav])):
              thy ,= ax.plot(X,np.array(data['XF'][flav][i]),color='red',alpha=0.3,zorder=2)
              axL.      plot(X,np.array(data['XF'][flav][i]),color='red',alpha=0.3,zorder=2)
    
      #--plot average and standard deviation
      if mode==1:
          thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)
          axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)

      #--plot other PPDF sets
      for SET in SETS:
          _SET, _color, _alpha = SETS[SET][0], SETS[SET][1], SETS[SET][2]

          if SET=='DSSV':
              _X=10**np.linspace(-4,-1,200)
              _X=np.append(_X,np.linspace(0.1,0.99,200))
              ub = _SET.xfxQ2(-2,_X,Q2)
              db = _SET.xfxQ2(-1,_X,Q2)
              ppdf = ub - db
              mean = ppdf[0]
              std  = 0
              for i in range(1,20):
                  std += (ppdf[i] - ppdf[-i])**2
              std = np.sqrt(std)/2.0
              hand[SET] = ax.fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=1)

          if SET=='NNPDF' or SET=='JAM17':
              ppdf = _SET.get_xpdf(flav,X,Q2)
              hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=1)


  for i in range(N):
        axs[i+1].set_xlim(8e-3,0.1)
        axs[i+1].semilogx()

        axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=6)
        axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=3)
        axs[i+1].set_xticks([0.01,0.1])
        axs[i+1].set_xticklabels([r'$0.01$',r'$0.1$'])

        axLs[i+1].set_xlim(0.1,1.0)

        axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=6)
        axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=3)
        axLs[i+1].set_xticks([0.3,0.5,0.7])
        axLs[i+1].set_xticklabels([r'$0.3$',r'$0.5$',r'$0.7$'])

  axs[1].set_ylim(-0.05,0.11)

  minorLocator = MultipleLocator(0.02)
  majorLocator = MultipleLocator(0.04)
  axs[1].yaxis.set_minor_locator(minorLocator)
  axs[1].yaxis.set_major_locator(majorLocator)

  for i in range(N):
      axs [i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axLs[i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
      axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

  axLs[1].set_xlabel(r'\boldmath$x$',size=30)   
  axLs[1].xaxis.set_label_coords(0.95,0.00)

  axs[1].text(0.07,0.83,r'\boldmath{$x (\Delta \bar{u} - \Delta \bar{d})$}', transform=axs[1].transAxes,size=40)

  if Q2 == 1.27**2: axs[1].text(0.07,0.07,r'$Q^2 = m_c^2$',                            transform=axs[1].transAxes,size=40)
  else:             axs[1].text(0.07,0.07,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axs[1].transAxes,size=40)

  handles,labels=[],[]
  if len(SETS) > 0:
      for _ in SETS:
          handles.append(hand[_])
          labels.append(SETS[_][3])
    
  axLs[1].legend(handles,labels,loc='upper right',fontsize=25,frameon=0,handletextpad=0.3,handlelength=1.0)
  
  if mode==0:
      sm   = py.cm.ScalarMappable(cmap=cmap)
      sm.set_array([])
      cax = fig.add_axes([0.67,0.87,0.25,0.03])
      cax.tick_params(axis='both',which='both',labelsize=20,direction='in')
      cax.xaxis.set_label_coords(0.65,-2.5)
      cbar = py.colorbar(sm,cax=cax,orientation='horizontal',ticks=[0.2,0.4,0.6,0.8])
      cbar.set_label(r'\boldmath${\rm scaled}~\chi^2_{\rm red}$',size=30)
 
  py.tight_layout()

  filename = '%s/gallery/ppdfs-asymmetry-Q2=%3.5f'%(wdir,Q2)
  if mode==1: filename+='-bands'
  filename+='.png'

  checkdir('%s/gallery'%wdir)
  py.savefig(filename)
  print ('Saving figure to %s'%filename)

def plot_s_asymmetry(wdir,Q2=None,mode=0,SETS={}):
  #--mode 0: plot each replica
  #--mode 1: plot average and standard deviation of replicas 


  nrows,ncols=1,1
  N = nrows*ncols
  fig = py.figure(figsize=(ncols*9,nrows*5))
  axs,axLs = {},{}
  for i in range(N):
      axs[i+1] = py.subplot(nrows,ncols,i+1)
      divider = make_axes_locatable(axs[i+1])
      axLs[i+1] = divider.append_axes("right",size=3.5,pad=0,sharey=axs[i+1])
      axLs[i+1].set_xlim(0.1,0.9)
      axLs[i+1].spines['left'].set_visible(False)
      axLs[i+1].yaxis.set_ticks_position('right')
      py.setp(axLs[i+1].get_xticklabels(),visible=True)

      axs[i+1].spines['right'].set_visible(False)

  thy  = {}
  hand = {}
  load_config('%s/input.py'%wdir)
  if Q2==None: Q2 = conf['Q20']
  istep=core.get_istep()
  if 'ppdf' not in conf['steps'][istep]['active distributions']:
      if 'ppdf' not in conf['steps'][istep]['passive distributions']:
          print('ppdf-proton not an active or passive distribution')
          return

  filename = '%s/data/ppdf-Q2=%3.5f.dat'%(wdir,Q2)
  try:
      data=load(filename)
  except:
      gen_xf(wdir,Q2)
      data=load(filename)

  replicas=core.get_replicas(wdir)
  cluster,colors,nc,cluster_order = classifier.get_clusters(wdir,istep,kc) 
  best_cluster=cluster_order[0]

  scale = classifier.get_scale(wdir)
  
  X=data['X']

  for flav in data['XF']:
      mean = np.mean(data['XF'][flav],axis=0)
      std = np.std(data['XF'][flav],axis=0)

      if flav=='sm' : ax,axL = axs[1],axLs[1]
      else: continue


      #--plot each replica
      if mode==0:
          for i in range(len(data['XF'][flav])):
              thy ,= ax.plot(X,np.array(data['XF'][flav][i]),color='red',alpha=0.3,zorder=2)
              axL.      plot(X,np.array(data['XF'][flav][i]),color='red',alpha=0.3,zorder=2)
    
      #--plot average and standard deviation
      if mode==1:
          thy = ax.fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)
          axL.     fill_between(X,(mean-std),(mean+std),color='red',alpha=0.6,zorder=2)

      #--plot other PPDF sets
      for SET in SETS:
          _SET, _color, _alpha = SETS[SET][0], SETS[SET][1], SETS[SET][2]

          if SET=='DSSV':
              _X=10**np.linspace(-4,-1,200)
              _X=np.append(_X,np.linspace(0.1,0.99,200))
              s  = _SET.xfxQ2( 3,_X,Q2)
              sb = _SET.xfxQ2(-3,_X,Q2)
              ppdf = s - sb
              mean = ppdf[0]
              std  = 0
              for i in range(1,20):
                  std += (ppdf[i] - ppdf[-i])**2
              std = np.sqrt(std)/2.0
              hand[SET] = ax.fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(_X,mean-std,mean+std,color=_color,alpha=_alpha,zorder=1)

          if SET=='NNPDF' or SET=='JAM17':
              ppdf = _SET.get_xpdf(flav,X,Q2)
              hand[SET] = ax.fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=1)
              axL.           fill_between(X,ppdf['xfmin'],ppdf['xfmax'],color=_color,alpha=_alpha,zorder=1)


  for i in range(N):
        axs[i+1].set_xlim(8e-3,0.1)
        axs[i+1].semilogx()

        axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=6)
        axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=3)
        axs[i+1].set_xticks([0.01,0.1])
        axs[i+1].set_xticklabels([r'$0.01$',r'$0.1$'])

        axLs[i+1].set_xlim(0.1,1.0)

        axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=6)
        axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=3)
        axLs[i+1].set_xticks([0.3,0.5,0.7])
        axLs[i+1].set_xticklabels([r'$0.3$',r'$0.5$',r'$0.7$'])

  axs[1].set_ylim(-0.30,0.30)

  minorLocator = MultipleLocator(0.05)
  majorLocator = MultipleLocator(0.10)
  axs[1].yaxis.set_minor_locator(minorLocator)
  axs[1].yaxis.set_major_locator(majorLocator)

  for i in range(N):
      axs [i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axLs[i+1].axhline(0  ,color='k',linestyle='--',alpha=0.5)
      axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
      axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

  axLs[1].set_xlabel(r'\boldmath$x$',size=30)   
  axLs[1].xaxis.set_label_coords(0.95,0.00)

  axs[1].text(0.07,0.83,r'\boldmath{$x (\Delta s - \Delta \bar{s})$}', transform=axs[1].transAxes,size=40)

  if Q2 == 1.27**2: axs[1].text(0.07,0.07,r'$Q^2 = m_c^2$',                            transform=axs[1].transAxes,size=40)
  else:             axs[1].text(0.07,0.07,r'$Q^2 = %s$~'%Q2 + r'\textrm{GeV}'+r'$^2$', transform=axs[1].transAxes,size=40)

  handles,labels=[],[]
  if len(SETS) > 0:
      for _ in SETS:
          handles.append(hand[_])
          labels.append(SETS[_][3])
    
  axLs[1].legend(handles,labels,loc='upper right',fontsize=25,frameon=0,handletextpad=0.3,handlelength=1.0)
  
  if mode==0:
      sm   = py.cm.ScalarMappable(cmap=cmap)
      sm.set_array([])
      cax = fig.add_axes([0.67,0.87,0.25,0.03])
      cax.tick_params(axis='both',which='both',labelsize=20,direction='in')
      cax.xaxis.set_label_coords(0.65,-2.5)
      cbar = py.colorbar(sm,cax=cax,orientation='horizontal',ticks=[0.2,0.4,0.6,0.8])
      cbar.set_label(r'\boldmath${\rm scaled}~\chi^2_{\rm red}$',size=30)
 
  py.tight_layout()

  filename = '%s/gallery/ppdfs-s-asymmetry-Q2=%3.5f'%(wdir,Q2)
  if mode==1: filename+='-bands'

  filename+='.png'

  checkdir('%s/gallery'%wdir)
  py.savefig(filename)
  print ('Saving figure to %s'%filename)
 
def plot_xf(wdir,Q2=None,mode=0,sets=False,reweight=None):

  #--get PPDF sets for comparison
  SETS = {}
  if sets:
      SETS['NNPDF'] = (QPDCALC('NNPDFpol11_100',ismc=True),'gold'    , 0.5 ,r'\textrm{\textbf{NNPDFpol1.1}}')
      SETS['JAM17'] = (QPDCALC('JAM17_PPDF_nlo',ismc=True),'darkcyan', 0.5 ,r'\textrm{\textbf{JAM17}}')
      SETS['DSSV']  = (DSSV(),                             'darkblue', 0.6 ,r'\textrm{\textbf{DSSV08}}')

  plot_xf_main    (wdir,Q2=Q2,mode=mode,SETS=SETS) 
  plot_asymmetry  (wdir,Q2=Q2,mode=mode,SETS=SETS,reweight=reweight) 
  plot_lowx       (wdir,Q2=Q2,mode=mode,SETS=SETS) 
  plot_s_asymmetry(wdir,Q2=Q2,mode=mode,SETS=SETS) 




#--ppdf moments (truncated)
def gen_moments_trunc_old(wdir, Q2 = 4, flavors = ['u','d','ub','db','s','sb','g'], mom = 1, xmin = 1e-2, xmax = 0.98):
    ## get truncated second moment integrated from a range of x_min to 1
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        if 'ppdf' not in conf['steps'][istep]['passive distributions']:
            print('ppdf-proton not an active or passive distribution')
            return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep]
    ## make sure 'parman' uses the same order for active distributions as all the replicas do

    ppdf = conf['ppdf']

    ## setup kinematics
    xs = np.geomspace(xmin,0.1,100)
    xs = np.append(xs, np.linspace(0.1, xmax, 100))
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating moment %d for ppdf-proton from %s at Q2 = %3.2f from %3.6f to %3.6f' % (mom, wdir, Q2, xmin, xmax))

    ## compute moments for all replicas
    moments = {}
    n_replicas = len(replicas)

    power = mom - 2

    for i in range(n_replicas): ## using 'scipy.integrate.cumtrapz' takes about 9.984 seconds for 100 x points, 4 flavors and 516 replicas
        lprint('%d/%d' % (i + 1, n_replicas))

        parman.order = copy.copy(replicas[i]['order'][istep])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        for flavor in flavors:
            if flavor not in moments:
                moments[flavor] = []
            if flavor == 'quark':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'u') + ppdf.get_xF(x, Q2, 'ub') + ppdf.get_xF(x, Q2, 'd') + ppdf.get_xF(x, Q2, 'db') + \
                                            ppdf.get_xF(x, Q2, 's') + ppdf.get_xF(x, Q2, 'sb') + ppdf.get_xF(x, Q2, 'c') + ppdf.get_xF(x, Q2, 'cb'))
            elif flavor == 'db-ub':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'db') - ppdf.get_xF(x, Q2, 'ub'))
            elif flavor == 'sp':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 's')  + ppdf.get_xF(x, Q2, 'sb'))
            else:
                func = lambda x:  x**power* ppdf.get_xF(x, Q2, flavor)

            function_values = [func(_) for _ in xs]
            moment_temp = cumtrapz(function_values, xs, initial = 0.0)
            moment_temp = np.array(moment_temp)
            moment_max = moment_temp[-1]
            moments[flavor].append(moment_max - moment_temp)
    print()

    checkdir('%s/data' % wdir)
    filename = '%s/data/ppdf-moment-trunc-%d-Q2=%3.5f-xmin=%3.5f-xmax=%3.5f.dat'% (wdir,mom,Q2,xmin,xmax)
    save({'X': xs, 'Q2': Q2, 'moments': moments, 'xmin': xmin, 'xmax': xmax}, filename)
    print('Saving truncated moments to %s'%filename)

def gen_moments_trunc(wdir, Q2 = 4, flavors = ['u','d','ub','db','sp','g'], mom = 1, xmin = 1e-2, xmax = 0.98):
    ## get truncated second moment integrated from a range of x_min to 1
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0])
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        if 'ppdf' not in conf['steps'][istep]['passive distributions']:
            print('ppdf-proton not an active or passive distribution')
            return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep]
    ## make sure 'parman' uses the same order for active distributions as all the replicas do

    ppdf = conf['ppdf']

    ## setup kinematics
    xs = np.geomspace(xmin,0.1,100)
    xs = np.append(xs, np.linspace(0.1, xmax, 100))
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating moment %d for ppdf-proton from %s at Q2 = %3.2f from %3.6f to %3.6f' % (mom, wdir, Q2, xmin, xmax))

    ## compute moments for all replicas
    moments = {}
    n_replicas = len(replicas)

    power = mom - 2

    for i in range(n_replicas): ## using 'scipy.integrate.cumtrapz' takes about 9.984 seconds for 100 x points, 4 flavors and 516 replicas
        lprint('%d/%d' % (i + 1, n_replicas))

        parman.order = copy.copy(replicas[i]['order'][istep])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        for flavor in flavors:
            if flavor not in moments:
                moments[flavor] = []
            if flavor == 'quark':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'u') + ppdf.get_xF(x, Q2, 'ub') + ppdf.get_xF(x, Q2, 'd') + ppdf.get_xF(x, Q2, 'db') + \
                                            ppdf.get_xF(x, Q2, 's') + ppdf.get_xF(x, Q2, 'sb') + ppdf.get_xF(x, Q2, 'c') + ppdf.get_xF(x, Q2, 'cb'))
            elif flavor == 'db-ub':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'db') - ppdf.get_xF(x, Q2, 'ub'))
            elif flavor == 'sp':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 's')  + ppdf.get_xF(x, Q2, 'sb'))
            elif flavor == 'gA':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'up')  - ppdf.get_xF(x, Q2, 'dp'))
            elif flavor == 'a8':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'up')  + ppdf.get_xF(x, Q2, 'dp') - 2*(ppdf.get_xF(x, Q2, 'sp')))
            elif flavor == 'Sigma':
                func = lambda x:  x**power*(ppdf.get_xF(x, Q2, 'up')  + ppdf.get_xF(x, Q2, 'dp') + ppdf.get_xF(x, Q2, 'sp'))
            else:
                func = lambda x:  x**power* ppdf.get_xF(x, Q2, flavor)

            function_values = [func(_) for _ in xs]
            moment_temp = cumulative_trapezoid(function_values, xs, initial = 0.0)
            moment_temp = np.array(moment_temp)
            moment_max = moment_temp[-1]
            moments[flavor].append(moment_max - moment_temp)

    checkdir('%s/data' % wdir)
    save({'X': xs, 'Q2': Q2, 'moments': moments}, '%s/data/ppdf-moment-trunc-%d-Q2=%3.5f-xmin=%3.5f-xmax=%3.5f.dat' % (wdir,mom,Q2,xmin,xmax))


#--ppdf full moments
def gen_moments(wdir, Q2 = 10, flavors = ['up','dp','u','d','ub','db','sp','g'], mom = 1):
    ## get truncated second moment integrated from a range of x_min to 1
    load_config('%s/input.py' % wdir)
    istep = core.get_istep()

    replicas = core.get_replicas(wdir)
    core.mod_conf(istep,replicas[0])
    ## 'conf' will be modified for each replica individually later in the loop over 'replicas'
    ## the reason for doing this is that 'fix parameters' has to be set correctly for each replica

    if 'ppdf' not in conf['steps'][istep]['active distributions']:
        if 'ppdf' not in conf['steps'][istep]['passive distributions']:
            print('ppdf-proton not an active or passive distribution')
            return

    resman = RESMAN(nworkers = 1, parallel = False, datasets = False)
    parman = resman.parman
    parman.order = replicas[0]['order'][istep]
    ## make sure 'parman' uses the same order for active distributions as all the replicas do

    ppdf = conf['ppdf-mom']

    ## setup kinematics
    if Q2 == None: Q2 = conf['Q20']
    print('\ngenerating full moment %d for ppdf-proton from %s at Q2 = %3.2f' % (mom, wdir, Q2))

    ## compute moments for all replicas
    moments = {}
    n_replicas = len(replicas)

    for i in range(n_replicas):
        lprint('%d/%d' % (i + 1, n_replicas))

        parman.order = copy.copy(replicas[i]['order'][istep])
        parman.set_new_params(replicas[i]['params'][istep], initial = True)

        ppdf.evolve(Q2)
        storage = ppdf.storage[Q2]

        for flavor in flavors:
            if flavor not in moments:
                moments[flavor] = []
            if flavor == 'quark':
                func = storage['u'] + storage['ub'] + storage['d'] + storage['db'] + \
                       storage['s'] + storage['sb'] + storage['c'] + storage['cb']
            elif flavor == 'db-ub':
                func = storage['db'] - storage['ub']
            elif flavor == 'sp':
                func = storage['s']  + storage['sb']
            elif flavor == 'gA':
                func = storage['u'] + storage['ub'] - storage['d'] - storage['db']
            elif flavor == 'a8':
                func = storage['u'] + storage['ub'] + storage['d'] + storage['db'] - 2*storage['s'] - 2*storage['sb']
            elif flavor == 'Sigma':
                func = storage['u'] + storage['ub'] + storage['d'] + storage['db'] + storage['s'] + storage['sb']
            else:
                func = storage[flavor]

            func = np.real(func[mom-1])

            moments[flavor].append(func)

    print()

    checkdir('%s/data' % wdir)
    save({'Q2': Q2, 'moments': moments}, '%s/data/ppdf-moment-%d-Q2=%3.5f.dat' % (wdir,mom,Q2))


