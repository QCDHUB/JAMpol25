#!/usr/bin/env python
import sys,os
import numpy as np
import copy
from subprocess import Popen, PIPE, STDOUT
import pandas as pd

path='/w/jam-sciwork24/ccocuzza/legacy'
sys.path.insert(0,path)
os.environ['FITPACK']=path

#--matplotlib
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import MultipleLocator
import pylab as py

#--from tools
from tools.tools     import load,save,checkdir,lprint,isnumeric
from tools.config    import conf,load_config

#--from fitlib
from fitlib.resman import RESMAN

#--from local
from analysis.corelib import core
from analysis.corelib import classifier

import scipy.stats
import math
import scipy
from sympy import *

from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('Agg')
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Times New Roman"
})
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors
import pylab as py

regen = True
#regen = False

ext = '.png'
#ext = '.pdf'

wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/HT4'

def get_res(wdir):

    replicas=core.get_replicas(wdir)
    load_config(wdir + '/input.py')
    istep=core.get_istep()
    core.mod_conf(istep,replicas[0])
    conf['predict'] = True
    conf['bootstrap']=False
    resman=RESMAN(nworkers=1,parallel=False,datasets=False)
    parman=resman.parman
    order=parman.order
  
    residuals = resman.lqcd_itd_res

    value = residuals.data['value']
    cov   = residuals.cov
    dim   = residuals.dim

    inv_cov=np.linalg.inv(cov)
    eigenvalues, eigenvectors = np.linalg.eig(inv_cov)
    sqrt_eigenvalues=np.sqrt(eigenvalues)
    II=np.argsort(sqrt_eigenvalues)
    diag=np.diag(eigenvalues)

    nrep = len(replicas) 
    cnt =0
    R=[]
    res = []
    for replica in replicas:
        lprint('Processing: [%s/%s]'%(cnt+1,nrep))
        parman.set_new_params(replica['params'][istep],initial=True)
        res_lattice,exp,thy = residuals.get_residuals()
        res.append(res_lattice)
        r=value[0]-thy
        r=np.einsum('ij,j',eigenvectors.T,r)
        r*=sqrt_eigenvalues
        R.append(r)
        cnt +=1

    res = np.array(res)

    chi2 = np.mean(np.sum(res**2,axis=1))/dim

    R=np.abs(np.mean(R,axis=0)[II])

    filename = './data/lqcd_%s.npy'%wdir.split('/')[-1]
    np.save(filename,R)
    print('Saving to %s'%filename)
    #return chi2, R

if __name__=="__main__":

    wdir_pos = wdir+'/pos_g'
    wdir_neg = wdir+'/neg_g'

    checkdir('data')
    filename_pos = './data/lqcd_pos_g.npy'
    filename_neg = './data/lqcd_neg_g.npy'
    if regen:
        get_res(wdir_pos)
        get_res(wdir_neg)
    try:
        R_pos = np.load(filename_pos)
        R_neg = np.load(filename_neg)
    except:
        get_res(wdir_pos)
        get_res(wdir_neg)
        R_pos = np.load(filename_pos)
        R_neg = np.load(filename_neg)

    ncols = 1
    nrows = 1
    fig,axs = py.subplots(nrows,ncols,figsize=(ncols*8,nrows*6))
    
    axs.scatter(np.arange(0,len(R_pos)),R_pos**2,color = 'r')
    axs.scatter(np.arange(0,len(R_neg)),R_neg**2,color = 'b')

    axs.set_xlim(-1,50)
    axs.set_xticks([0,10,20,30,40])
    axs.xaxis.set_minor_locator(MultipleLocator(5))
    axs.set_ylim(-1,22)
    #axs.set_yticks([0,10,20])
    axs.set_yticks([0,5,10,15,20])
    axs.yaxis.set_minor_locator(MultipleLocator(1))
    axs.set_xlabel(r'\rm \bf eigendirection',fontsize = 35)
    axs.tick_params(axis='both', which='major', top=True, right = True, direction='in',labelsize=30,length=8)
    axs.tick_params(axis='both', which='minor', top=True, right = True, direction='in',labelsize=30,length=4)
    
    # axs[0].set_ylabel(r'\rm \boldmath{$|r_{\Delta g < 0}^*| - |r_{\Delta g > 0}^*|$}',fontsize = 30)
    axs.set_ylabel(r'\rm \boldmath{$\chi^2_i$}',fontsize = 35)
    
    #axs[0].set_title(r'\boldmath{$\rm baseline + LQCD$}',size=35)
    #axs[1].set_title(r'\rm \bf + high-\boldmath{$x$} DIS',size=35)
    
    axs.axhline(1,ls = '--',color = 'k')
    axs.text(0.65,0.89,r'\boldmath{$\Delta g > 0$}',transform=axs.transAxes,size=35,color = 'r')
    axs.text(0.65,0.75,r'\boldmath{$\Delta g < 0$}',transform=axs.transAxes,size=35,color = 'b')
    # axs[1].text(0.4,0.78,r'\boldmath{$\chi^2_{dof}(\Delta g > 0) = 0.58$}',transform=axs[1].transAxes,size=20,color = 'r')
    # axs[1].text(0.4,0.9,r'\boldmath{$\chi^2_{dof}(\Delta g < 0) = 3.92$}',transform=axs[1].transAxes,size=20,color = 'b')
    
    
    checkdir('plots')
    filename = 'plots/fig_lattice'
    filename += ext
    py.tight_layout()
    py.subplots_adjust(top=0.99,right=0.99,hspace=0.01,wspace=0.01)
    
    axs.set_rasterized(True)
    
    py.savefig(filename)
    py.close()
    print('Saving figure to %s'%filename)






