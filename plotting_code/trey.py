#!/usr/bin/env python
import sys, os
import numpy as np
import copy
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

import scipy

#--matplotlib
import matplotlib
matplotlib.use('Agg')
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm
import pylab as py

## from fitpack tools
from tools.tools     import load, save, checkdir, lprint
from tools.config    import conf, load_config

## from fitpack fitlib
from fitlib.resman import RESMAN

## from fitpack analysis
from analysis.corelib import core, classifier, jar, summary
from analysis.qpdlib import ppdf as _ppdf

import kmeanconf as kc

import argparse

wdir = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/results/simul2/pos_g'

cwd = 'plotting_code'


#--force functions to be generated again
regen = {}
#regen['chris'] = True
regen['chris'] = False

regen['trey'] = False

#--generation codes
def gen_predictions(params,resman,filename):


    parman=resman.parman
    order=parman.order
    datasets = copy.copy(conf['datasets'])
       
    obsres={}
    if 'idis'     in conf['datasets'] : obsres['idis']     = resman.idis_res
    if 'pidis'    in conf['datasets'] : obsres['pidis']    = resman.pidis_res
    if 'SU23'     in conf['datasets'] : obsres['SU23']     = resman.SU23_res
    if 'sidis'    in conf['datasets'] : obsres['sidis']    = resman.sidis_res
    if 'psidis'   in conf['datasets'] : obsres['psidis']   = resman.psidis_res
    if 'dy'       in conf['datasets'] : obsres['dy']       = resman.dy_res
    if 'wzrv'     in conf['datasets'] : obsres['wzrv']     = resman.wzrv_res
    if 'wasym'    in conf['datasets'] : obsres['wasym']    = resman.wasym_res
    if 'zrap'     in conf['datasets'] : obsres['zrap']     = resman.zrap_res
    if 'sia'      in conf['datasets'] : obsres['sia']      = resman.sia_res
    if 'qpdf'     in conf['datasets'] : obsres['qpdf']     = resman.qpdf_res
    if 'dy-pion'  in conf['datasets'] : obsres['dy-pion']  = resman.dy_pion_res
    if 'pion_qT'  in conf['datasets'] : obsres['pion_qT']  = resman.pion_qTres
    if 'ln'       in conf['datasets'] : obsres['ln']       = resman.ln_res
    if 'jet'      in conf['datasets'] : obsres['jet']      = resman.jet_res
    if 'pjet'     in conf['datasets'] : obsres['pjet']     = resman.pjet_res
    if 'wc'       in conf['datasets'] : obsres['wc']       = resman.wc_res 
    #--setup big table to store all we want
    data={}
    data['order']=order
    data['params']=params
    data['reactions']={}
    data['res']=[]
    data['rres']=[]
    data['nres']=[]

    for _ in obsres:
        tabs=copy.copy(obsres[_].tabs)
        #--create a space to store all the predictions from replicas
        for idx in tabs:
            tabs[idx]['prediction-rep']=[]
            tabs[idx]['residuals-rep']=[]
            tabs[idx]['shift-rep']=[]
            tabs[idx]['rres-rep']=[]
        data['reactions'][_]=tabs

    print('generating predictions for %s'%filename)

    #--total replica count
    nrep = params.shape[0]

    for i in range(nrep):
        lprint('progress: %d/%d'%(i+1,nrep))

        #--compute residuals (==theory)
        res,rres,nres=resman.get_residuals(params[i],initial=True)
        data['res'].append(res)
        data['rres'].append(rres)
        data['nres'].append(nres)

        #--save predictions of the current step and current replica at data
        for _ in obsres:
            for idx in data['reactions'][_]:
                prediction= copy.copy(obsres[_].tabs[idx]['prediction'])
                residuals = copy.copy(obsres[_].tabs[idx]['residuals'])
                shift     = copy.copy(obsres[_].tabs[idx]['shift'])
                rres      = obsres[_]._get_rres(idx)
                data['reactions'][_][idx]['prediction-rep'].append(prediction)
                data['reactions'][_][idx]['residuals-rep'].append(residuals)
                data['reactions'][_][idx]['shift-rep'].append(shift)
                data['reactions'][_][idx]['rres-rep'].append(rres)


    #--convert tables to numpy array before saving
    for _ in ['res','rres','nres']:
        data[_]=np.array(data[_])

    np.save(filename,data)
    print 


def gen_pdf(params,filename,Q2=10):
    nrep = params.shape[0]
    #--generate pdfs
    X = np.geomspace(1e-6,0.1,200)
    X = np.append(X,np.linspace(0.1,1))
    pdf = conf['pdf']

    flavs = ['u','d','g','uv','dv','ub','db','s','sb']
    data = {}
    data['X']  = X
    data['Q2'] = Q2
    data['XF'] = {}
    for flav in flavs: data['XF'][flav] = []
    for i in range(nrep):
        lprint('generating unpolarized PDFs [%s/%s]'%(i+1,nrep))
        parman.set_new_params(params[i],initial=True)
        for flav in flavs: data['XF'][flav].append([pdf.get_xF(x,Q2,flav) for x in X])

    for flav in flavs: data['XF'][flav] = np.array(data['XF'][flav])

    np.save(filename,data)
    print('Saving unpolarized PDF data to %s'%filename)

def gen_ff(had,params,filename,Q2=100):
    nrep = params.shape[0]
    #--generate pdfs
    Z = np.geomspace(1e-6,0.1,200)
    Z = np.append(Z,np.linspace(0.1,1))
    ff = conf['ff%s'%had]

    flavs = ['u','d','s','c','b','ub','db','sb','cb','bb','g']
    data = {}
    data['Z']  = Z
    data['Q2'] = Q2
    data['ZF'] = {}
    for flav in flavs: data['ZF'][flav] = []
    for i in range(nrep):
        lprint('generating FF %s [%s/%s]'%(had,i+1,nrep))
        parman.set_new_params(params[i],initial=True)
        for flav in flavs: data['ZF'][flav].append([ff.get_xF(z,Q2,flav) for z in Z])

    for flav in flavs: data['ZF'][flav] = np.array(data['ZF'][flav])

    np.save(filename,data)
    print('Saving FF %s data to %s'%(had,filename))


#--summary codes
def get_z_score(chi2_red,dof):
    x = chi2_red
    cdf = scipy.stats.chi2.cdf(x,dof)
    p_val = 1-cdf
    return -scipy.stats.norm.ppf(p_val)

def get_norm(params,order): 
    tab={}
    for par in params: 
        for i in range(len(order)):
            if order[i][0]==2:
                reaction=order[i][1]
                idx=order[i][2]
                if reaction not in tab: tab[reaction]={}
                if idx not in tab[reaction]: tab[reaction][idx]=[] 
                tab[reaction][idx].append(par[i])

        for k in conf['datasets']:
            for kk in conf['datasets'][k]['norm']:
                if conf['datasets'][k]['norm'][kk]['fixed'] == True:  continue
                if conf['datasets'][k]['norm'][kk]['fixed'] == False: continue
                reference_norm = conf['datasets'][k]['norm'][kk]['fixed']
                if k  not in tab: tab[k]={}
                if kk not in tab[k]: tab[k][kk]=[] 
                tab[k][kk].append(tab[k][reference_norm])
                       
    for reaction in tab:
        for idx in tab[reaction]:
            norm=tab[reaction][idx][:]
            tab[reaction][idx]={}
            tab[reaction][idx]['mean']=np.mean(norm)
            tab[reaction][idx]['std']=np.std(norm)
    return tab

def get_chi2(predictions): 

    data=predictions['reactions']
    tab={}
    for reaction in data:
        if len(list(data[reaction])) == 0: continue
        if reaction not in tab: 
            tab[reaction]={}
            tab[reaction]['chi2'] = 0
            tab[reaction]['npts'] = 0
        for idx in data[reaction]:
            if idx not in tab[reaction]: tab[reaction][idx]={}
      
            value=data[reaction][idx]['value']
            alpha=data[reaction][idx]['alpha']
            if 'rres-rep' in data[reaction][idx].keys():
                rres=np.mean(data[reaction][idx]['rres-rep'],axis=0)
            else:
                rres=0.
            if np.isnan(rres).any(): rres=0.0
            thy=np.mean(data[reaction][idx]['prediction-rep'],axis=0)
            col=data[reaction][idx]['col'][0]
            if 'target' in            data[reaction][idx]: tar=data[reaction][idx]['target'][0]
            elif 'tar' in             data[reaction][idx]: tar=data[reaction][idx]['tar'][0]
            elif 'particles-in' in    data[reaction][idx]: tar=data[reaction][idx]['particles-in'][0]
            elif 'reaction' in        data[reaction][idx]: tar=data[reaction][idx]['reaction'][0]
            else:                                          tar='-' 
            if tar=='deuteron': tar='d'
            if 'hadron' in            data[reaction][idx]: had=data[reaction][idx]['hadron'][0]
            else:                                          had='-'
            obs=data[reaction][idx]['obs'][0]
            #--rewrite some observables to be more compact
            if reaction=='jet':
                obs = obs.replace('<','').replace('>','').replace('_over_','/').replace('d2_sigma','dsig')\
                         .replace('d_y','dy').replace('d_pt','dpt').replace('2_pi','2pi').replace('_',' ')
            if reaction=='pjet':
                obs = obs.replace('<','').replace('>','')
            res=(value-thy)/alpha
            chi2=np.sum(res**2)#+np.sum(rres**2)
            npts=res.size
            chi2_npts=chi2/npts
            zscore = get_z_score(chi2,npts)

            tab[reaction][idx]['col']       = col
            tab[reaction][idx]['tar']       = tar
            tab[reaction][idx]['had']       = had
            tab[reaction][idx]['obs']       = obs
            tab[reaction][idx]['chi2']      = chi2
            tab[reaction][idx]['npts']      = npts
            tab[reaction][idx]['chi2_npts'] = chi2_npts
            tab[reaction][idx]['zscore']    = zscore
            tab[reaction]['chi2']          += chi2
            tab[reaction]['npts']          += npts
        tab[reaction]['chi2_npts'] = tab[reaction]['chi2']/tab[reaction]['npts']
        tab[reaction]['zscore']    = get_z_score(tab[reaction]['chi2'],tab[reaction]['npts'])

    return tab

def print_summary(DATA):
    params, order, predictions = DATA['params'], DATA['order'], DATA['predictions'] 
    norm_tab=get_norm(params,order)
    chi2_tab=get_chi2(predictions)

    L = params.shape[0] 

    #--adjust width of col and obs columns
    lobs,lcol = [],[]
    for reaction in chi2_tab:
        for idx in chi2_tab[reaction]:
            if idx in ['chi2','npts','chi2_npts','zscore']: continue
            lcol.append(len(chi2_tab[reaction][idx]['col']))
            lobs.append(len(chi2_tab[reaction][idx]['obs']))
    lcol = np.max(lcol)-2
    lobs = np.max(lobs)-2

    msg1='%10s '
    msg1+='%15s '
    msg1+=' '*lcol
    msg1+='%s '
    msg1+='%10s '
    msg1+='%10s '
    msg1+=' '*lobs
    msg1+='%s '
    msg1+='%5s '
    msg1+='%10s '
    msg1+='%10s '
    msg1+='%10s '
    msg1+='%12s '
    msg1=msg1%('reaction','idx','col','tar','had','obs','npts','chi2','chi2/npts','zscore','norm')
    print(msg1)
    chi2_tot=0
    npts_tot=0
    for reaction in chi2_tab:
        l = list(chi2_tab[reaction])
        for _ in ['chi2','npts','chi2_npts','zscore']: l.remove(_)
        for idx in sorted(l):
            col=chi2_tab[reaction][idx]['col']
            tar=chi2_tab[reaction][idx]['tar']
            had=chi2_tab[reaction][idx]['had']
            obs=chi2_tab[reaction][idx]['obs']
            npts=chi2_tab[reaction][idx]['npts']
            chi2=chi2_tab[reaction][idx]['chi2']
            chi2_npts=chi2_tab[reaction][idx]['chi2_npts']
            zscore=chi2_tab[reaction][idx]['zscore']
            if reaction in norm_tab:
                if idx in norm_tab[reaction]:
                    mean=norm_tab[reaction][idx]['mean']
                    std =norm_tab[reaction][idx]['std']
                    std =str(int(np.round(std*1000)))
                    norm = '%4.3f'%mean
                    if len(std)==1: norm+=' '
                    norm += '('+std+')'
                else:
                    norm='N/A'
            else: norm = 'N/A'

            chi2_tot+=chi2
            npts_tot+=npts

            msg2 ='%10s '
            msg2+='%15s '
            msg2+=' '*(lcol-len(col)+3)
            msg2+='%s '
            msg2+='%10s '
            msg2+='%10s '
            msg2+=' '*(lobs-len(obs)+3)
            msg2+='%s '
            msg2+='%5d '
            msg2+='%10.2f '
            msg2+='%10.2f '
            msg2+='%10.2f '
            if norm=='N/A': msg2 += '%12s ' 
            else:           msg2 += '%12s ' 
            print(msg2%(reaction,idx,col,tar,had,obs,npts,chi2,chi2_npts,zscore,norm))


    print("-"*len(msg1))
    #--print chi2 per experiment
    for reaction in chi2_tab:
        npts      = chi2_tab[reaction]['npts']
        chi2      = chi2_tab[reaction]['chi2']
        chi2_npts = chi2_tab[reaction]['chi2_npts']
        zscore    = chi2_tab[reaction]['zscore']

        msg3 ='%10s '
        msg3+='%15s '
        msg3+=' '*(lcol+2)
        msg3+='%s '
        msg3+='%8s '
        msg3+='%s '
        msg3+='%10s '
        msg3+=' '*(lobs+2)
        msg3+='%s '
        msg3+='%5d '
        msg3+='%10.2f '
        msg3+='%10.2f '
        msg3+='%10.2f '
        print(msg3%(reaction,' ',' ',' ',' ',' ',' ',npts,chi2,chi2_npts,zscore))

    chi2_npts_tot=chi2_tot/npts_tot
    zscore_tot = get_z_score(chi2_tot,npts_tot)

    zscore_tot = get_z_score(chi2_tot,npts_tot)
    print("-"*len(msg1))
    msg4 ='%10s '
    msg4+='%15s '
    msg4+=' '*lcol
    msg4+='%s '
    msg4+='%8s '
    msg4+='%s '
    msg4+='%10s '
    msg4+=' '*(lobs+4)
    msg4+='%s '
    msg4+='%5d '
    msg4+='%10.2f '
    msg4+='%10.3f '
    msg4+='%10.3f '
    msg4+='%12s'
    print(msg4%('total',' ',' ',' ',' ',' ',' ',npts_tot,chi2_tot,chi2_npts_tot,zscore_tot,' '))


#--plotting codes
def plot_pdf(DATA,Q2=10,mode=0):
    #--mode 0: plot each replica
    #--mode 1: plot average and standard deviation of replicas 

    nrows,ncols=3,2
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*8,nrows*4))
    axs,axLs = {},{}
    for i in range(N):
        axs[i+1] = py.subplot(nrows,ncols,i+1)
        divider = make_axes_locatable(axs[i+1])
        axLs[i+1] = divider.append_axes("right",size=3.10,pad=0,sharey=axs[i+1])
        axLs[i+1].set_xlim(0.1,0.9)
        axLs[i+1].spines['left'].set_visible(False)
        axLs[i+1].yaxis.set_ticks_position('right')
        py.setp(axLs[i+1].get_xticklabels(),visible=True)

        axs[i+1].spines['right'].set_visible(False)
    nrows,ncols=3,2

    hand = {}
    thy  = {}

    flavs = ['uv','dv','g','db+ub','db-ub','s+sb','Rs']

    for name in DATA:
        data = DATA[name]['pdf']
        load_config(DATA[name]['input'])
        X=data['X']

        XF = data['XF']

        if name=='chris': color='red'
        if name=='trey':  color='blue'
        for flav in flavs:

            if flav=='uv' or flav=='dv': ax,axL = axs[1],axLs[1]
            elif flav=='g':              ax,axL = axs[2],axLs[2]
            elif flav=='db+ub':          ax,axL = axs[3],axLs[3]
            elif flav=='db-ub':          ax,axL = axs[4],axLs[4]
            elif flav=='s+sb':           ax,axL = axs[5],axLs[5]
            elif flav=='Rs':             ax,axL = axs[6],axLs[6]

            if   flav=='db+ub': pdf = XF['ub'] + XF['db']
            elif flav=='db-ub': pdf = XF['db'] - XF['ub']
            elif flav=='s+sb' : pdf = XF['s']  + XF['sb']
            elif flav=='Rs'   : pdf = (XF['s'] + XF['sb'])/(XF['ub']  + XF['db'])
            elif flav=='g'    : pdf = XF['g']/10
            else:               pdf = XF[flav]
            mean = np.mean(pdf,axis=0)
            std  = np.std (pdf,axis=0)

            #--plot each replica
            if mode==0:
                for i in range(len(pdf)):
                    hand[name] ,= ax .plot(X,np.array(pdf[i]),color=color,alpha=0.10)
                    hand[name] ,= axL.plot(X,np.array(pdf[i]),color=color,alpha=0.10)
 
            #--plot average and standard deviation
            if mode==1:

                where = [1 for i in range(len(X))]
                if flav=='Rs':
                    where = []
                    for x in X:
                        if x < 0.2: where.append(1)
                        if x > 0.2: where.append(0)

                hand[name] = ax .fill_between(X,(mean-std),(mean+std),color=color,alpha=0.7,zorder=5,where=np.array(where))
                hand[name] = axL.fill_between(X,(mean-std),(mean+std),color=color,alpha=0.7,zorder=5,where=np.array(where))


    for i in range(N):
        axs[i+1].set_xlim(1e-4,0.1)
        axs[i+1].semilogx()

        axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
        axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
        axs[i+1].tick_params(axis='x',    which='major', pad = 8)
        axs[i+1].set_xticks([0.0001,0.001,0.01,0.1])
        axs[i+1].set_xticklabels([r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$',r'$0.1$'])

        axLs[i+1].set_xlim(0.1,1.0)

        axLs[i+1].tick_params(axis='both', which='major', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=8)
        axLs[i+1].tick_params(axis='both', which='minor', top=True, right=True, left=False, labelright=False, direction='in',labelsize=30,length=4)
        axLs[i+1].tick_params(axis='x',    which='major', pad = 8)
        axLs[i+1].set_xticks([0.3,0.5,0.7])
        axLs[i+1].set_xticklabels([r'$0.3$',r'$0.5$',r'$0.7$'])

    for i in [1,2,3,4]:
        axs[i] .tick_params(labelbottom=False)
        axLs[i].tick_params(labelbottom=False)

    axs[1].set_ylim(0,0.7)
    axs[2].set_ylim(0,1.2)
    axs[3].set_ylim(-0.05,1.5)
    axs[4].set_ylim(-0.09,0.09)
    axs[5].set_ylim(0,1.5)
    axs[6].set_ylim(0,1.2)

    axs[1].set_yticks([0.2,0.4,0.6])
    axs[2].set_yticks([0.2,0.4,0.6,0.8,1.0])
    axs[3].set_yticks([0,0.2,0.4,0.6,0.8,1.0,1.2,1.4])
    axs[4].set_yticks([-0.08,-0.04,0,0.04,0.08])
    axs[5].set_yticks([0.2,0.4,0.6,0.8,1.0,1.2,1.4])
    axs[6].set_yticks([0.5,1.0])

    for i in range(N):
        axs [i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)
        axLs[i+1].axvline(0.1,color='k',linestyle=':' ,alpha=0.5)

    axs [3].axhline(0.0,ls='--',color='black',alpha=0.5)
    axs [4].axhline(0.0,ls='--',color='black',alpha=0.5)
    axLs[3].axhline(0.0,ls='--',color='black',alpha=0.5)
    axLs[4].axhline(0.0,ls='--',color='black',alpha=0.5)

    axLs[5].set_xlabel(r'\boldmath$x$',size=40)
    axLs[6].set_xlabel(r'\boldmath$x$',size=40)   
    axLs[5].xaxis.set_label_coords(0.92,0.00)
    axLs[6].xaxis.set_label_coords(0.92,0.00)

    axLs[1].text(0.50 ,0.50  ,r'\boldmath{$xu_{v}$}'            , transform=axLs[1].transAxes,size=30)
    axLs[1].text(0.01 ,0.20  ,r'\boldmath{$xd_{v}$}'            , transform=axLs[1].transAxes,size=30)
    axLs[2].text(0.05 ,0.25  ,r'\boldmath{$xg/10$}'             , transform=axLs[2].transAxes,size=30)
    axs[3] .text(0.10 ,0.10  ,r'\boldmath{$x(\bar{d}+\bar{u})$}', transform=axs[3] .transAxes,size=30)
    axLs[4].text(0.20 ,0.10  ,r'\boldmath{$x(\bar{d}-\bar{u})$}', transform=axLs[4].transAxes,size=30)
    axs[5] .text(0.10 ,0.10  ,r'\boldmath{$x(s+\bar{s})$}',       transform=axs[5] .transAxes,size=30)
    axs[6] .text(0.10 ,0.10  ,r'\boldmath{$R_s$}',                transform=axs[6] .transAxes,size=30)

    if Q2 == 1.27**2: axs[2].text(0.05,0.08,r'$Q^2 = m_c^2$'                                  , transform=axs[2].transAxes,size=30)
    else:             axs[2].text(0.05,0.08,r'$Q^2 = %s$'%Q2 + ' ' + r'\textrm{GeV}' + r'$^2$', transform=axs[2].transAxes,size=30)

    minorLocator = MultipleLocator(0.05)
    axs[1].yaxis.set_minor_locator(minorLocator)
    axs[2].yaxis.set_minor_locator(minorLocator)
    axs[3].yaxis.set_minor_locator(minorLocator)
    axs[5].yaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(0.005)
    axs[4].yaxis.set_minor_locator(minorLocator)
    minorLocator = MultipleLocator(0.1)
    axs[6].yaxis.set_minor_locator(minorLocator)

    handles,labels = [],[]
    handles.append(hand['chris'])
    handles.append(hand['trey'])
    labels.append(r'\textrm{\textbf{Chris}}')
    labels.append(r'\textrm{\textbf{Trey}}')
    axs[1].legend(handles,labels,loc='upper left', fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.subplots_adjust(hspace = 0, wspace = 0.20)

    filename = '%s/gallery_trey/pdfs-Q2=%3.5f'%(cwd,pdf_Q2)
    if mode==1: filename += '-bands'

    filename+='.png'

    py.savefig(filename)
    py.clf()
    print ('Saving figure to %s'%filename)

def plot_ffpion(DATA,Q2=100,mode=0):
   
    nrows,ncols=3,3
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*9,nrows*5))
    axs = {}
    for i in range(N):
        axs[i+1] = py.subplot(nrows,ncols,i+1)

    hand = {}
    j = 0

    flavs = ['u','d','s','ub','db','sb','c','b','g']

    for name in DATA:

        data = DATA[name]['ffpion']
        load_config(DATA[name]['input'])
        Z=data['Z']

        ZF = data['ZF']

        if name=='chris': color='red'
        if name=='trey':  color='blue'

        for flav in flavs:

            if   flav=='u':  ax = axs[1]
            elif flav=='d':  ax = axs[2]
            elif flav=='s':  ax = axs[3]
            elif flav=='ub': ax = axs[4]
            elif flav=='db': ax = axs[5]
            elif flav=='sb': ax = axs[6]
            elif flav=='c':  ax = axs[7]
            elif flav=='b':  ax = axs[8]
            elif flav=='g':  ax = axs[9]

            ff = ZF[flav]

            mean = np.mean(ff,axis=0)
            std  = np.std (ff,axis=0)

            #--plot each replica
            if mode==0:
                for i in range(len(ff)):
                    hand[name] ,= ax.plot(Z,ff[i],color=color,alpha=0.10)
 
            #--plot average and standard deviation
            if mode==1:

                where = [1 for i in range(len(Z))]

                hand[name]  = ax.fill_between(Z,(mean-std),(mean+std),color=color,alpha=0.7,zorder=5,where=np.array(where))

   
    for i in range(N):
          axs[i+1].set_xlim(0.2,0.9)

          axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
          axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
          axs[i+1].set_xticks([0.4,0.6,0.8])

          axs[i+1].axhline(0,0,1,color='black')
          minorLocator = MultipleLocator(0.05)
          axs[i+1].xaxis.set_minor_locator(minorLocator)

    axs [1].tick_params(labelbottom=False)
    axs [2].tick_params(labelbottom=False)
    axs [3].tick_params(labelbottom=False)
    axs [4].tick_params(labelbottom=False)
    axs [5].tick_params(labelbottom=False)
    axs [6].tick_params(labelbottom=False)


    axs [2].tick_params(labelleft=False)
    axs [3].tick_params(labelleft=False)
    axs [5].tick_params(labelleft=False)
    axs [6].tick_params(labelleft=False)
    axs [8].tick_params(labelleft=False)
    axs [9].tick_params(labelleft=False)


    for i in [1,2,3]:
        axs[i].set_ylim(-0.2,0.70)

        axs[i].set_yticks([0.2,0.4,0.6])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)
   
    for i in [4,5,6]:
        axs[i].set_ylim(-0.2,0.70)

        axs[i].set_yticks([0.2,0.4,0.6])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)

    for i in [7,8,9]:
        axs[i].set_ylim(-0.2,0.70)

        axs[i].set_yticks([0.2,0.4,0.6])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)

    axs[7].set_xlabel(r'\boldmath$z$',size=40)   
    axs[7].xaxis.set_label_coords(0.95,0.00)
    axs[8].set_xlabel(r'\boldmath$z$',size=40)   
    axs[8].xaxis.set_label_coords(0.95,0.00)
    axs[9].set_xlabel(r'\boldmath$z$',size=40)   
    axs[9].xaxis.set_label_coords(0.95,0.00)

    axs[1].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{u}$'       , transform=axs[1].transAxes,size=40)
    axs[2].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{d}$'       , transform=axs[2].transAxes,size=40)
    axs[3].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{s}$'       , transform=axs[3].transAxes,size=40)
    axs[4].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{\bar{u}}$' , transform=axs[4].transAxes,size=40)
    axs[5].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{\bar{d}}$' , transform=axs[5].transAxes,size=40)
    axs[6].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{\bar{s}}$' , transform=axs[6].transAxes,size=40)
    axs[7].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{c}$'       , transform=axs[7].transAxes,size=40)
    axs[8].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{b}$'       , transform=axs[8].transAxes,size=40)
    axs[9].text(0.80,0.83,r'\boldmath$z D^{\pi^+}_{g}$'       , transform=axs[9].transAxes,size=40)

    axs[4].text(0.07,0.85,r'$Q^2 = %d$~'%Q2  + r'\textrm{GeV}'+r'$^2$', transform=axs[4].transAxes,size=25)


    blank ,= axs[1].plot(0,0,color='white',alpha=0.0)

    handles,labels = [],[]
    handles.append(hand['chris'])
    handles.append(hand['trey'])
    labels.append(r'\textrm{\textbf{Chris}}')
    labels.append(r'\textrm{\textbf{Trey}}')
    axs[3].legend(handles,labels,loc='upper left', fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.subplots_adjust(hspace=0.02,wspace=0.02,top=0.99,right=0.99,left=0.03)

    filename = '%s/gallery_trey/ffs-pion-Q2=%3.5f'%(cwd,Q2)
    if mode==1: filename+='-bands'
    filename+='.png'
    for i in range(N):
        axs[i+1] .set_rasterized(True)

    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def plot_ffkaon(DATA,Q2=100,mode=0):
   
    wdir1 = 'results/star/final'

    WDIR = [wdir1]

    nrows,ncols=3,3
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*9,nrows*5))
    axs = {}
    for i in range(N):
        axs[i+1] = py.subplot(nrows,ncols,i+1)

    hand = {}
    j = 0

    Q2 = 100

    flavs = ['u','d','s','ub','db','sb','c','b','g']

    for name in DATA:

        data = DATA[name]['ffkaon']
        load_config(DATA[name]['input'])
        Z=data['Z']

        ZF = data['ZF']

        if name=='chris': color='red'
        if name=='trey':  color='blue'

        for flav in flavs:

            if   flav=='u':  ax = axs[1]
            elif flav=='d':  ax = axs[2]
            elif flav=='s':  ax = axs[3]
            elif flav=='ub': ax = axs[4]
            elif flav=='db': ax = axs[5]
            elif flav=='sb': ax = axs[6]
            elif flav=='c':  ax = axs[7]
            elif flav=='b':  ax = axs[8]
            elif flav=='g':  ax = axs[9]

            ff = ZF[flav]

            mean = np.mean(ff,axis=0)
            std  = np.std (ff,axis=0)

            #--plot each replica
            if mode==0:
                for i in range(len(ff)):
                    hand[name] ,= ax.plot(Z,ff[i],color=color,alpha=0.10)
 
            #--plot average and standard deviation
            if mode==1:

                where = [1 for i in range(len(Z))]

                hand[name]  = ax.fill_between(Z,(mean-std),(mean+std),color=color,alpha=0.7,zorder=5,where=np.array(where))



   
    for i in range(N):
          axs[i+1].set_xlim(0.2,0.9)

          axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
          axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
          axs[i+1].set_xticks([0.4,0.6,0.8])

          minorLocator = MultipleLocator(0.05)
          axs[i+1].xaxis.set_minor_locator(minorLocator)

    axs [1].tick_params(labelbottom=False)
    axs [2].tick_params(labelbottom=False)
    axs [3].tick_params(labelbottom=False)
    axs [4].tick_params(labelbottom=False)
    axs [5].tick_params(labelbottom=False)
    axs [6].tick_params(labelbottom=False)


    axs [2].tick_params(labelleft=False)
    axs [3].tick_params(labelleft=False)
    axs [5].tick_params(labelleft=False)
    axs [6].tick_params(labelleft=False)
    axs [8].tick_params(labelleft=False)
    axs [9].tick_params(labelleft=False)


    for i in [1,2,3]:
        axs[i].set_ylim(0,0.40)

        axs[i].set_yticks([0.1,0.2,0.3])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)
   
    for i in [4,5,6]:
        axs[i].set_ylim(0,0.50)

        axs[i].set_yticks([0.1,0.2,0.3,0.4])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)

    for i in [7,8,9]:
        axs[i].set_ylim(0,0.30)

        axs[i].set_yticks([0.1,0.2])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)

    axs[7].set_xlabel(r'\boldmath$z$',size=40)   
    axs[7].xaxis.set_label_coords(0.95,0.00)
    axs[8].set_xlabel(r'\boldmath$z$',size=40)   
    axs[8].xaxis.set_label_coords(0.95,0.00)
    axs[9].set_xlabel(r'\boldmath$z$',size=40)   
    axs[9].xaxis.set_label_coords(0.95,0.00)

    axs[1].text(0.80,0.83,r'\boldmath$z D^{K^+}_{u}$'       , transform=axs[1].transAxes,size=40)
    axs[2].text(0.80,0.83,r'\boldmath$z D^{K^+}_{d}$'       , transform=axs[2].transAxes,size=40)
    axs[3].text(0.80,0.83,r'\boldmath$z D^{K^+}_{s}$'       , transform=axs[3].transAxes,size=40)
    axs[4].text(0.80,0.83,r'\boldmath$z D^{K^+}_{\bar{u}}$' , transform=axs[4].transAxes,size=40)
    axs[5].text(0.80,0.83,r'\boldmath$z D^{K^+}_{\bar{d}}$' , transform=axs[5].transAxes,size=40)
    axs[6].text(0.80,0.83,r'\boldmath$z D^{K^+}_{\bar{s}}$' , transform=axs[6].transAxes,size=40)
    axs[7].text(0.80,0.83,r'\boldmath$z D^{K^+}_{c}$'       , transform=axs[7].transAxes,size=40)
    axs[8].text(0.80,0.83,r'\boldmath$z D^{K^+}_{b}$'       , transform=axs[8].transAxes,size=40)
    axs[9].text(0.80,0.83,r'\boldmath$z D^{K^+}_{g}$'       , transform=axs[9].transAxes,size=40)

    axs[4].text(0.07,0.85,r'$Q^2 = %d$~'%Q2  + r'\textrm{GeV}'+r'$^2$', transform=axs[4].transAxes,size=25)

    blank ,= axs[1].plot(0,0,color='white',alpha=0.0)

    handles,labels = [],[]
    handles.append(hand['chris'])
    handles.append(hand['trey'])
    labels.append(r'\textrm{\textbf{Chris}}')
    labels.append(r'\textrm{\textbf{Trey}}')
    axs[1].legend(handles,labels,loc='upper left', fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.subplots_adjust(hspace=0.02,wspace=0.02,top=0.99,right=0.99,left=0.03)

    filename = '%s/gallery_trey/ffs-kaon-Q2=%3.5f'%(cwd,Q2)
    if mode==1: filename+='-bands'
    filename+='.png'
    for i in range(N):
        axs[i+1] .set_rasterized(True)

    py.savefig(filename)
    print ('Saving figure to %s'%filename)

def plot_ffhadron(DATA,Q2=100,mode=0):
   
    wdir1 = 'results/star/final'

    WDIR = [wdir1]

    nrows,ncols=3,3
    N = nrows*ncols
    fig = py.figure(figsize=(ncols*9,nrows*5))
    axs = {}
    for i in range(N):
        axs[i+1] = py.subplot(nrows,ncols,i+1)

    hand = {}
    j = 0

    Q2 = 100

    flavs = ['u','d','s','ub','db','sb','c','b','g']

    for name in DATA:

        data = DATA[name]['ffhadron']
        load_config(DATA[name]['input'])
        Z=data['Z']

        ZF = data['ZF']

        if name=='chris': color='red'
        if name=='trey':  color='blue'

        for flav in flavs:

            if   flav=='u':  ax = axs[1]
            elif flav=='d':  ax = axs[2]
            elif flav=='s':  ax = axs[3]
            elif flav=='ub': ax = axs[4]
            elif flav=='db': ax = axs[5]
            elif flav=='sb': ax = axs[6]
            elif flav=='c':  ax = axs[7]
            elif flav=='b':  ax = axs[8]
            elif flav=='g':  ax = axs[9]

            ff = ZF[flav]

            mean = np.mean(ff,axis=0)
            std  = np.std (ff,axis=0)

            #--plot each replica
            if mode==0:
                for i in range(len(ff)):
                    hand[name] ,= ax.plot(Z,ff[i],color=color,alpha=0.10)
 
            #--plot average and standard deviation
            if mode==1:

                where = [1 for i in range(len(Z))]

                hand[name]  = ax.fill_between(Z,(mean-std),(mean+std),color=color,alpha=0.7,zorder=5,where=np.array(where))



   
    for i in range(N):
          axs[i+1].set_xlim(0.2,0.9)

          axs[i+1].tick_params(axis='both', which='major', top=True, direction='in',labelsize=30,length=8)
          axs[i+1].tick_params(axis='both', which='minor', top=True, direction='in',labelsize=30,length=4)
          axs[i+1].set_xticks([0.4,0.6,0.8])

          minorLocator = MultipleLocator(0.05)
          axs[i+1].xaxis.set_minor_locator(minorLocator)

    axs [1].tick_params(labelbottom=False)
    axs [2].tick_params(labelbottom=False)
    axs [3].tick_params(labelbottom=False)
    axs [4].tick_params(labelbottom=False)
    axs [5].tick_params(labelbottom=False)
    axs [6].tick_params(labelbottom=False)


    axs [2].tick_params(labelleft=False)
    axs [3].tick_params(labelleft=False)
    axs [5].tick_params(labelleft=False)
    axs [6].tick_params(labelleft=False)
    axs [8].tick_params(labelleft=False)
    axs [9].tick_params(labelleft=False)


    for i in [1,2,3]:
        axs[i].set_ylim(0,0.90)

        axs[i].set_yticks([0.2,0.4,0.6,0.8])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)
   
    for i in [4,5,6]:
        axs[i].set_ylim(0,0.85)

        axs[i].set_yticks([0.2,0.4,0.6,0.8])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)

    for i in [7,8,9]:
        axs[i].set_ylim(0,0.85)

        axs[i].set_yticks([0.2,0.4,0.6,0.8])
        minorLocator = MultipleLocator(0.05)
        axs[i].yaxis.set_minor_locator(minorLocator)

    axs[7].set_xlabel(r'\boldmath$z$',size=40)   
    axs[7].xaxis.set_label_coords(0.95,0.00)
    axs[8].set_xlabel(r'\boldmath$z$',size=40)   
    axs[8].xaxis.set_label_coords(0.95,0.00)
    axs[9].set_xlabel(r'\boldmath$z$',size=40)   
    axs[9].xaxis.set_label_coords(0.95,0.00)

    axs[1].text(0.80,0.83,r'\boldmath$z D^{h^+}_{u}$'       , transform=axs[1].transAxes,size=40)
    axs[2].text(0.80,0.83,r'\boldmath$z D^{h^+}_{d}$'       , transform=axs[2].transAxes,size=40)
    axs[3].text(0.80,0.83,r'\boldmath$z D^{h^+}_{s}$'       , transform=axs[3].transAxes,size=40)
    axs[4].text(0.80,0.83,r'\boldmath$z D^{h^+}_{\bar{u}}$' , transform=axs[4].transAxes,size=40)
    axs[5].text(0.80,0.83,r'\boldmath$z D^{h^+}_{\bar{d}}$' , transform=axs[5].transAxes,size=40)
    axs[6].text(0.80,0.83,r'\boldmath$z D^{h^+}_{\bar{s}}$' , transform=axs[6].transAxes,size=40)
    axs[7].text(0.80,0.83,r'\boldmath$z D^{h^+}_{c}$'       , transform=axs[7].transAxes,size=40)
    axs[8].text(0.80,0.83,r'\boldmath$z D^{h^+}_{b}$'       , transform=axs[8].transAxes,size=40)
    axs[9].text(0.80,0.83,r'\boldmath$z D^{h^+}_{g}$'       , transform=axs[9].transAxes,size=40)

    axs[4].text(0.07,0.85,r'$Q^2 = %d$~'%Q2  + r'\textrm{GeV}'+r'$^2$', transform=axs[4].transAxes,size=25)

    blank ,= axs[1].plot(0,0,color='white',alpha=0.0)

    handles,labels = [],[]
    handles.append(hand['chris'])
    handles.append(hand['trey'])
    labels.append(r'\textrm{\textbf{Chris}}')
    labels.append(r'\textrm{\textbf{Trey}}')
    axs[8].legend(handles,labels,loc='upper left', fontsize = 30, frameon = 0, handletextpad = 0.3, handlelength = 1.0)

    py.tight_layout()
    py.subplots_adjust(hspace=0.02,wspace=0.02,top=0.99,right=0.99,left=0.03)

    filename = '%s/gallery_trey/ffs-hadron-Q2=%3.5f'%(cwd,Q2)
    if mode==1: filename+='-bands'
    filename+='.png'
    for i in range(N):
        axs[i+1] .set_rasterized(True)

    py.savefig(filename)
    print ('Saving figure to %s'%filename)

if __name__=="__main__":

    checkdir('%s/data'%cwd)
    checkdir('%s/gallery_trey'%cwd)

    DATA = {}
    #names = ['chris','trey']
    names = ['trey','chris']
    for name in names: DATA[name] = {}

    trey = np.load('/w/jam-sciwork24/share/replicas/treyand24/treyand-24.npy',allow_pickle=True).item()
    DATA['trey']['input'] = '/w/jam-sciwork24/share/replicas/treyand24/input.py'

    DATA['trey']['params'] = trey['params']
    DATA['trey']['order']  = trey['order']


    chris = np.load('%s/replicas.npy'%wdir,allow_pickle=True).item()
    input_path = '/w/jam-sciwork24/ccocuzza/pol-high-x-Chris/plotting_code/data/chris_input2.py'
    DATA['chris']['input'] = input_path
    print('Loading input file %s'%input_path)

    DATA['chris']['params'] = chris['params']
    DATA['chris']['order']  = chris['order']

    #--load names of all data files
    data_files = os.listdir('%s/data'%cwd)

    pdf_Q2 = 10
    ff_Q2  = 100

    for name in names:
        print('Generating and/or loading results for %s replicas'%name)
        #--set up resman
        load_config(DATA[name]['input'])
        conf['predict'] = True
        conf['bootstrap']=False
        filename = '%s/data/predictions-%s.npy'%(cwd,name)
        if (filename.split('/')[-1] not in data_files) or regen[name]: datasets = True
        else:                                                          datasets = False 
        resman = RESMAN(nworkers=1,datasets=datasets,parallel=False)

        order, params = DATA[name]['order'],DATA[name]['params']

        parman = resman.parman
        parman.order = order

        #--generate predictions
        if (filename.split('/')[-1] not in data_files) or regen[name]: gen_predictions(params,resman,filename)
        DATA[name]['predictions'] = np.load(filename,allow_pickle=True).item()

        #--generate unpolarized PDFs
        filename = '%s/data/pdf-%s-Q2=%3.5f.npy'%(cwd,name,pdf_Q2)
        if (filename.split('/')[-1] not in data_files) or regen[name]: gen_pdf(params,filename,pdf_Q2)
        DATA[name]['pdf'] = np.load(filename,allow_pickle=True).item()
        print('Loaded unpolarized PDF data %s'%filename)

        #--generate fragmentation functions
        hads = ['pion','kaon','hadron']
        for had in hads:
            filename = '%s/data/ff-%s-%s-Q2=%3.5f.npy'%(cwd,had,name,ff_Q2)
            if (filename.split('/')[-1] not in data_files) or regen[name]: gen_ff(had,params,filename,ff_Q2)
            DATA[name]['ff%s'%had] = np.load(filename,allow_pickle=True).item()
            print('Loaded FF %s data %s'%(had,filename))
        

    #--print summary
    print('SUMMARY FOR CHRIS:')
    print_summary(DATA['chris'])
    print() 
    print() 
    print() 
    print('SUMMARY FOR TREY:')
    print_summary(DATA['trey']) 

    #--plot unpolarized PDFs
    plot_pdf     (DATA,Q2=pdf_Q2,mode=0)
    plot_pdf     (DATA,Q2=pdf_Q2,mode=1)

    plot_ffpion  (DATA,Q2=ff_Q2,mode=0)
    plot_ffpion  (DATA,Q2=ff_Q2,mode=1)

    plot_ffkaon  (DATA,Q2=ff_Q2,mode=0)
    plot_ffkaon  (DATA,Q2=ff_Q2,mode=1)

    plot_ffhadron(DATA,Q2=ff_Q2,mode=0)
    plot_ffhadron(DATA,Q2=ff_Q2,mode=1)






