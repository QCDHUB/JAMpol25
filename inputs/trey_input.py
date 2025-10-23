conf={}
conf['bootstrap'] = True
conf['ht']        = False
conf['offpdf']    = False
conf['resample version'] = 0
conf['grid path'] = '/work/JAM/share/fitpack_grids'
conf['ftol']=1e-6
conf['flat par']=False
conf['dglap mode']='truncated'
conf['alphaSmode']='backward'
conf['order'] = 'NLO'
conf['Q20']   = 1.27**2
conf['nuc']     = True
conf['tmc']     = 'AOT'
conf['ht type'] = 'mult'
conf['ww']  = False
conf['dis version'] = 1

conf['datasets']={}

#--lepton-hadron reactions

Q2cut=1.3**2
W2cut=10.0

##--IDIS
conf['datasets']['idis']={}
conf['datasets']['idis']['filters']=[]
conf['datasets']['idis']['filters'].append("Q2>%f"%Q2cut)
conf['datasets']['idis']['filters'].append("W2>%f"%W2cut)
conf['datasets']['idis']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC
conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC
#conf['datasets']['idis']['xlsx'][10003]='idis/expdata/10003.xlsx' # proton   | sigma red     | JLab Hall C (E00-106)
conf['datasets']['idis']['xlsx'][10026]='idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)
conf['datasets']['idis']['xlsx'][10027]='idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)
conf['datasets']['idis']['xlsx'][10028]='idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)
conf['datasets']['idis']['xlsx'][10029]='idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)
conf['datasets']['idis']['xlsx'][10030]='idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-
conf['datasets']['idis']['xlsx'][10031]='idis/expdata/10031.xlsx' # proton   | sigma red     | HERA II CC e+
conf['datasets']['idis']['xlsx'][10032]='idis/expdata/10032.xlsx' # proton   | sigma red     | HERA II CC e-
#conf['datasets']['idis']['xlsx'][10007]='idis/expdata/10007.xlsx' # proton   | sigma red     | HERMES
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC
#conf['datasets']['idis']['xlsx'][10006]='idis/expdata/10006.xlsx' # deuteron | F2            | HERMES
#conf['datasets']['idis']['xlsx'][10002]='idis/expdata/10002.xlsx' # deuteron | F2            | JLab Hall C (E00-106)
#conf['datasets']['idis']['xlsx'][10033]='idis/expdata/10033.xlsx' # n/d      | F2n/F2d       | BONUS
#------------------------------------------------------------------------------------------------------------------
#conf['datasets']['idis']['xlsx'][10041]='idis/expdata/10041.xlsx' # h/d      | F2h/F2d       | JLab Hall C (E03-103)
#conf['datasets']['idis']['xlsx'][10050]='idis/expdata/10050.xlsx' # d/p      | F2d/F2p       | MARATHON
#conf['datasets']['idis']['xlsx'][10051]='idis/expdata/10051.xlsx' # h/t      | F2h/F2t       | MARATHON
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['idis']['norm']={}
#conf['datasets']['idis']['norm'][10002]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10003]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10006]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10007]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10010]={'value':    1.04352e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10011]={'value':    1.04141e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10016]={'value':    9.89544e-01, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10017]={'value':    1.01306e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10020]={'value':    1.02003e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10021]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10033]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10041]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10050]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#conf['datasets']['idis']['norm'][10051]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

conf['datasets']['dy']={}
conf['datasets']['dy']['filters']=[]
conf['datasets']['dy']['filters'].append("Q2>0")
conf['datasets']['dy']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['dy']['xlsx'][10001]='dy/expdata/10001.xlsx'
conf['datasets']['dy']['xlsx'][20001]='dy/expdata/20001.xlsx'
conf['datasets']['dy']['xlsx'][20002]='dy/expdata/20002.xlsx'
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['dy']['norm']={}
conf['datasets']['dy']['norm'][10001]={'value':    1,'fixed':False,'min':   0.5,'max':    1.5}
conf['datasets']['dy']['norm'][20001]={'value':    1,'fixed':False,'min':   0.9,'max':    1.1}
conf['datasets']['dy']['norm'][20002]={'value':    1,'fixed':False,'min':   0.9,'max':    1.1}

##--charge asymmetry
conf['datasets']['wzrv']={}
conf['datasets']['wzrv']['filters']=[]
conf['datasets']['wzrv']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wzrv']['xlsx'][2010]='wzrv/expdata/2010.xlsx'
conf['datasets']['wzrv']['xlsx'][2011]='wzrv/expdata/2011.xlsx'
conf['datasets']['wzrv']['xlsx'][2012]='wzrv/expdata/2012.xlsx'
conf['datasets']['wzrv']['xlsx'][2013]='wzrv/expdata/2013.xlsx'
conf['datasets']['wzrv']['xlsx'][2014]='wzrv/expdata/2014.xlsx'
conf['datasets']['wzrv']['xlsx'][2016]='wzrv/expdata/2016.xlsx'
conf['datasets']['wzrv']['xlsx'][2017]='wzrv/expdata/2017.xlsx'
conf['datasets']['wzrv']['xlsx'][2020]='wzrv/expdata/2020.xlsx'
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wzrv']['norm']={}
#------------------------------------------------------------------------------------------------------------------

##--W asymmetry
conf['datasets']['wasym']={}
conf['datasets']['wasym']['filters']=[]
conf['datasets']['wasym']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wasym']['xlsx'][1000]='wasym/expdata/1000.xlsx'
conf['datasets']['wasym']['xlsx'][1001]='wasym/expdata/1001.xlsx'
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wasym']['norm']={}
#------------------------------------------------------------------------------------------------------------------

##--Z rapidity
conf['datasets']['zrap']={}
conf['datasets']['zrap']['filters']=[]
conf['datasets']['zrap']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['zrap']['xlsx'][1000]='zrap/expdata/1000.xlsx'
conf['datasets']['zrap']['xlsx'][1001]='zrap/expdata/1001.xlsx'
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['zrap']['norm']={}
conf['datasets']['zrap']['norm'][1000]={'value':    1,'fixed':False,'min':   0.5,'max':    1.5}
#------------------------------------------------------------------------------------------------------------------

## jet in unpolarized pp collisions
jet_pt_cut = 8.0 ## GeV
jet_x_cut  = 1 #0.85
conf['datasets']['jet'] = {}
conf['datasets']['jet']['filters'] = []
conf['datasets']['jet']['filters'].append('pT>%f' % jet_pt_cut)
conf['datasets']['jet']['filters'].append('x<%f' % jet_x_cut)
conf['datasets']['jet']['xlsx'] = {}
conf['datasets']['jet']['norm'] = {}
#---------------------------------------------------------------------------------------
conf['datasets']['jet']['xlsx'][10001] = 'jet/expdata/10001.xlsx' ## D0        dataset |
conf['datasets']['jet']['xlsx'][10002] = 'jet/expdata/10002.xlsx' ## CDF       dataset |
conf['datasets']['jet']['xlsx'][10003] = 'jet/expdata/10003.xlsx' ## STAR 2003 dataset |
conf['datasets']['jet']['xlsx'][10004] = 'jet/expdata/10004.xlsx' ## STAR 2004 dataset |
#---------------------------------------------------------------------------------------
conf['datasets']['jet']['norm'][10001] = {'value': 1.0, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['jet']['norm'][10002] = {'value': 1.0, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['jet']['norm'][10003] = {'value': 1.0, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['jet']['norm'][10004] = {'value': 1.0, 'fixed': False, 'min': 0.5, 'max': 1.5}

##--WC
conf['datasets']['wc']={}
conf['datasets']['wc']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wc']['xlsx'][10001]='wc/expdata/10001.xlsx'  # pp | sig | ATLAS7Wm
conf['datasets']['wc']['xlsx'][10002]='wc/expdata/10002.xlsx'  # pp | sig | ATLAS7Wp
conf['datasets']['wc']['xlsx'][20001]='wc/expdata/20001.xlsx'  # pp | Wp+Wm | CMS7WpWm_add
conf['datasets']['wc']['xlsx'][20002]='wc/expdata/20002.xlsx'  # pp | Wp/Wm | CMS7WpWm_ratio
conf['datasets']['wc']['xlsx'][20003]='wc/expdata/20003.xlsx'  # pp | Wp+Wm | CMS13WpWm_add
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wc']['norm']={}
conf['datasets']['wc']['norm'][10001] ={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['wc']['norm'][10002] ={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['wc']['norm'][20003] ={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

#--lepton-lepton reactions

conf['datasets']['sia']={}
conf['datasets']['sia']['filters']=[]
conf['datasets']['sia']['filters'].append('z>0.2 and z<0.9')
conf['datasets']['sia']['xlsx']={}
conf['datasets']['sia']['norm']={}

conf['datasets']['sia']['xlsx'][1001]='sia/expdata/1001.xlsx'  # hadron: pion exp: TASSO
conf['datasets']['sia']['xlsx'][1002]='sia/expdata/1002.xlsx'  # hadron: pion exp: TASSO
conf['datasets']['sia']['xlsx'][1003]='sia/expdata/1003.xlsx'  # hadron: pion exp: TASSO
conf['datasets']['sia']['xlsx'][1004]='sia/expdata/1004.xlsx'  # hadron: pion exp: TASSO
conf['datasets']['sia']['xlsx'][1005]='sia/expdata/1005.xlsx'  # hadron: pion exp: TASSO
conf['datasets']['sia']['xlsx'][1006]='sia/expdata/1006.xlsx'  # hadron: pion exp: TASSO
conf['datasets']['sia']['xlsx'][1007]='sia/expdata/1007.xlsx'  # hadron: pion exp: TPC
conf['datasets']['sia']['xlsx'][1008]='sia/expdata/1008.xlsx'  # hadron: pion exp: TPC
conf['datasets']['sia']['xlsx'][1012]='sia/expdata/1012.xlsx'  # hadron: pion exp: HRS
conf['datasets']['sia']['xlsx'][1013]='sia/expdata/1013.xlsx'  # hadron: pion exp: TOPAZ
conf['datasets']['sia']['xlsx'][1014]='sia/expdata/1014.xlsx'  # hadron: pion exp: SLD
conf['datasets']['sia']['xlsx'][1018]='sia/expdata/1018.xlsx'  # hadron: pion exp: ALEPH
conf['datasets']['sia']['xlsx'][1019]='sia/expdata/1019.xlsx'  # hadron: pion exp: OPAL
conf['datasets']['sia']['xlsx'][1025]='sia/expdata/1025.xlsx'  # hadron: pion exp: DELPHI
conf['datasets']['sia']['xlsx'][1028]='sia/expdata/1028.xlsx'  # hadron: pion exp: BABAR
conf['datasets']['sia']['xlsx'][1029]='sia/expdata/1029.xlsx'  # hadron: pion exp: BELL
conf['datasets']['sia']['xlsx'][1030]='sia/expdata/1030.xlsx'  # hadron: pion exp: ARGUS
conf['datasets']['sia']['norm'][1001]={'value':    1.10478e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1002]={'value':    9.82581e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1003]={'value':    1.03054e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1005]={'value':    1.03419e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1006]={'value':    9.79162e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1014]={'value':    9.97770e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1018]={'value':    1.02378e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1028]={'value':    9.76001e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1029]={'value':    8.68358e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1030]={'value':    1.01509e+00,'fixed':False,'min':0.5,'max':1.5}

##--SIA pion HQ
conf['datasets']['sia']['xlsx'][1010]='sia/expdata/1010.xlsx'  # hadron: pion exp: TPC(c)
conf['datasets']['sia']['xlsx'][1011]='sia/expdata/1011.xlsx'  # hadron: pion exp: TPC(b)
conf['datasets']['sia']['xlsx'][1016]='sia/expdata/1016.xlsx'  # hadron: pion exp: SLD(c)
conf['datasets']['sia']['xlsx'][1017]='sia/expdata/1017.xlsx'  # hadron: pion exp: SLD(b)
conf['datasets']['sia']['xlsx'][1023]='sia/expdata/1023.xlsx'  # hadron: pion exp: OPAL(c)
conf['datasets']['sia']['xlsx'][1024]='sia/expdata/1024.xlsx'  # hadron: pion exp: OPAL(b)
conf['datasets']['sia']['xlsx'][1027]='sia/expdata/1027.xlsx'  # hadron: pion exp: DELPHI(b)
conf['datasets']['sia']['norm'][1016]={'value':    1.18920e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1017]={'value':    1.00345e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1023]={'value':    1.33434e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1024]={'value':    1.19151e+00,'fixed':False,'min':0.5,'max':1.5}

##--SIA kaon
conf['datasets']['sia']['xlsx'][2030]='sia/expdata/2030.xlsx'  # hadron: kaon exp: ARGUS
conf['datasets']['sia']['xlsx'][2028]='sia/expdata/2028.xlsx'  # hadron: kaon exp: BABAR
conf['datasets']['sia']['xlsx'][2029]='sia/expdata/2029.xlsx'  # hadron: kaon exp: BELL
conf['datasets']['sia']['xlsx'][2001]='sia/expdata/2001.xlsx'  # hadron: kaon exp: TASSO
conf['datasets']['sia']['xlsx'][2002]='sia/expdata/2002.xlsx'  # hadron: kaon exp: TASSO
conf['datasets']['sia']['xlsx'][2003]='sia/expdata/2003.xlsx'  # hadron: kaon exp: TASSO
conf['datasets']['sia']['xlsx'][2004]='sia/expdata/2004.xlsx'  # hadron: kaon exp: TASSO
conf['datasets']['sia']['xlsx'][2005]='sia/expdata/2005.xlsx'  # hadron: kaon exp: TASSO
conf['datasets']['sia']['xlsx'][2006]='sia/expdata/2006.xlsx'  # hadron: kaon exp: TASSO
conf['datasets']['sia']['xlsx'][2007]='sia/expdata/2007.xlsx'  # hadron: kaon exp: TPC
conf['datasets']['sia']['xlsx'][2008]='sia/expdata/2008.xlsx'  # hadron: kaon exp: TPC
conf['datasets']['sia']['xlsx'][2012]='sia/expdata/2012.xlsx'  # hadron: kaon exp: HRS
conf['datasets']['sia']['xlsx'][2013]='sia/expdata/2013.xlsx'  # hadron: kaon exp: TOPAZ
conf['datasets']['sia']['xlsx'][2014]='sia/expdata/2014.xlsx'  # hadron: kaon exp: SLD
conf['datasets']['sia']['xlsx'][2018]='sia/expdata/2018.xlsx'  # hadron: kaon exp: ALEPH
conf['datasets']['sia']['xlsx'][2019]='sia/expdata/2019.xlsx'  # hadron: kaon exp: OPAL
conf['datasets']['sia']['xlsx'][2025]='sia/expdata/2025.xlsx'  # hadron: kaon exp: DELPHI
conf['datasets']['sia']['xlsx'][2031]='sia/expdata/2031.xlsx'  # hadron: kaon exp: DELPHI
conf['datasets']['sia']['norm'][2030]={'value':    1.00482e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2028]={'value':    9.97435e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2029]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2002]={'value':    9.83394e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2003]={'value':    9.94421e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2005]={'value':    9.92876e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2014]={'value':    9.12186e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2018]={'value':    8.88453e-01,'fixed':False,'min':0.5,'max':1.5}

##--SIA kaon HQ
conf['datasets']['sia']['xlsx'][2016]='sia/expdata/2016.xlsx'  # hadron: kaon exp: SLD(c)
conf['datasets']['sia']['xlsx'][2017]='sia/expdata/2017.xlsx'  # hadron: kaon exp: SLD(b)
conf['datasets']['sia']['xlsx'][2023]='sia/expdata/2023.xlsx'  # hadron: kaon exp: OPAL(c)
conf['datasets']['sia']['xlsx'][2024]='sia/expdata/2024.xlsx'  # hadron: kaon exp: OPAL(b)
conf['datasets']['sia']['xlsx'][2027]='sia/expdata/2027.xlsx'  # hadron: kaon exp: DELPHI(b)
conf['datasets']['sia']['norm'][2016]={'value':    1.03996e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2017]={'value':    1.00015e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2023]={'value':    1.35922e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2024]={'value':    1.38163e+00,'fixed':False,'min':0.5,'max':1.5}

##--SIA hadrons
conf['datasets']['sia']['xlsx'][4000]='sia/expdata/4000.xlsx'  # hadron: hadrons exp: ALEPH
conf['datasets']['sia']['xlsx'][4001]='sia/expdata/4001.xlsx'  # hadron: hadrons exp: DELPHI
conf['datasets']['sia']['xlsx'][4013]='sia/expdata/4013.xlsx'  # hadron: hadrons exp: DELPHI(b)
conf['datasets']['sia']['xlsx'][4002]='sia/expdata/4002.xlsx'  # hadron: hadrons exp: SLD
conf['datasets']['sia']['xlsx'][4014]='sia/expdata/4014.xlsx'  # hadron: hadrons exp: SLD(c)
conf['datasets']['sia']['xlsx'][4015]='sia/expdata/4015.xlsx'  # hadron: hadrons exp: SLD(b)
#conf['datasets']['sia']['xlsx'][4003]='sia/expdata/4003.xlsx'  # hadron: hadrons exp: TASSO
#conf['datasets']['sia']['xlsx'][4008]='sia/expdata/4008.xlsx'  # hadron: hadrons exp: TASSO
#conf['datasets']['sia']['xlsx'][4009]='sia/expdata/4009.xlsx'  # hadron: hadrons exp: TASSO
#conf['datasets']['sia']['xlsx'][4010]='sia/expdata/4010.xlsx'  # hadron: hadrons exp: TASSO
conf['datasets']['sia']['xlsx'][4011]='sia/expdata/4011.xlsx'  # hadron: hadrons exp: TASSO
conf['datasets']['sia']['xlsx'][4012]='sia/expdata/4012.xlsx'  # hadron: hadrons exp: TASSO
conf['datasets']['sia']['xlsx'][4004]='sia/expdata/4004.xlsx'  # hadron: hadrons exp: TPC
conf['datasets']['sia']['xlsx'][4005]='sia/expdata/4005.xlsx'  # hadron: hadrons exp: OPAL(b)
conf['datasets']['sia']['xlsx'][4006]='sia/expdata/4006.xlsx'  # hadron: hadrons exp: OPAL(c)
conf['datasets']['sia']['xlsx'][4007]='sia/expdata/4007.xlsx'  # hadron: hadrons exp: OPAL
conf['datasets']['sia']['norm'][4000]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4002]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sia']['norm'][4003]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sia']['norm'][4008]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sia']['norm'][4009]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sia']['norm'][4010]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4011]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4012]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4014]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4015]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}

##--SIDIS
W2cut_SIDIS = 20.0
conf['datasets']['sidis']={}
conf['datasets']['sidis']['filters']=[]
conf['datasets']['sidis']['filters'].append("Q2>%f"%Q2cut)
conf['datasets']['sidis']['filters'].append("W2>%f"%W2cut)
conf['datasets']['sidis']['filters'].append("W2SIDIS>%f"%W2cut_SIDIS)
conf['datasets']['sidis']['filters'].append('Z>0.2 and Z<0.8')
conf['datasets']['sidis']['xlsx']={}
conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx' # deuteron , mult , pi+ , COMPASS
conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx' # deuteron , mult , pi- , COMPASS
conf['datasets']['sidis']['xlsx'][2005]='sidis/expdata/2005.xlsx' # deuteron , mult , K+  , COMPASS
conf['datasets']['sidis']['xlsx'][2006]='sidis/expdata/2006.xlsx' # deuteron , mult , K-  , COMPASS
conf['datasets']['sidis']['xlsx'][3000]='sidis/expdata/3000.xlsx' # deuteron , mult , h+  , COMPASS
conf['datasets']['sidis']['xlsx'][3001]='sidis/expdata/3001.xlsx' # deuteron , mult , h-  , COMPASS
conf['datasets']['sidis']['norm']={}
#conf['datasets']['sidis']['norm'][1005]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sidis']['norm'][1006]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sidis']['norm'][2005]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sidis']['norm'][2006]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sidis']['norm'][3000]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}
#conf['datasets']['sidis']['norm'][3001]={'value':    1.0000000,'fixed':False,'min':0.5,'max':1.5}

#--parameters
conf['params']={}

#--pdf parameters
conf['params']['pdf']={}

conf['params']['pdf']['g1 N']    ={'value':    3.87592e-01   , 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['g1 a']    ={'value':   -6.23068169e-01, 'min':  -1.9, 'max':    10, 'fixed': False}
conf['params']['pdf']['g1 b']    ={'value':    9.25741583e+00, 'min':     0, 'max':    20, 'fixed': False}
conf['params']['pdf']['g1 c']    ={'value':    0.00000000e+00, 'min':  -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['g1 d']    ={'value':    0.00000000e+00, 'min':  -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['uv1 N']   ={'value':    3.47549e-01   , 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['uv1 a']   ={'value':   -1.21835956e-01, 'min':  -0.9, 'max':    10, 'fixed': False}
conf['params']['pdf']['uv1 b']   ={'value':    3.20766744e+00, 'min':     0, 'max':    20, 'fixed': False}
conf['params']['pdf']['uv1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['uv1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['dv1 N']   ={'value':    1.52089e-01   , 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['dv1 a']   ={'value':   -2.39874967e-01, 'min':  -0.9, 'max':    10, 'fixed': False}
conf['params']['pdf']['dv1 b']   ={'value':    3.83902620e+00, 'min':     0, 'max':    20, 'fixed': False}
conf['params']['pdf']['dv1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['dv1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

#conf['params']['pdf']['mix N']    ={'value':    0.00000000e+00   , 'min':  -0.1, 'max': 5.0, 'fixed': False}
#conf['params']['pdf']['mix a']    ={'value':    2.00000000e+00   , 'min':   1.0, 'max':10.0, 'fixed': False}

conf['params']['pdf']['db1 N']   ={'value':    3.67609928e-02, 'min':     0, 'max':      1, 'fixed': False}
conf['params']['pdf']['db1 a']   ={'value':   -8.41360631e-01, 'min':     -1, 'max':     10, 'fixed': False}
conf['params']['pdf']['db1 b']   ={'value':    5.31285539e+00, 'min':      0, 'max':     20, 'fixed': False}
conf['params']['pdf']['db1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['db1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['ub1 N']   ={'value':    1.95464789e-02, 'min':     0, 'max':      1, 'fixed': False}
conf['params']['pdf']['ub1 a']   ={'value':   -9.93659187e-01, 'min':     -1, 'max':     10, 'fixed': False}
conf['params']['pdf']['ub1 b']   ={'value':    8.38905814e+00, 'min':      0, 'max':     20, 'fixed': False}
conf['params']['pdf']['ub1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['ub1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['s1 N']    ={'value':       0.00000e+00, 'min':     0, 'max':     1, 'fixed': True }
conf['params']['pdf']['s1 a']    ={'value':    1.34706224e-01, 'min':  -0.9, 'max':    10, 'fixed': False}
conf['params']['pdf']['s1 b']    ={'value':    6.00759596e+00, 'min':     0, 'max':    20, 'fixed': False}
#conf['params']['pdf']['s1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
#conf['params']['pdf']['s1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['sb1 N']   ={'value':    7.46109845e-07, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['sb1 a']   ={'value':    3.83495317e-01, 'min':  -0.9, 'max':    10, 'fixed': False}
conf['params']['pdf']['sb1 b']   ={'value':    4.61209808e+00, 'min':     0, 'max':    20, 'fixed': False}
#conf['params']['pdf']['sb1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
#conf['params']['pdf']['sb1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['sea1 N']   ={'value':    5.71081196e-03, 'min':    0,  'max':      1, 'fixed': False}
conf['params']['pdf']['sea1 a']   ={'value':   -1.36329697e+00, 'min':  -1.9,  'max':     10, 'fixed': False}
conf['params']['pdf']['sea1 b']   ={'value':    4.74721050e+00, 'min':     0,  'max':     20, 'fixed': False}
conf['params']['pdf']['sea1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['sea1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['sea2 N']   ={'value':    5.71081196e-03, 'min':    0,  'max':      1, 'fixed': False}
conf['params']['pdf']['sea2 a']   ={'value':   -1.36329697e+00, 'min':  -1.9,  'max':     10, 'fixed': False}
conf['params']['pdf']['sea2 b']   ={'value':    4.74721050e+00, 'min':     0,  'max':     20, 'fixed': False}
#conf['params']['pdf']['sea2 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
#conf['params']['pdf']['sea2 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

#--kaon fragmentation
conf['params']['ffpion']={}
conf['ffpion parametrization']=0

conf['params']['ffpion']['g1 N']  ={'value':   2.33320e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['g1 a']  ={'value':   1.48737e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['g1 b']  ={'value':   9.62755e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['u1 N']  ={'value':   4.03672e-02, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['u1 a']  ={'value':   1.26356e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['u1 b']  ={'value':   1.62596e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['d1 N']  ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['d1 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['d1 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['s1 N']  ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['s1 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['s1 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['c1 N']  ={'value':   7.76923e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['c1 a']  ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['c1 b']  ={'value':   2.50452e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['b1 N']  ={'value':   5.66971e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['b1 a']  ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['b1 b']  ={'value':   3.41727e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['ub1 N'] ={'value':   4.03672e-02, 'min':   0  , 'max':     1,'fixed':'d1 N'}
conf['params']['ffpion']['ub1 a'] ={'value':   2.26356e+00, 'min': -1.9 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffpion']['ub1 b'] ={'value':   1.62596e+00, 'min':  0   , 'max':    20,'fixed':'d1 b'}

conf['params']['ffpion']['db1 N'] ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':'u1 N'}
conf['params']['ffpion']['db1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffpion']['db1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'u1 b'}

conf['params']['ffpion']['sb1 N'] ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':'s1 N'}
conf['params']['ffpion']['sb1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'s1 a'}
conf['params']['ffpion']['sb1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'s1 b'}

conf['params']['ffpion']['cb1 N'] ={'value':   7.76923e-01, 'min':   0  , 'max':     1,'fixed':'c1 N'}
conf['params']['ffpion']['cb1 a'] ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':'c1 a'}
conf['params']['ffpion']['cb1 b'] ={'value':   2.50452e+00, 'min':  0   , 'max':    20,'fixed':'c1 b'}

conf['params']['ffpion']['bb1 N'] ={'value':   5.66971e-01, 'min':   0  , 'max':     1,'fixed':'b1 N'}
conf['params']['ffpion']['bb1 a'] ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':'b1 a'}
conf['params']['ffpion']['bb1 b'] ={'value':   3.41727e+00, 'min':  0   , 'max':    20,'fixed':'b1 b'}

conf['params']['ffpion']['u2 N']  ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['u2 a']  ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['u2 b']  ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['d2 N']  ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffpion']['d2 a']  ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffpion']['d2 b']  ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffpion']['ub2 N'] ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':'d2 N'}
conf['params']['ffpion']['ub2 a'] ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffpion']['ub2 b'] ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':'d2 b'}

conf['params']['ffpion']['db2 N'] ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':'u2 N'}
conf['params']['ffpion']['db2 a'] ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':'u2 a'}
conf['params']['ffpion']['db2 b'] ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':'u2 b'}

#--kaon fragmentation
conf['params']['ffkaon']={}
conf['ffkaon parametrization']=0

conf['params']['ffkaon']['g1 N']  ={'value':   2.33320e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['g1 a']  ={'value':   1.48737e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['g1 b']  ={'value':   9.62755e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['u1 N']  ={'value':   4.03672e-02, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['u1 a']  ={'value':   1.26356e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['u1 b']  ={'value':   1.62596e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['d1 N']  ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['d1 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['d1 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['s1 N']  ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':'d1 N'}
conf['params']['ffkaon']['s1 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffkaon']['s1 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    20,'fixed':'d1 b'}

conf['params']['ffkaon']['c1 N']  ={'value':   7.76923e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['c1 a']  ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['c1 b']  ={'value':   2.50452e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['b1 N']  ={'value':   5.66971e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['b1 a']  ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['b1 b']  ={'value':   3.41727e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['ub1 N'] ={'value':   4.03672e-02, 'min':   0  , 'max':     1,'fixed':'d1 N'}
conf['params']['ffkaon']['ub1 a'] ={'value':   2.26356e+00, 'min': -1.9 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffkaon']['ub1 b'] ={'value':   1.62596e+00, 'min':  0   , 'max':    20,'fixed':'d1 b'}

conf['params']['ffkaon']['db1 N'] ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':'d1 N'}
conf['params']['ffkaon']['db1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffkaon']['db1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'d1 b'}

conf['params']['ffkaon']['sb1 N'] ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['sb1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['sb1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['cb1 N'] ={'value':   7.76923e-01, 'min':   0  , 'max':     1,'fixed':'c1 N'}
conf['params']['ffkaon']['cb1 a'] ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':'c1 a'}
conf['params']['ffkaon']['cb1 b'] ={'value':   2.50452e+00, 'min':  0   , 'max':    20,'fixed':'c1 b'}

conf['params']['ffkaon']['bb1 N'] ={'value':   5.66971e-01, 'min':   0  , 'max':     1,'fixed':'b1 N'}
conf['params']['ffkaon']['bb1 a'] ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':'b1 a'}
conf['params']['ffkaon']['bb1 b'] ={'value':   3.41727e+00, 'min':  0   , 'max':    20,'fixed':'b1 b'}

conf['params']['ffkaon']['u2 N']  ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['u2 a']  ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['u2 b']  ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['d2 N']  ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':False}
conf['params']['ffkaon']['d2 a']  ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['d2 b']  ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffkaon']['s2 N']  ={'value':   2.51671e-01, 'min':   0  , 'max':     1,'fixed':'d2 N'}
conf['params']['ffkaon']['s2 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffkaon']['s2 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    20,'fixed':'d2 b'}

conf['params']['ffkaon']['ub2 N'] ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':'d2 N'}
conf['params']['ffkaon']['ub2 a'] ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffkaon']['ub2 b'] ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':'d2 b'}

conf['params']['ffkaon']['db2 N'] ={'value':   0.00000e+00, 'min':   0  , 'max':     1,'fixed':'d2 N'}
conf['params']['ffkaon']['db2 a'] ={'value':   1.76877e-01, 'min': -1.9 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffkaon']['db2 b'] ={'value':   4.81521e+00, 'min':  0   , 'max':    20,'fixed':'d2 b'}

#--unidentified hadrons
conf['params']['ffhadron']={}
conf['ffhadron parametrization']='sum2'

conf['params']['ffhadron']['g1 N']  ={'value':   2.33320e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['g1 a']  ={'value':   1.48737e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['g1 b']  ={'value':   9.62755e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffhadron']['u1 N']  ={'value':   4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['u1 a']  ={'value':   1.26356e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['u1 b']  ={'value':   1.62596e+00, 'min':  0   , 'max':    20,'fixed':False}
conf['params']['ffhadron']['u1 c']  ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':False}
conf['params']['ffhadron']['u1 d']  ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':False}

conf['params']['ffhadron']['d1 N']  ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'u1 N'}
conf['params']['ffhadron']['d1 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffhadron']['d1 b']  ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'u1 b'}
conf['params']['ffhadron']['d1 c']  ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 c'}
conf['params']['ffhadron']['d1 d']  ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 d'}

conf['params']['ffhadron']['s1 N']  ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'u1 N'}
conf['params']['ffhadron']['s1 a']  ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffhadron']['s1 b']  ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'u1 b'}
conf['params']['ffhadron']['s1 c']  ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 c'}
conf['params']['ffhadron']['s1 d']  ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 d'}

conf['params']['ffhadron']['c1 N']  ={'value':   7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['c1 a']  ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['c1 b']  ={'value':   2.50452e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffhadron']['b1 N']  ={'value':   5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['b1 a']  ={'value':  -1.80000e+00, 'min': -1.9 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['b1 b']  ={'value':   3.41727e+00, 'min':  0   , 'max':    20,'fixed':False}

conf['params']['ffhadron']['ub1 N'] ={'value':   4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['ub1 a'] ={'value':   2.26356e+00, 'min': -1.9 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffhadron']['ub1 b'] ={'value':   1.62596e+00, 'min':  0   , 'max':    20,'fixed':False}
conf['params']['ffhadron']['ub1 c'] ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 c'}
conf['params']['ffhadron']['ub1 d'] ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 d'}

conf['params']['ffhadron']['db1 N'] ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'ub1 N'}
conf['params']['ffhadron']['db1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffhadron']['db1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'ub1 b'}
conf['params']['ffhadron']['db1 c'] ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 c'}
conf['params']['ffhadron']['db1 d'] ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 d'}

conf['params']['ffhadron']['sb1 N'] ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'ub1 N'}
conf['params']['ffhadron']['sb1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffhadron']['sb1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'ub1 b'}
conf['params']['ffhadron']['sb1 c'] ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 c'}
conf['params']['ffhadron']['sb1 d'] ={'value':   0.00000e+00, 'min': -100 , 'max':   100,'fixed':'u1 d'}

conf['params']['ffhadron']['cb1 N'] ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'c1 N'}
conf['params']['ffhadron']['cb1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'c1 a'}
conf['params']['ffhadron']['cb1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'c1 b'}

conf['params']['ffhadron']['bb1 N'] ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'b1 N'}
conf['params']['ffhadron']['bb1 a'] ={'value':  -1.43444e+00, 'min': -1.9 , 'max':    10,'fixed':'b1 a'}
conf['params']['ffhadron']['bb1 b'] ={'value':   1.65143e+00, 'min':  0   , 'max':    20,'fixed':'b1 b'}
