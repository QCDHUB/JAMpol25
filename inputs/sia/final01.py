import os
conf={}

conf['enforce limits'] = False

conf['beta version'] = 1

conf['sia version'] = 1
data_version = 1

#--fitting setups
conf['bootstrap']=True
conf['flat par'] =True

conf['grid path'] = '/w/jam-sciwork24/share/fitpack_grids'
conf['dis_version'] = 1

#--setups for DGLAP
conf['dglap mode']='truncated'
conf['alphaSmode']='backward'
conf['order'] = 'NLO'
conf['Q20']   = 1.27**2

#--nuclear smearing
conf['nuc']   = True

#--datasets

conf['datasets']={}

#--lepton-hadron reactions

Q2cut=1.3**2
W2cut=3.5

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
conf['datasets']['idis']['xlsx'][10003]='idis/expdata/10003.xlsx' # proton   | sigma red     | JLab Hall C (E00-106)
conf['datasets']['idis']['xlsx'][10026]='idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)
conf['datasets']['idis']['xlsx'][10027]='idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)
conf['datasets']['idis']['xlsx'][10028]='idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)
conf['datasets']['idis']['xlsx'][10029]='idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)
conf['datasets']['idis']['xlsx'][10030]='idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-
conf['datasets']['idis']['xlsx'][10031]='idis/expdata/10031.xlsx' # proton   | sigma red     | HERA II CC e+
conf['datasets']['idis']['xlsx'][10032]='idis/expdata/10032.xlsx' # proton   | sigma red     | HERA II CC e-
conf['datasets']['idis']['xlsx'][10007]='idis/expdata/10007.xlsx' # proton   | sigma red     | HERMES
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC
conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS
conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC
conf['datasets']['idis']['xlsx'][10006]='idis/expdata/10006.xlsx' # deuteron | F2            | HERMES
conf['datasets']['idis']['xlsx'][10002]='idis/expdata/10002.xlsx' # deuteron | F2            | JLab Hall C (E00-106)
conf['datasets']['idis']['xlsx'][10033]='idis/expdata/10033.xlsx' # n/d      | F2n/F2d       | BONUS
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['idis']['xlsx'][10041]='idis/expdata/10041.xlsx' # h/d      | F2h/F2d       | JLab Hall C (E03-103)
conf['datasets']['idis']['xlsx'][10050]='idis/expdata/10050.xlsx' # d/p      | F2d/F2p       | MARATHON
conf['datasets']['idis']['xlsx'][10051]='idis/expdata/10051.xlsx' # h/t      | F2h/F2t       | MARATHON
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['idis']['norm']={}
conf['datasets']['idis']['norm'][10002]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10003]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10010]={'value':    1.04352e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10011]={'value':    1.04141e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10016]={'value':    9.89544e-01, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10017]={'value':    1.01306e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10020]={'value':    1.02003e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10021]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10033]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10041]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10050]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['idis']['norm'][10051]={'value':    1.00000e+00, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

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

##--PIDIS
pQ2cut=1.3**2
pW2cut=10.0

conf['datasets']['pidis']={}
conf['datasets']['pidis']['filters']=[]
conf['datasets']['pidis']['filters'].append("Q2>%f"%pQ2cut) 
conf['datasets']['pidis']['filters'].append("W2>%f"%pW2cut) 
conf['datasets']['pidis']['xlsx']={}
conf['datasets']['pidis']['norm']={}
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10002]='pidis/expdata/10002.xlsx' # 10002 | proton   | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10003]='pidis/expdata/10003.xlsx' # 10003 | proton   | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10004]='pidis/expdata/10004.xlsx' # 10004 | proton   | A1   | EMC             |          |
conf['datasets']['pidis']['xlsx'][10005]='pidis/expdata/10005.xlsx' # 10005 | neutron  | A1   | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10007]='pidis/expdata/10007.xlsx' # 10007 | proton   | Apa  | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10008]='pidis/expdata/10008.xlsx' # 10008 | proton   | A2   | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10017]='pidis/expdata/10017.xlsx' # 10017 | proton   | Apa  | JLabHB(EG1DVCS) |          |
conf['datasets']['pidis']['xlsx'][10022]='pidis/expdata/10022.xlsx' # 10022 | proton   | Apa  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10023]='pidis/expdata/10023.xlsx' # 10023 | proton   | Ape  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10028]='pidis/expdata/10028.xlsx' # 10028 | proton   | Ape  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10029]='pidis/expdata/10029.xlsx' # 10029 | proton   | Apa  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10031]='pidis/expdata/10031.xlsx' # 10031 | proton   | Atpe | SLAC(E155x)     |          |
conf['datasets']['pidis']['xlsx'][10032]='pidis/expdata/10032.xlsx' # 10032 | proton   | Apa  | SLACE80E130     |          |
conf['datasets']['pidis']['xlsx'][10035]='pidis/expdata/10035.xlsx' # 10035 | proton   | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10036]='pidis/expdata/10036.xlsx' # 10036 | proton   | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10043]='pidis/expdata/10043.xlsx' # 10043 | proton   | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['datasets']['pidis']['xlsx'][10044]='pidis/expdata/10044.xlsx' # 10044 | proton   | Apa  | JLabHB(EG1b)    | E =5 GeV |
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10001]='pidis/expdata/10001.xlsx' # 10001 | deuteron | A1   | COMPASS         |          |
conf['datasets']['pidis']['xlsx'][10006]='pidis/expdata/10006.xlsx' # 10006 | deuteron | Apa  | HERMES          |          |
conf['datasets']['pidis']['xlsx'][10016]='pidis/expdata/10016.xlsx' # 10016 | deuteron | Apa  | JLabHB(EG1DVCS) |          |
conf['datasets']['pidis']['xlsx'][10020]='pidis/expdata/10020.xlsx' # 10020 | deuteron | Ape  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10021]='pidis/expdata/10021.xlsx' # 10021 | deuteron | Apa  | SLAC(E143)      |          |
conf['datasets']['pidis']['xlsx'][10026]='pidis/expdata/10026.xlsx' # 10026 | deuteron | Ape  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10027]='pidis/expdata/10027.xlsx' # 10027 | deuteron | Apa  | SLAC(E155)      |          |
conf['datasets']['pidis']['xlsx'][10030]='pidis/expdata/10030.xlsx' # 10030 | deuteron | Atpe | SLAC(E155x)     |          |
conf['datasets']['pidis']['xlsx'][10033]='pidis/expdata/10033.xlsx' # 10033 | deuteron | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10034]='pidis/expdata/10034.xlsx' # 10034 | deuteron | A1   | SMC             |          |
conf['datasets']['pidis']['xlsx'][10039]='pidis/expdata/10039.xlsx' # 10039 | deuteron | Apa  | JLabHB(EG1b)    | E =4 GeV |
conf['datasets']['pidis']['xlsx'][10040]='pidis/expdata/10040.xlsx' # 10040 | deuteron | Apa  | JLabHB(EG1b)    | E =5 GeV |
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][10009]='pidis/expdata/10009.xlsx' # 10009 | helium   | Apa  | JLabHA(E01-012) |          |
conf['datasets']['pidis']['xlsx'][10010]='pidis/expdata/10010.xlsx' # 10010 | helium   | Apa  | JLabHA(E06-014) |          |
conf['datasets']['pidis']['xlsx'][10011]='pidis/expdata/10011.xlsx' # 10011 | helium   | Ape  | JLabHA(E06-014) |          |
conf['datasets']['pidis']['xlsx'][10012]='pidis/expdata/10012.xlsx' # 10012 | helium   | Apa  | JLabHA(E97-103) | wrong?   |
conf['datasets']['pidis']['xlsx'][10013]='pidis/expdata/10013.xlsx' # 10013 | helium   | Ape  | JLabHA(E97-103) | wrong?   |
conf['datasets']['pidis']['xlsx'][10014]='pidis/expdata/10014.xlsx' # 10014 | helium   | Apa  | JLabHA(E99-117) |          |
conf['datasets']['pidis']['xlsx'][10015]='pidis/expdata/10015.xlsx' # 10015 | helium   | Ape  | JLabHA(E99-117) |          |
conf['datasets']['pidis']['xlsx'][10018]='pidis/expdata/10018.xlsx' # 10018 | helium   | A1   | SLAC(E142)      |          |
conf['datasets']['pidis']['xlsx'][10019]='pidis/expdata/10019.xlsx' # 10019 | helium   | A2   | SLAC(E142)      |          |
conf['datasets']['pidis']['xlsx'][10024]='pidis/expdata/10024.xlsx' # 10024 | helium   | Ape  | SLAC(E154)      |          |
conf['datasets']['pidis']['xlsx'][10025]='pidis/expdata/10025.xlsx' # 10025 | helium   | Apa  | SLAC(E154)      |          |
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['xlsx'][20001]='pidis/expdata/20001.xlsx' #  | SANE    |          |
conf['datasets']['pidis']['xlsx'][20002]='pidis/expdata/20002.xlsx' #  | SANE    |          |
conf['datasets']['pidis']['xlsx'][20003]='pidis/expdata/20003.xlsx' #  | SANE    |          |
conf['datasets']['pidis']['xlsx'][20004]='pidis/expdata/20004.xlsx' #  | SANE    |          |
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['norm'][10002]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10003]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10004]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10005]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10007]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10022]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10023]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10029]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10031]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10032]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['norm'][10006]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10020]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10021]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10001]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10027]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10030]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['norm'][10018]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10025]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['norm'][10008]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10010]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10011]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10014]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10015]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10019]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10024]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
#--------------------------------------------------------------------------------------------------------------------------
conf['datasets']['pidis']['norm'][10016]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10017]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10040]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10044]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10039]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['pidis']['norm'][10043]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}

#--sum rules
conf['datasets']['SU23']={}
conf['datasets']['SU23']['filters']=[]
conf['datasets']['SU23']['xlsx']={}
conf['datasets']['SU23']['norm']={}
#---------------------------------------------------------------------------------------------------------------------------
conf['datasets']['SU23']['xlsx'][10001]='SU23/expdata/10001.xlsx' #  SU2   | nominal 
conf['datasets']['SU23']['xlsx'][10002]='SU23/expdata/10002.xlsx' #  SU3   | nominal 
conf['datasets']['SU23']['xlsx'][20001]='SU23/expdata/20001.xlsx' #  SU2   | JAM17
conf['datasets']['SU23']['xlsx'][20002]='SU23/expdata/20002.xlsx' #  SU3   | JAM17



#--hadron-hadron reactions

##--DY 
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
#------------------------------------------------------------------------------------------------------------------

##--charge asymmetry 
conf['datasets']['wzrv']={}
conf['datasets']['wzrv']['filters']=[]
conf['datasets']['wzrv']['xlsx']={}
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['wzrv']['xlsx'][1000]='wzrv/expdata/1000.xlsx'
conf['datasets']['wzrv']['xlsx'][1020]='wzrv/expdata/1020.xlsx'
conf['datasets']['wzrv']['xlsx'][1021]='wzrv/expdata/1021.xlsx'
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
conf['datasets']['wzrv']['norm'][1000]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['wzrv']['norm'][1020]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
conf['datasets']['wzrv']['norm'][1021]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
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

#--jets
conf['datasets']['jet'] = {}
conf['datasets']['jet']['filters'] = []
conf['datasets']['jet']['filters'].append("pT>8")
conf['datasets']['jet']['xlsx'] = {}
conf['datasets']['jet']['xlsx'][10001] = 'jet/expdata/10001.xlsx' ## D0 dataset
conf['datasets']['jet']['xlsx'][10002] = 'jet/expdata/10002.xlsx' ## CDF dataset
conf['datasets']['jet']['xlsx'][10003] = 'jet/expdata/10003.xlsx' ## STAR MB dataset
conf['datasets']['jet']['xlsx'][10004] = 'jet/expdata/10004.xlsx' ## STAR HT dataset
conf['datasets']['jet']['norm'] = {}
conf['datasets']['jet']['norm'][10001]={'value':    1,'fixed':False,'min':   0.0,'max':    2.0}
conf['datasets']['jet']['norm'][10002]={'value':    1,'fixed':False,'min':   0.0,'max':    2.0}
conf['datasets']['jet']['norm'][10003]={'value':    1,'fixed':False,'min':   0.0,'max':    2.0}
conf['datasets']['jet']['norm'][10004]={'value':    1,'fixed':False,'min':   0.0,'max':    2.0}
conf['channels'] = {'q_qp': 0, 'q_qbp': 1, 'q_q': 2, 'q_qb': 3, 'q_g': 4, 'g_g': 5}
conf['qr'] = {}
conf['qr']['parameter'] = {'epsilon': 1e-11, 'block_size': 10, 'power_scheme': 1}
conf['jet_qr_fit'] = {'method': 'fixed', 'f_scale': 0.25, 'r_scale': 0.25}
conf['parton_to_hadron'] = True

#--polarized jets
conf['datasets']['pjet'] = {}
conf['datasets']['pjet']['filters'] = []
conf['datasets']['pjet']['filters'].append("pT>8")
conf['datasets']['pjet']['xlsx'] = {}
conf['datasets']['pjet']['xlsx'][20001] = 'pjet/expdata/20001.xlsx' ## STAR 2006 paper on 2003 and 2004 data
conf['datasets']['pjet']['xlsx'][20002] = 'pjet/expdata/20002.xlsx' ## STAR 2012 paper on 2005 data
conf['datasets']['pjet']['xlsx'][20003] = 'pjet/expdata/20003.xlsx' ## STAR 2012 paper on 2006 data
conf['datasets']['pjet']['xlsx'][20004] = 'pjet/expdata/20004.xlsx' ## STAR 2015 paper on 2009 data
conf['datasets']['pjet']['xlsx'][20005] = 'pjet/expdata/20005.xlsx' ## PHENIX 2011 paper on 2005 data
conf['datasets']['pjet']['xlsx'][20006] = 'pjet/expdata/20006.xlsx' ## STAR 2019 paper on 2012 data
conf['datasets']['pjet']['xlsx'][20007] = 'pjet/expdata/20007.xlsx' ## STAR 2021 paper on 2015 data
conf['datasets']['pjet']['xlsx'][20008] = 'pjet/expdata/20008.xlsx' ## STAR 2021 paper on 2013 data
conf['datasets']['pjet']['norm'] = {}
conf['datasets']['pjet']['norm'][20002] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['pjet']['norm'][20003] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['pjet']['norm'][20004] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['pjet']['norm'][20005] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['pjet']['norm'][20006] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['pjet']['norm'][20007] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['datasets']['pjet']['norm'][20008] = {'value': 1, 'fixed': False, 'min': 0.5, 'max': 1.5}
conf['channels'] = {'q_qp': 0, 'q_qbp': 1, 'q_q': 2, 'q_qb': 3, 'q_g': 4, 'g_g': 5}
conf['qr'] = {}
conf['qr']['parameter'] = {'epsilon': 1e-11, 'block_size': 10, 'power_scheme': 1}
conf['pjet_qr_fit'] = {'method': 'fixed', 'f_scale': 0.25, 'r_scale': 0.25}

##--PSIDIS 
conf['datasets']['psidis']={}
conf['datasets']['psidis']['filters']=[]
conf['datasets']['psidis']['filters'].append("Q2>%f"%pQ2cut) 
conf['datasets']['psidis']['filters'].append("W2SIDIS>20") 
conf['datasets']['psidis']['filters'].append('Z>0.2 and Z<0.8')
conf['datasets']['psidis']['xlsx']={}
conf['datasets']['psidis']['xlsx'][20004]='psidis/expdata/20004.xlsx' # 20004 | proton   | A1pi+  | HERMES  
conf['datasets']['psidis']['xlsx'][20005]='psidis/expdata/20005.xlsx' # 20005 | proton   | A1pi-  | HERMES  
conf['datasets']['psidis']['xlsx'][20008]='psidis/expdata/20008.xlsx' # 20008 | deuteron | A1pi+  | HERMES  
conf['datasets']['psidis']['xlsx'][20009]='psidis/expdata/20009.xlsx' # 20009 | deuteron | A1pi-  | HERMES  
conf['datasets']['psidis']['xlsx'][20012]='psidis/expdata/20012.xlsx' # 20012 | deuteron | A1K+   | HERMES  
conf['datasets']['psidis']['xlsx'][20013]='psidis/expdata/20013.xlsx' # 20013 | deuteron | A1K-   | HERMES  
conf['datasets']['psidis']['xlsx'][20014]='psidis/expdata/20014.xlsx' # 20014 | deuteron | A1Ksum | HERMES  
conf['datasets']['psidis']['xlsx'][20017]='psidis/expdata/20017.xlsx' # 20017 | proton   | A1pi+  | COMPASS 
conf['datasets']['psidis']['xlsx'][20018]='psidis/expdata/20018.xlsx' # 20018 | proton   | A1pi-  | COMPASS 
conf['datasets']['psidis']['xlsx'][20019]='psidis/expdata/20019.xlsx' # 20019 | proton   | A1K+   | COMPASS 
conf['datasets']['psidis']['xlsx'][20020]='psidis/expdata/20020.xlsx' # 20020 | proton   | A1K-   | COMPASS 
conf['datasets']['psidis']['xlsx'][20021]='psidis/expdata/20021.xlsx' # 20021 | deuteron | A1pi+  | COMPASS 
conf['datasets']['psidis']['xlsx'][20022]='psidis/expdata/20022.xlsx' # 20022 | deuteron | A1pi-  | COMPASS 
conf['datasets']['psidis']['xlsx'][20025]='psidis/expdata/20025.xlsx' # 20025 | deuteron | A1K+   | COMPASS 
conf['datasets']['psidis']['xlsx'][20026]='psidis/expdata/20026.xlsx' # 20026 | deuteron | A1K-   | COMPASS 
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['psidis']['xlsx'][20000]='psidis/expdata/20000.xlsx' # 20000 | proton   | A1h+   | SMC 
conf['datasets']['psidis']['xlsx'][20001]='psidis/expdata/20001.xlsx' # 20001 | proton   | A1h-   | SMC 
conf['datasets']['psidis']['xlsx'][20002]='psidis/expdata/20002.xlsx' # 20002 | deuteron | A1h+   | SMC 
conf['datasets']['psidis']['xlsx'][20003]='psidis/expdata/20003.xlsx' # 20003 | deuteron | A1h-   | SMC 
conf['datasets']['psidis']['xlsx'][20006]='psidis/expdata/20006.xlsx' # 20006 | proton   | A1h+   | HERMES 
conf['datasets']['psidis']['xlsx'][20007]='psidis/expdata/20007.xlsx' # 20007 | proton   | A1h-   | HERMES 
conf['datasets']['psidis']['xlsx'][20010]='psidis/expdata/20010.xlsx' # 20010 | deuteron | A1h+   | HERMES 
conf['datasets']['psidis']['xlsx'][20011]='psidis/expdata/20011.xlsx' # 20011 | deuteron | A1h-   | HERMES 
conf['datasets']['psidis']['xlsx'][20015]='psidis/expdata/20015.xlsx' # 20015 | helium   | A1h+   | HERMES 
conf['datasets']['psidis']['xlsx'][20016]='psidis/expdata/20016.xlsx' # 20016 | helium   | A1h-   | HERMES 
conf['datasets']['psidis']['xlsx'][20023]='psidis/expdata/20023.xlsx' # 20023 | deuteron | A1h+   | COMPASS 
conf['datasets']['psidis']['xlsx'][20024]='psidis/expdata/20024.xlsx' # 20024 | deuteron | A1h-   | COMPASS 
#------------------------------------------------------------------------------------------------------------------
conf['datasets']['psidis']['norm']={}

#--lepton-lepton reactions

conf['datasets']['sia']={}
conf['datasets']['sia']['filters']={}
#conf['datasets']['sia']['filters'].append('z>0.2 and z<0.9')
zmax  = 0.9
zmin_lowE = 0.15 #--zmin for low energy datasets
zmin_bigE = 0.02 #--zmin for high energy datasets
conf['datasets']['sia']['filters'][1001] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1002] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1003] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1004] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1005] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1006] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1007] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1008] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1010] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1011] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1013] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1014] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1015] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1016] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1017] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1018] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1019] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1023] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1024] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1025] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1026] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1027] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][1028] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1029] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][1030] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['xlsx']={}

#--PION
conf['datasets']['sia']['xlsx'][1001]='sia/expdata_V1/1001freenorm.xlsx'  # TASSO (12 GeV)
conf['datasets']['sia']['xlsx'][1002]='sia/expdata_V1/1002freenorm.xlsx'  # TASSO (14 GeV)
conf['datasets']['sia']['xlsx'][1003]='sia/expdata_V1/1003freenorm.xlsx'  # TASSO (22 GeV)
conf['datasets']['sia']['xlsx'][1004]='sia/expdata_V1/1004freenorm.xlsx'  # TASSO (30 GeV)
conf['datasets']['sia']['xlsx'][1005]='sia/expdata_V1/1005.xlsx'          # TASSO (34 GeV)
conf['datasets']['sia']['xlsx'][1006]='sia/expdata_V1/1006.xlsx'          # TASSO (44 GeV)
conf['datasets']['sia']['xlsx'][1007]='sia/expdata_V1/1007.xlsx'          # TPC 1984
conf['datasets']['sia']['xlsx'][1008]='sia/expdata_V1/1008.xlsx'          # TPC 1988
conf['datasets']['sia']['xlsx'][1010]='sia/expdata_V1/1010.xlsx'          # TPC(c)
conf['datasets']['sia']['xlsx'][1011]='sia/expdata_V1/1011.xlsx'          # TPC(b)
conf['datasets']['sia']['xlsx'][1012]='sia/expdata_V1/1012.xlsx'          # HRS
conf['datasets']['sia']['xlsx'][1013]='sia/expdata_V1/1013.xlsx'          # TOPAZ
conf['datasets']['sia']['xlsx'][1014]='sia/expdata_V1/1014.xlsx'          # SLD
conf['datasets']['sia']['xlsx'][1015]='sia/expdata_V1/1015.xlsx'          # SLD(uds)
conf['datasets']['sia']['xlsx'][1016]='sia/expdata_V1/1016.xlsx'          # SLD(c)
conf['datasets']['sia']['xlsx'][1017]='sia/expdata_V1/1017.xlsx'          # SLD(b)
conf['datasets']['sia']['xlsx'][1018]='sia/expdata_V1/1018.xlsx'          # ALEPH
conf['datasets']['sia']['xlsx'][1019]='sia/expdata_V1/1019.xlsx'          # OPAL
conf['datasets']['sia']['xlsx'][1023]='sia/expdata_V1/1023.xlsx'          # OPAL(c)
conf['datasets']['sia']['xlsx'][1024]='sia/expdata_V1/1024.xlsx'          # OPAL(b)
conf['datasets']['sia']['xlsx'][1025]='sia/expdata_V1/1025.xlsx'          # DELPHI
conf['datasets']['sia']['xlsx'][1026]='sia/expdata_V1/1026.xlsx'          # DELPHI(uds)
conf['datasets']['sia']['xlsx'][1027]='sia/expdata_V1/1027.xlsx'          # DELPHI(b)
conf['datasets']['sia']['xlsx'][1028]='sia/expdata_V1/1028corrmat.xlsx'   # BABAR
conf['datasets']['sia']['xlsx'][1029]='sia/expdata_V1/1029freenorm.xlsx'  # BELL
conf['datasets']['sia']['xlsx'][1030]='sia/expdata_V1/1030.xlsx'          # ARGUS
conf['datasets']['sia']['norm']={}
conf['datasets']['sia']['norm'][1001]={'value':    1.10478e+00,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1002]={'value':    9.82581e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1003]={'value':    1.03054e+00,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1004]={'value':    1.03054e+00,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1005]={'value':    1.03419e+00,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1006]={'value':    9.79162e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1014]={'value':    9.97770e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1015]={'value':    9.97770e-01,'fixed':1014 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1016]={'value':    1.18920e+00,'fixed':1014 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1017]={'value':    1.00345e+00,'fixed':1014 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1018]={'value':    1.02378e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1023]={'value':    1.33434e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1024]={'value':    1.19151e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1028]={'value':    9.76001e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][1029]={'value':    8.68358e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][1030]={'value':    1.01509e+00,'fixed':False,'min':0.5,'max':1.5}


#--KAON 
zmin_babar = 0.2 #--special cut for BABAR kaon
conf['datasets']['sia']['filters'][2001] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2002] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2003] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2004] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2005] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2006] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2007] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2008] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2013] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2014] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2015] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2016] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2017] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2018] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2019] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2023] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2024] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2025] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2026] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2027] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][2028] = ["z < %s and z > %s"%(zmax,zmin_babar)]
conf['datasets']['sia']['filters'][2029] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2030] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][2031] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['xlsx'][2001]='sia/expdata_V1/2001freenorm.xlsx'  # TASSO (12 GeV)
conf['datasets']['sia']['xlsx'][2002]='sia/expdata_V1/2002freenorm.xlsx'  # TASSO (14 GeV)
conf['datasets']['sia']['xlsx'][2003]='sia/expdata_V1/2003freenorm.xlsx'  # TASSO (22 GeV)
conf['datasets']['sia']['xlsx'][2004]='sia/expdata_V1/2004freenorm.xlsx'  # TASSO (30 GeV)
conf['datasets']['sia']['xlsx'][2005]='sia/expdata_V1/2005.xlsx'          # TASSO (34 GeV)
conf['datasets']['sia']['xlsx'][2006]='sia/expdata_V1/2006.xlsx'          # TASSO (44 Gev)
conf['datasets']['sia']['xlsx'][2007]='sia/expdata_V1/2007.xlsx'          # TPC 1984
conf['datasets']['sia']['xlsx'][2008]='sia/expdata_V1/2008.xlsx'          # TPC 1988
conf['datasets']['sia']['xlsx'][2012]='sia/expdata_V1/2012.xlsx'          # HRS
conf['datasets']['sia']['xlsx'][2013]='sia/expdata_V1/2013.xlsx'          # TOPAZ
conf['datasets']['sia']['xlsx'][2014]='sia/expdata_V1/2014.xlsx'          # SLD
conf['datasets']['sia']['xlsx'][2015]='sia/expdata_V1/2015.xlsx'          # SLD(uds)
conf['datasets']['sia']['xlsx'][2016]='sia/expdata_V1/2016.xlsx'          # SLD(c)
conf['datasets']['sia']['xlsx'][2017]='sia/expdata_V1/2017.xlsx'          # SLD(b)
conf['datasets']['sia']['xlsx'][2018]='sia/expdata_V1/2018.xlsx'          # ALEPH
conf['datasets']['sia']['xlsx'][2019]='sia/expdata_V1/2019.xlsx'          # OPAL
conf['datasets']['sia']['xlsx'][2023]='sia/expdata_V1/2023.xlsx'          # OPAL(c)
conf['datasets']['sia']['xlsx'][2024]='sia/expdata_V1/2024.xlsx'          # OPAL(b)
conf['datasets']['sia']['xlsx'][2025]='sia/expdata_V1/2025.xlsx'          # DELPHI
conf['datasets']['sia']['xlsx'][2026]='sia/expdata_V1/2026.xlsx'          # DELPHI(uds)
conf['datasets']['sia']['xlsx'][2027]='sia/expdata_V1/2027.xlsx'          # DELPHI(b)
conf['datasets']['sia']['xlsx'][2028]='sia/expdata_V1/2028corrmat.xlsx'   # BABAR
conf['datasets']['sia']['xlsx'][2029]='sia/expdata_V1/2029freenorm.xlsx'  # BELL
conf['datasets']['sia']['xlsx'][2030]='sia/expdata_V1/2030.xlsx'          # ARGUS
conf['datasets']['sia']['xlsx'][2031]='sia/expdata_V1/2031.xlsx'          # DELPHI 1995
conf['datasets']['sia']['norm'][2001]={'value':    9.83394e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][2002]={'value':    9.83394e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][2003]={'value':    9.94421e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][2004]={'value':    9.94421e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][2005]={'value':    9.92876e-01,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][2014]={'value':    9.12186e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2015]={'value':    9.12186e-01,'fixed':2014 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2016]={'value':    1.03996e+00,'fixed':2014 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2017]={'value':    1.00015e+00,'fixed':2014 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2018]={'value':    8.88453e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2023]={'value':    1.35922e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2024]={'value':    1.38163e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2028]={'value':    9.97435e-01,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][2029]={'value':    1.00000e+00,'fixed':False,'min':0.1,'max':2.0}
conf['datasets']['sia']['norm'][2030]={'value':    1.00482e+00,'fixed':False,'min':0.5,'max':1.5}

#--HADRONS
conf['datasets']['sia']['filters'][4000] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4001] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4002] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4003] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4004] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4005] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4006] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4007] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4008] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4009] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4010] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4011] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4012] = ["z < %s and z > %s"%(zmax,zmin_lowE)]
conf['datasets']['sia']['filters'][4013] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4014] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4015] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4016] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['filters'][4017] = ["z < %s and z > %s"%(zmax,zmin_bigE)]
conf['datasets']['sia']['xlsx'][4000]='sia/expdata_V1/4000.xlsx'     # ALEPH
conf['datasets']['sia']['xlsx'][4001]='sia/expdata_V1/4001.xlsx'     # DELPHI
conf['datasets']['sia']['xlsx'][4002]='sia/expdata_V1/4002.xlsx'     # SLD
conf['datasets']['sia']['xlsx'][4003]='sia/expdata_V1/4003.xlsx'     # TASSO (12 GeV)
conf['datasets']['sia']['xlsx'][4004]='sia/expdata_V1/4004.xlsx'     # TPC
conf['datasets']['sia']['xlsx'][4005]='sia/expdata_V1/4005.xlsx'     # OPAL(b)
conf['datasets']['sia']['xlsx'][4006]='sia/expdata_V1/4006.xlsx'     # OPAL(c)
conf['datasets']['sia']['xlsx'][4007]='sia/expdata_V1/4007.xlsx'     # OPAL
conf['datasets']['sia']['xlsx'][4008]='sia/expdata_V1/4008.xlsx'     # TASSO (14 GeV)
conf['datasets']['sia']['xlsx'][4009]='sia/expdata_V1/4009.xlsx'     # TASSO (22 GeV)
conf['datasets']['sia']['xlsx'][4010]='sia/expdata_V1/4010.xlsx'     # TASSO (30 GeV)
conf['datasets']['sia']['xlsx'][4011]='sia/expdata_V1/4011.xlsx'     # TASSO (35 GeV)
conf['datasets']['sia']['xlsx'][4012]='sia/expdata_V1/4012.xlsx'     # TASSO (44 GeV)
conf['datasets']['sia']['xlsx'][4013]='sia/expdata_V1/4013.xlsx'     # DELPHI(b)
conf['datasets']['sia']['xlsx'][4014]='sia/expdata_V1/4014.xlsx'     # SLD(c)
conf['datasets']['sia']['xlsx'][4015]='sia/expdata_V1/4015.xlsx'     # SLD(b)
conf['datasets']['sia']['xlsx'][4016]='sia/expdata_V1/4016.xlsx'     # SLD(uds)
conf['datasets']['sia']['xlsx'][4017]='sia/expdata_V1/4017.xlsx'     # DELPHI(uds)
conf['datasets']['sia']['norm'][4000]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4002]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4003]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4008]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4009]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4010]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4011]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4012]={'value':    1.00000e+00,'fixed':False,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4014]={'value':    1.00000e+00,'fixed':4002 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4015]={'value':    1.00000e+00,'fixed':4002 ,'min':0.5,'max':1.5}
conf['datasets']['sia']['norm'][4016]={'value':    1.00000e+00,'fixed':4002 ,'min':0.5,'max':1.5}



##--SIDIS 
conf['datasets']['sidis']={}
conf['datasets']['sidis']['filters']=[]
conf['datasets']['sidis']['filters'].append("Q2>%f"%Q2cut) 
conf['datasets']['sidis']['filters'].append("W2SIDIS>20") 
conf['datasets']['sidis']['filters'].append('Z>0.2 and Z<0.8')
conf['datasets']['sidis']['xlsx']={}
conf['datasets']['sidis']['xlsx'][1005]='sidis/expdata/1005.xlsx' # deuteron , mult , pi+ , COMPASS
conf['datasets']['sidis']['xlsx'][1006]='sidis/expdata/1006.xlsx' # deuteron , mult , pi- , COMPASS
conf['datasets']['sidis']['xlsx'][2005]='sidis/expdata/2005.xlsx' # deuteron , mult , K+  , COMPASS
conf['datasets']['sidis']['xlsx'][2006]='sidis/expdata/2006.xlsx' # deuteron , mult , K-  , COMPASS
conf['datasets']['sidis']['xlsx'][3000]='sidis/expdata/3000.xlsx' # deuteron , mult , h+  , COMPASS
conf['datasets']['sidis']['xlsx'][3001]='sidis/expdata/3001.xlsx' # deuteron , mult , h-  , COMPASS
conf['datasets']['sidis']['norm']={}

#--parameters
conf['params']={}

#--pdf parameters
conf['params']['pdf']={}

conf['params']['pdf']['g1 N']    ={'value':    3.87592e-01   , 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['g1 a']    ={'value':   -6.23068169e-01, 'min':  -1.9, 'max':     1, 'fixed': False}
conf['params']['pdf']['g1 b']    ={'value':    9.25741583e+00, 'min':     0, 'max':    20, 'fixed': False}
conf['params']['pdf']['g1 c']    ={'value':    0.00000000e+00, 'min':  -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['g1 d']    ={'value':    0.00000000e+00, 'min':  -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['uv1 N']   ={'value':    3.47549e-01   , 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['uv1 a']   ={'value':   -1.21835956e-01, 'min':  -0.9, 'max':     1, 'fixed': False}
conf['params']['pdf']['uv1 b']   ={'value':    3.20766744e+00, 'min':     0, 'max':    10, 'fixed': False}
conf['params']['pdf']['uv1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['uv1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['dv1 N']   ={'value':    1.52089e-01   , 'min':  None, 'max':  None, 'fixed': True }
conf['params']['pdf']['dv1 a']   ={'value':   -2.39874967e-01, 'min':  -0.9, 'max':     1, 'fixed': False}
conf['params']['pdf']['dv1 b']   ={'value':    3.83902620e+00, 'min':     0, 'max':    10, 'fixed': False}
conf['params']['pdf']['dv1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['dv1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['mix N']    ={'value':    0.00000000e+00   , 'min':  -0.1, 'max': 5.0, 'fixed': False}
conf['params']['pdf']['mix a']    ={'value':    2.00000000e+00   , 'min':   1.0, 'max':10.0, 'fixed': False}

conf['params']['pdf']['db1 N']   ={'value':    3.67609928e-02, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['db1 a']   ={'value':   -8.41360631e-01, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['pdf']['db1 b']   ={'value':    5.31285539e+00, 'min':     0, 'max':    50, 'fixed': False}
conf['params']['pdf']['db1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['db1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['ub1 N']   ={'value':    1.95464789e-02, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['ub1 a']   ={'value':   -9.93659187e-01, 'min':    -1, 'max':     1, 'fixed': False}
conf['params']['pdf']['ub1 b']   ={'value':    8.38905814e+00, 'min':     0, 'max':    20, 'fixed': False}
conf['params']['pdf']['ub1 c']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}
conf['params']['pdf']['ub1 d']   ={'value':    0.00000000e+00, 'min':   -100, 'max':    100, 'fixed': False}

conf['params']['pdf']['s1 N']    ={'value':       0.00000e+00, 'min':     0, 'max':     1, 'fixed': True }
conf['params']['pdf']['s1 a']    ={'value':    1.34706224e-01, 'min':  -0.5, 'max':     1, 'fixed': False}
conf['params']['pdf']['s1 b']    ={'value':    6.00759596e+00, 'min':     0, 'max':   100, 'fixed': False}

conf['params']['pdf']['sb1 N']   ={'value':    7.46109845e-07, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['sb1 a']   ={'value':    3.83495317e-01, 'min':  -0.9, 'max':     1, 'fixed': False}
conf['params']['pdf']['sb1 b']   ={'value':    4.61209808e+00, 'min':     0, 'max':    10, 'fixed': False}

conf['params']['pdf']['sea1 N']  ={'value':    5.71081196e-03, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['sea1 a']  ={'value':   -1.36329697e+00, 'min':  -1.9, 'max':    -1, 'fixed': False}
conf['params']['pdf']['sea1 b']  ={'value':    4.74721050e+00, 'min':     0, 'max':   100, 'fixed': False}

conf['params']['pdf']['sea2 N']  ={'value':       2.08640e-02, 'min':     0, 'max':     1, 'fixed': False}
conf['params']['pdf']['sea2 a']  ={'value':      -1.500000000, 'min':  -1.9, 'max':    -1, 'fixed': False}
conf['params']['pdf']['sea2 b']  ={'value':       1.00000e+01, 'min':     0, 'max':   100, 'fixed': False}

#--ppdf parameters
conf['ppdf_choice'] = 'valence'
conf['params']['ppdf'] = {}

conf['params']['ppdf']['g1 N']    = {'value':    1.54709e+00, 'min':   -10, 'max':    30, 'fixed': False}
conf['params']['ppdf']['g1 a']    = {'value':   -7.13518e-01, 'min': -0.99, 'max':     5, 'fixed': False}
conf['params']['ppdf']['g1 b']    = {'value':    9.45902e-01, 'min':     0, 'max':    30, 'fixed': False}
conf['params']['ppdf']['g1 d']    = {'value':    0.00000e+00, 'min': -1000, 'max':    10, 'fixed': False}

conf['params']['ppdf']['uv1 N']   = {'value':    3.47549e-01, 'min':   -10, 'max':    10, 'fixed': False}
conf['params']['ppdf']['uv1 a']   = {'value':    2.92830e-01, 'min': -0.99, 'max':     5, 'fixed': False}
conf['params']['ppdf']['uv1 b']   = {'value':    2.95841e+00, 'min':     0, 'max':    30, 'fixed': False}
conf['params']['ppdf']['uv1 d']   = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': False}

conf['params']['ppdf']['dv1 N']   = {'value':    1.52089e-01, 'min':   -10, 'max':    10, 'fixed': False}
conf['params']['ppdf']['dv1 a']   = {'value':    5.34011e-01, 'min': -0.99, 'max':     5, 'fixed': False}
conf['params']['ppdf']['dv1 b']   = {'value':    2.01125e+00, 'min':     0, 'max':    30, 'fixed': False}
conf['params']['ppdf']['dv1 d']   = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': False}

conf['params']['ppdf']['ub1 N']   = {'value':   -1.21705e-01, 'min':   -10, 'max':    10, 'fixed': False}
conf['params']['ppdf']['ub1 a']   = {'value':    1.36941e+00, 'min': -0.99, 'max':     5, 'fixed': False}
conf['params']['ppdf']['ub1 b']   = {'value':    9.99976e+00, 'min':     0, 'max':    30, 'fixed': False}
conf['params']['ppdf']['ub1 d']   = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': False}

conf['params']['ppdf']['db1 N']   = {'value':   -1.21705e-01, 'min':   -10, 'max':    10, 'fixed': False}
conf['params']['ppdf']['db1 a']   = {'value':    1.36941e+00, 'min': -0.99, 'max':     5, 'fixed': False}
conf['params']['ppdf']['db1 b']   = {'value':    9.99976e+00, 'min':     0, 'max':    30, 'fixed': False}
conf['params']['ppdf']['db1 d']   = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': False}

conf['params']['ppdf']['s1 N']    = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf']['s1 a']    = {'value':    0.000000000, 'min': -0.99, 'max':     5, 'fixed': True}
conf['params']['ppdf']['s1 b']    = {'value':    0.000000000, 'min':     0, 'max':    30, 'fixed': True}
conf['params']['ppdf']['s1 d']    = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': True}

conf['params']['ppdf']['sb1 N']   = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': True}
conf['params']['ppdf']['sb1 a']   = {'value':    0.000000000, 'min': -0.99, 'max':     5, 'fixed': True}
conf['params']['ppdf']['sb1 b']   = {'value':    0.000000000, 'min':     0, 'max':    30, 'fixed': True}
conf['params']['ppdf']['sb1 d']   = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': True}

conf['params']['ppdf']['sea1 N']  = {'value':   -1.21705e-01, 'min':   -10, 'max':    10, 'fixed': False}
conf['params']['ppdf']['sea1 a']  = {'value':    1.36941e+00, 'min': -0.99, 'max':     5, 'fixed': False}
conf['params']['ppdf']['sea1 b']  = {'value':    9.99976e+00, 'min':     0, 'max':    30, 'fixed': False}
conf['params']['ppdf']['sea1 d']  = {'value':    0.000000000, 'min':   -10, 'max':    10, 'fixed': False}

# gXres params
conf['params']['gXres'] = {}

conf['params']['gXres']['g1p N0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1p a0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1p b0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1p c0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1p d0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}

conf['params']['gXres']['g2p N0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2p a0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2p b0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2p c0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2p d0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}

conf['params']['gXres']['g1n N0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1n a0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1n b0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1n c0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g1n d0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}

conf['params']['gXres']['g2n N0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2n a0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2n b0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2n c0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}
conf['params']['gXres']['g2n d0'] = {'value':   0, 'min':   -1.5, 'max':    1.5, 'fixed': False}



#--ht params (NOTE: What's being fitted here is the (M^2/Q^2) and (m^2/Q^2) correction H. To get the correction to FX, take H/Q^2)

conf['params']['ht4']={}

conf['params']['ht4']['F2p N']  ={'value':   0.0,   'min': -30.0,  'max': 30.0, 'fixed':False}
conf['params']['ht4']['F2p a']  ={'value':   1.5,   'min':  0.0,   'max': 4,    'fixed':False}
conf['params']['ht4']['F2p b']  ={'value':   0.0,   'min':  None,  'max': None, 'fixed':True}
conf['params']['ht4']['F2p c']  ={'value':   0.0,   'min': -10.0,  'max': 10,   'fixed':True}
conf['params']['ht4']['F2p d']  ={'value':  -5.0,   'min': -10.0,  'max': -1.0, 'fixed':False}

conf['params']['ht4']['FLp N']  ={'value':   0.0,   'min': -30.0,  'max': 30.0, 'fixed':'F2p N'}
conf['params']['ht4']['FLp a']  ={'value':   1.5,   'min':  0.0,   'max': 4,    'fixed':'F2p a'}
conf['params']['ht4']['FLp b']  ={'value':   0.0,   'min':  None,  'max': None, 'fixed':'F2p b'}
conf['params']['ht4']['FLp c']  ={'value':   0.0,   'min': -10.0,  'max': 10,   'fixed':'F2p c'}
conf['params']['ht4']['FLp d']  ={'value':  -5.0,   'min': -10.0,  'max': -1.0, 'fixed':'F2p d'}

conf['params']['ht4']['F3p N']  ={'value':   0.0,   'min': -30.0,  'max': 30.0, 'fixed':'F2p N'}
conf['params']['ht4']['F3p a']  ={'value':   1.5,   'min':  0.0,   'max': 4,    'fixed':'F2p a'}
conf['params']['ht4']['F3p b']  ={'value':   0.0,   'min':  None,  'max': None, 'fixed':'F2p b'}
conf['params']['ht4']['F3p c']  ={'value':   0.0,   'min': -10.0,  'max': 10,   'fixed':'F2p c'}
conf['params']['ht4']['F3p d']  ={'value':  -5.0,   'min': -10.0,  'max': -1.0, 'fixed':'F2p d'}

#--for p != n
conf['params']['ht4']['F2n N']  ={'value':   0.0,   'min': -30.0,  'max': 30.0, 'fixed':False}
conf['params']['ht4']['F2n a']  ={'value':   1.5,   'min':  0.0,   'max': 4,    'fixed':False}
conf['params']['ht4']['F2n b']  ={'value':   0.0,   'min':  None,  'max': None, 'fixed':True}
conf['params']['ht4']['F2n c']  ={'value':   0.0,   'min': -10.0,  'max': 10,   'fixed':True}
conf['params']['ht4']['F2n d']  ={'value':  -5.0,   'min':  -10.0, 'max': -1.0, 'fixed':False}

conf['params']['ht4']['FLn N']  ={'value':   0.0,   'min': -30.0,  'max': 30.0, 'fixed':'F2n N'}
conf['params']['ht4']['FLn a']  ={'value':   1.5,   'min':  0.0,   'max': 4,    'fixed':'F2n a'}
conf['params']['ht4']['FLn b']  ={'value':   0.0,   'min':  None,  'max': None, 'fixed':'F2n b'}
conf['params']['ht4']['FLn c']  ={'value':   0.0,   'min': -10.0,  'max': 10,   'fixed':'F2n c'}
conf['params']['ht4']['FLn d']  ={'value':  -5.0,   'min':  -10.0, 'max': -1.0, 'fixed':'F2n d'}

#--offshell paramaterization
conf['params']['off pdf']={}

conf['params']['off pdf']['uv1 N']   ={'value':   0.0, 'min': -10.0, 'max':  10.0, 'fixed': False}
conf['params']['off pdf']['uv1 a']   ={'value':   0.0, 'min':  -0.9, 'max':     5, 'fixed': False}
conf['params']['off pdf']['uv1 b']   ={'value':   3.0, 'min':   0.0, 'max':    20, 'fixed': False}
conf['params']['off pdf']['uv1 c']   ={'value':   0.0, 'min': -10.0, 'max':  10.0, 'fixed': True}
conf['params']['off pdf']['uv1 d']   ={'value':   0.0, 'min':  None, 'max':  None, 'fixed': True}

conf['params']['off pdf']['dv1 N']   ={'value':   0.0, 'min': -10.0, 'max':  10.0, 'fixed': False}
conf['params']['off pdf']['dv1 a']   ={'value':   0.0, 'min':  -0.9, 'max':     5, 'fixed': False}
conf['params']['off pdf']['dv1 b']   ={'value':   3.0, 'min':   0.0, 'max':    20, 'fixed': False}
conf['params']['off pdf']['dv1 c']   ={'value':   0.0, 'min': -10.0, 'max':  10.0, 'fixed': True}
conf['params']['off pdf']['dv1 d']   ={'value':   0.0, 'min':  None, 'max':  None, 'fixed': True}


#--pion fragmentation
conf['ffpion parametrization']=0
conf['params']['ffpion']={}

conf['params']['ffpion']['g1 N']  ={'value':    2.95437e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['g1 a']  ={'value':    1.00469e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['g1 b']  ={'value':    6.85766e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['g1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['g1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                
conf['params']['ffpion']['u1 N']  ={'value':    2.67821e-02, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['u1 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['u1 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['u1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['u1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffpion']['d1 N']  ={'value':    2.99974e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['d1 a']  ={'value':   -6.89477e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['d1 b']  ={'value':    4.79992e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['d1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['d1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffpion']['s1 N']  ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['s1 a']  ={'value':    3.00305e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['s1 b']  ={'value':    1.83178e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['s1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['s1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                  
conf['params']['ffpion']['c1 N']  ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['c1 a']  ={'value':   -5.05798e-02, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['c1 b']  ={'value':    3.19952e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['c1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['c1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffpion']['b1 N']  ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['b1 a']  ={'value':   -1.59541e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['b1 b']  ={'value':    4.50102e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['b1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['b1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffpion']['ub1 N'] ={'value':    2.99974e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['ffpion']['ub1 a'] ={'value':   -6.89477e-01, 'min': -1.8 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffpion']['ub1 b'] ={'value':    4.79992e+00, 'min':  0   , 'max':    30,'fixed':'d1 b'}
conf['params']['ffpion']['ub1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 c'}
conf['params']['ffpion']['ub1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 d'}
                                                                                                  
conf['params']['ffpion']['db1 N'] ={'value':    2.67821e-02, 'min':  0.0 , 'max':     1,'fixed':'u1 N'}
conf['params']['ffpion']['db1 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'u1 a'}
conf['params']['ffpion']['db1 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'u1 b'}
conf['params']['ffpion']['db1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'u1 c'}
conf['params']['ffpion']['db1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'u1 d'}
                                                                                                  
conf['params']['ffpion']['sb1 N'] ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':'s1 N'}
conf['params']['ffpion']['sb1 a'] ={'value':    3.00305e-01, 'min': -1.8 , 'max':    10,'fixed':'s1 a'}
conf['params']['ffpion']['sb1 b'] ={'value':    1.83178e+00, 'min':  0   , 'max':    30,'fixed':'s1 b'}
conf['params']['ffpion']['sb1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s1 c'}
conf['params']['ffpion']['sb1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s1 d'}
                                                                                                  
conf['params']['ffpion']['cb1 N'] ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':'c1 N'}
conf['params']['ffpion']['cb1 a'] ={'value':   -5.05798e-02, 'min': -1.8 , 'max':    10,'fixed':'c1 a'}
conf['params']['ffpion']['cb1 b'] ={'value':    3.19952e+00, 'min':  0   , 'max':    30,'fixed':'c1 b'}
conf['params']['ffpion']['cb1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c1 c'}
conf['params']['ffpion']['cb1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c1 d'}
                                                                                                  
conf['params']['ffpion']['bb1 N'] ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':'b1 N'}
conf['params']['ffpion']['bb1 a'] ={'value':   -1.59541e+00, 'min': -1.8 , 'max':    10,'fixed':'b1 a'}
conf['params']['ffpion']['bb1 b'] ={'value':    4.50102e+00, 'min':  0   , 'max':    30,'fixed':'b1 b'}
conf['params']['ffpion']['bb1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b1 c'}
conf['params']['ffpion']['bb1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b1 d'}
                          
        
conf['params']['ffpion']['g2 N']  ={'value':    2.95437e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['g2 a']  ={'value':    1.00469e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['g2 b']  ={'value':    6.85766e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['g2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['g2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
 
conf['params']['ffpion']['u2 N']  ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['u2 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['u2 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['u2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['u2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                
conf['params']['ffpion']['d2 N']  ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['d2 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['d2 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['d2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['d2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffpion']['ub2 N'] ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':'d2 N'}
conf['params']['ffpion']['ub2 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffpion']['ub2 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'d2 b'}
conf['params']['ffpion']['ub2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d2 c'}
conf['params']['ffpion']['ub2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d2 d'}
                                                                       
conf['params']['ffpion']['db2 N'] ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':'u2 N'}
conf['params']['ffpion']['db2 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'u2 a'}
conf['params']['ffpion']['db2 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'u2 b'}
conf['params']['ffpion']['db2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'u2 c'}
conf['params']['ffpion']['db2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'u2 d'}
                                                                                           

conf['params']['ffpion']['s2 N']  ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['s2 a']  ={'value':    3.00305e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['s2 b']  ={'value':    1.83178e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['s2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['s2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                               
conf['params']['ffpion']['c2 N']  ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['c2 a']  ={'value':   -5.05798e-02, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['c2 b']  ={'value':    3.19952e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['c2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['c2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                 
conf['params']['ffpion']['b2 N']  ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffpion']['b2 a']  ={'value':   -1.59541e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffpion']['b2 b']  ={'value':    4.50102e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffpion']['b2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffpion']['b2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}

conf['params']['ffpion']['sb2 N'] ={'value':    1.54863e-01, 'min':  0.0 , 'max':     1,'fixed':'s2 N'}
conf['params']['ffpion']['sb2 a'] ={'value':    3.00305e-01, 'min': -1.8 , 'max':    10,'fixed':'s2 a'}
conf['params']['ffpion']['sb2 b'] ={'value':    1.83178e+00, 'min':  0   , 'max':    30,'fixed':'s2 b'}
conf['params']['ffpion']['sb2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s2 c'}
conf['params']['ffpion']['sb2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s2 d'}

conf['params']['ffpion']['cb2 N'] ={'value':    1.84550e-01, 'min':  0.0 , 'max':     1,'fixed':'c2 N'}
conf['params']['ffpion']['cb2 a'] ={'value':   -5.05798e-02, 'min': -1.8 , 'max':    10,'fixed':'c2 a'}
conf['params']['ffpion']['cb2 b'] ={'value':    3.19952e+00, 'min':  0   , 'max':    30,'fixed':'c2 b'}
conf['params']['ffpion']['cb2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c2 c'}
conf['params']['ffpion']['cb2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c2 d'}

conf['params']['ffpion']['bb2 N'] ={'value':    3.74125e-01, 'min':  0.0 , 'max':     1,'fixed':'b2 N'}
conf['params']['ffpion']['bb2 a'] ={'value':   -1.59541e+00, 'min': -1.8 , 'max':    10,'fixed':'b2 a'}
conf['params']['ffpion']['bb2 b'] ={'value':    4.50102e+00, 'min':  0   , 'max':    30,'fixed':'b2 b'}
conf['params']['ffpion']['bb2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b2 c'}
conf['params']['ffpion']['bb2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b2 d'}
                                                                                           

#--kaon fragmentation
conf['params']['ffkaon']={}
conf['ffkaon parametrization']=0

conf['params']['ffkaon']['g1 N']  ={'value':    2.33320e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['g1 a']  ={'value':    1.48737e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['g1 b']  ={'value':    9.62755e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['g1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['g1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                
conf['params']['ffkaon']['u1 N']  ={'value':    4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['u1 a']  ={'value':    1.26356e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['u1 b']  ={'value':    1.62596e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['u1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['u1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffkaon']['d1 N']  ={'value':    2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['d1 a']  ={'value':   -1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['d1 b']  ={'value':    5.65143e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['d1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['d1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffkaon']['s1 N']  ={'value':    2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['ffkaon']['s1 a']  ={'value':  - 1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffkaon']['s1 b']  ={'value':    5.65143e+00, 'min':  0   , 'max':    30,'fixed':'d1 b'}
conf['params']['ffkaon']['s1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 c'}
conf['params']['ffkaon']['s1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 d'}
                                                                                                   
conf['params']['ffkaon']['c1 N']  ={'value':    7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['c1 a']  ={'value':   -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['c1 b']  ={'value':    2.50452e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['c1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['c1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffkaon']['b1 N']  ={'value':    5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['b1 a']  ={'value':   -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['b1 b']  ={'value':    3.41727e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['b1 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['b1 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffkaon']['ub1 N'] ={'value':    4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['ffkaon']['ub1 a'] ={'value':    2.26356e+00, 'min': -1.8 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffkaon']['ub1 b'] ={'value':    1.62596e+00, 'min':  0   , 'max':    30,'fixed':'d1 b'}
conf['params']['ffkaon']['ub1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 c'}
conf['params']['ffkaon']['ub1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 d'}
                                                                                               
conf['params']['ffkaon']['db1 N'] ={'value':    2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':'d1 N'}
conf['params']['ffkaon']['db1 a'] ={'value':   -1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':'d1 a'}
conf['params']['ffkaon']['db1 b'] ={'value':    1.65143e+00, 'min':  0   , 'max':    30,'fixed':'d1 b'}
conf['params']['ffkaon']['db1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 c'}
conf['params']['ffkaon']['db1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d1 d'}
                                                                                               
conf['params']['ffkaon']['sb1 N'] ={'value':    2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['sb1 a'] ={'value':   -1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['sb1 b'] ={'value':    1.65143e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['sb1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['sb1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                  
conf['params']['ffkaon']['cb1 N'] ={'value':    7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':'c1 N'}
conf['params']['ffkaon']['cb1 a'] ={'value':   -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':'c1 a'}
conf['params']['ffkaon']['cb1 b'] ={'value':    2.50452e+00, 'min':  0   , 'max':    30,'fixed':'c1 b'}
conf['params']['ffkaon']['cb1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c1 c'}
conf['params']['ffkaon']['cb1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c1 d'}
                                                                                                  
conf['params']['ffkaon']['bb1 N'] ={'value':    5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':'b1 N'}
conf['params']['ffkaon']['bb1 a'] ={'value':   -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':'b1 a'}
conf['params']['ffkaon']['bb1 b'] ={'value':    3.41727e+00, 'min':  0   , 'max':    30,'fixed':'b1 b'}
conf['params']['ffkaon']['bb1 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b1 c'}
conf['params']['ffkaon']['bb1 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b1 d'}


conf['params']['ffkaon']['g2 N']  ={'value':    2.33320e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['g2 a']  ={'value':    1.48737e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['g2 b']  ={'value':    9.62755e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['g2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['g2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}

conf['params']['ffkaon']['u2 N']  ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['u2 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['u2 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['u2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['u2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                              
conf['params']['ffkaon']['d2 N']  ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['d2 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['d2 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['d2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['d2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}

conf['params']['ffkaon']['s2 N']  ={'value':    2.51671e-01, 'min':   0.0 , 'max':     1,'fixed':'d2 N'}
conf['params']['ffkaon']['s2 a']  ={'value':   -1.43444e+00, 'min':  -1.8 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffkaon']['s2 b']  ={'value':    5.65143e+00, 'min':   0   , 'max':    30,'fixed':'d2 b'}
conf['params']['ffkaon']['s2 c']  ={'value':    4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'d2 c'}
conf['params']['ffkaon']['s2 d']  ={'value':    4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'d2 d'}
                                                                                                   
conf['params']['ffkaon']['c2 N']  ={'value':    7.76923e-01, 'min':   0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['c2 a']  ={'value':   -1.80000e+00, 'min':  -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['c2 b']  ={'value':    2.50452e+00, 'min':   0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['c2 c']  ={'value':    4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['c2 d']  ={'value':    4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffkaon']['b2 N']  ={'value':    5.66971e-01, 'min':   0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['b2 a']  ={'value':   -1.80000e+00, 'min':  -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['b2 b']  ={'value':    3.41727e+00, 'min':   0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['b2 c']  ={'value':    4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['b2 d']  ={'value':    4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
                                                                                                 
conf['params']['ffkaon']['ub2 N'] ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':'d2 N'}
conf['params']['ffkaon']['ub2 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffkaon']['ub2 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'d2 b'}
conf['params']['ffkaon']['ub2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d2 c'}
conf['params']['ffkaon']['ub2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d2 d'}
                                                                     
conf['params']['ffkaon']['db2 N'] ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':'d2 N'}
conf['params']['ffkaon']['db2 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'d2 a'}
conf['params']['ffkaon']['db2 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'d2 b'}
conf['params']['ffkaon']['db2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d2 c'}
conf['params']['ffkaon']['db2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'d2 d'}

conf['params']['ffkaon']['sb2 N'] ={'value':    2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffkaon']['sb2 a'] ={'value':   -1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffkaon']['sb2 b'] ={'value':    1.65143e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffkaon']['sb2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffkaon']['sb2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                  
conf['params']['ffkaon']['cb2 N'] ={'value':    7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':'c2 N'}
conf['params']['ffkaon']['cb2 a'] ={'value':   -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':'c2 a'}
conf['params']['ffkaon']['cb2 b'] ={'value':    2.50452e+00, 'min':  0   , 'max':    30,'fixed':'c2 b'}
conf['params']['ffkaon']['cb2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c2 c'}
conf['params']['ffkaon']['cb2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c2 d'}
                                                                                                  
conf['params']['ffkaon']['bb2 N'] ={'value':    5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':'b2 N'}
conf['params']['ffkaon']['bb2 a'] ={'value':   -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':'b2 a'}
conf['params']['ffkaon']['bb2 b'] ={'value':    3.41727e+00, 'min':  0   , 'max':    30,'fixed':'b2 b'}
conf['params']['ffkaon']['bb2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b2 c'}
conf['params']['ffkaon']['bb2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b2 d'}

#--unidentified hadrons
conf['params']['ffhadron']={}
conf['ffhadron parametrization']='sum2'

conf['params']['ffhadron']['g1 N']  ={'value':   2.33320e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['g1 a']  ={'value':   1.48737e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['g1 b']  ={'value':   9.62755e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['g1 c']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['g1 d']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                
conf['params']['ffhadron']['u1 N']  ={'value':   4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['u1 a']  ={'value':   1.26356e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['u1 b']  ={'value':   1.62596e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['u1 c']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['u1 d']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['d1 N']  ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['d1 a']  ={'value':  -1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['d1 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['d1 c']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['d1 d']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['s1 N']  ={'value':   2.51671e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['s1 a']  ={'value':  -1.43444e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['s1 b']  ={'value':   5.65143e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['s1 c']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['s1 d']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['c1 N']  ={'value':   7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['c1 a']  ={'value':  -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['c1 b']  ={'value':   2.50452e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['c1 c']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['c1 d']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['b1 N']  ={'value':   5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['b1 a']  ={'value':  -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['b1 b']  ={'value':   3.41727e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['b1 c']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['b1 d']  ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['ub1 N'] ={'value':   4.03672e-02, 'min':  0.0 , 'max':     1,'fixed':'s1 N'}
conf['params']['ffhadron']['ub1 a'] ={'value':   2.26356e+00, 'min': -1.8 , 'max':    10,'fixed':'s1 a'}
conf['params']['ffhadron']['ub1 b'] ={'value':   1.62596e+00, 'min':  0   , 'max':    30,'fixed':'s1 b'}
conf['params']['ffhadron']['ub1 c'] ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s1 c'}
conf['params']['ffhadron']['ub1 d'] ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s1 d'}
                                                                                               
conf['params']['ffhadron']['db1 N'] ={'value':   2.51671e-01, 'min':   0.0 , 'max':     1,'fixed':'s1 N'}
conf['params']['ffhadron']['db1 a'] ={'value':  -1.43444e+00, 'min':  -1.8 , 'max':    10,'fixed':'s1 a'}
conf['params']['ffhadron']['db1 b'] ={'value':   1.65143e+00, 'min':   0   , 'max':    30,'fixed':'s1 b'}
conf['params']['ffhadron']['db1 c'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'s1 c'}
conf['params']['ffhadron']['db1 d'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'s1 d'}
                                                                                               
conf['params']['ffhadron']['sb1 N'] ={'value':   2.51671e-01, 'min':   0.0 , 'max':     1,'fixed':'s1 N'}
conf['params']['ffhadron']['sb1 a'] ={'value':  -1.43444e+00, 'min':  -1.8 , 'max':    10,'fixed':'s1 a'}
conf['params']['ffhadron']['sb1 b'] ={'value':   1.65143e+00, 'min':   0   , 'max':    30,'fixed':'s1 b'}
conf['params']['ffhadron']['sb1 c'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'s1 c'}
conf['params']['ffhadron']['sb1 d'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'s1 d'}
                                                                                                  
conf['params']['ffhadron']['cb1 N'] ={'value':   7.76923e-01, 'min':   0.0 , 'max':     1,'fixed':'c1 N'}
conf['params']['ffhadron']['cb1 a'] ={'value':  -1.80000e+00, 'min':  -1.8 , 'max':    10,'fixed':'c1 a'}
conf['params']['ffhadron']['cb1 b'] ={'value':   2.50452e+00, 'min':   0   , 'max':    30,'fixed':'c1 b'}
conf['params']['ffhadron']['cb1 c'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'c1 c'}
conf['params']['ffhadron']['cb1 d'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'c1 d'}
                                                                                                  
conf['params']['ffhadron']['bb1 N'] ={'value':   5.66971e-01, 'min':   0.0 , 'max':     1,'fixed':'b1 N'}
conf['params']['ffhadron']['bb1 a'] ={'value':  -1.80000e+00, 'min':  -1.8 , 'max':    10,'fixed':'b1 a'}
conf['params']['ffhadron']['bb1 b'] ={'value':   3.41727e+00, 'min':   0   , 'max':    30,'fixed':'b1 b'}
conf['params']['ffhadron']['bb1 c'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'b1 c'}
conf['params']['ffhadron']['bb1 d'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'b1 d'}



conf['params']['ffhadron']['g2 N']  ={'value':   2.33320e-01, 'min':   0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['g2 a']  ={'value':   1.48737e+00, 'min':  -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['g2 b']  ={'value':   9.62755e+00, 'min':   0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['g2 c']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['g2 d']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}

conf['params']['ffhadron']['u2 N']  ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['u2 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['u2 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['u2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['u2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
                                                              
conf['params']['ffhadron']['d2 N']  ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['d2 a']  ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['d2 b']  ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['d2 c']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['d2 d']  ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':False}

conf['params']['ffhadron']['s2 N']  ={'value':   2.51671e-01, 'min':   0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['s2 a']  ={'value':  -1.43444e+00, 'min':  -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['s2 b']  ={'value':   5.65143e+00, 'min':   0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['s2 c']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['s2 d']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['c2 N']  ={'value':   7.76923e-01, 'min':   0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['c2 a']  ={'value':  -1.80000e+00, 'min':  -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['c2 b']  ={'value':   2.50452e+00, 'min':   0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['c2 c']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['c2 d']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
                                                                                                   
conf['params']['ffhadron']['b2 N']  ={'value':   5.66971e-01, 'min':   0.0 , 'max':     1,'fixed':False}
conf['params']['ffhadron']['b2 a']  ={'value':  -1.80000e+00, 'min':  -1.8 , 'max':    10,'fixed':False}
conf['params']['ffhadron']['b2 b']  ={'value':   3.41727e+00, 'min':   0   , 'max':    30,'fixed':False}
conf['params']['ffhadron']['b2 c']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
conf['params']['ffhadron']['b2 d']  ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':False}
                                                                                                 
conf['params']['ffhadron']['ub2 N'] ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':'s2 N'}
conf['params']['ffhadron']['ub2 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'s2 a'}
conf['params']['ffhadron']['ub2 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'s2 b'}
conf['params']['ffhadron']['ub2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s2 c'}
conf['params']['ffhadron']['ub2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s2 d'}
                                                                     
conf['params']['ffhadron']['db2 N'] ={'value':    0.00000e+00, 'min':  0.0 , 'max':     1,'fixed':'s2 N'}
conf['params']['ffhadron']['db2 a'] ={'value':    1.76877e-01, 'min': -1.8 , 'max':    10,'fixed':'s2 a'}
conf['params']['ffhadron']['db2 b'] ={'value':    4.81521e+00, 'min':  0   , 'max':    30,'fixed':'s2 b'}
conf['params']['ffhadron']['db2 c'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s2 c'}
conf['params']['ffhadron']['db2 d'] ={'value':    4.81521e+00, 'min': -10  , 'max':    10,'fixed':'s2 d'}

conf['params']['ffhadron']['sb2 N'] ={'value':   2.51671e-01, 'min':  0.0  , 'max':     1,'fixed':'s2 N'}
conf['params']['ffhadron']['sb2 a'] ={'value':  -1.43444e+00, 'min':  -1.8 , 'max':    10,'fixed':'s2 a'}
conf['params']['ffhadron']['sb2 b'] ={'value':   1.65143e+00, 'min':   0   , 'max':    30,'fixed':'s2 b'}
conf['params']['ffhadron']['sb2 c'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'s2 c'}
conf['params']['ffhadron']['sb2 d'] ={'value':   4.81521e+00, 'min':  -10  , 'max':    10,'fixed':'s2 d'}
                                                                                                  
conf['params']['ffhadron']['cb2 N'] ={'value':   7.76923e-01, 'min':  0.0 , 'max':     1,'fixed':'c2 N'}
conf['params']['ffhadron']['cb2 a'] ={'value':  -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':'c2 a'}
conf['params']['ffhadron']['cb2 b'] ={'value':   2.50452e+00, 'min':  0   , 'max':    30,'fixed':'c2 b'}
conf['params']['ffhadron']['cb2 c'] ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c2 c'}
conf['params']['ffhadron']['cb2 d'] ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':'c2 d'}
                                                                                                  
conf['params']['ffhadron']['bb2 N'] ={'value':   5.66971e-01, 'min':  0.0 , 'max':     1,'fixed':'b2 N'}
conf['params']['ffhadron']['bb2 a'] ={'value':  -1.80000e+00, 'min': -1.8 , 'max':    10,'fixed':'b2 a'}
conf['params']['ffhadron']['bb2 b'] ={'value':   3.41727e+00, 'min':  0   , 'max':    30,'fixed':'b2 b'}
conf['params']['ffhadron']['bb2 c'] ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b2 c'}
conf['params']['ffhadron']['bb2 d'] ={'value':   4.81521e+00, 'min': -10  , 'max':    10,'fixed':'b2 d'}


#--steps
conf['steps']={}

istep=43
#--iterate on simultaneous fit (PDIS W2 > 10, Q2 > mc2, no Ape (w/ smearing, valence parameterization, OS+ for all, g2 = 0))
conf['ftol']=1e-3
conf['tmc']     = 'AOT'
conf['ht']      = True
conf['ht type'] = 'mult'
conf['offpdf']  = True

#conf['pidis nuc']      = True
#conf['pidis tmc']      = False
#conf['pidis gXres']    = False

#conf['lqcd_itd'] = True

#conf['fix_gA']   = True

conf['steps'][istep]={}
conf['steps'][istep]['dep']=[42]
conf['steps'][istep]['active distributions']=['ffpion','ffkaon','ffhadron']
conf['steps'][istep]['passive distributions']=['pdf','ht4','off pdf','ppdf']
#------------------------------------------------------------------------------------------------------------------
conf['steps'][istep]['datasets']={}
#conf['steps'][istep]['datasets']['pidis']=[]
#conf['steps'][istep]['datasets']['pidis'].append(10002) # 10002 | proton   | A1   | COMPASS         |          |
#conf['steps'][istep]['datasets']['pidis'].append(10003) # 10003 | proton   | A1   | COMPASS         |          |
#conf['steps'][istep]['datasets']['pidis'].append(10004) # 10004 | proton   | A1   | EMC             |          |
#conf['steps'][istep]['datasets']['pidis'].append(10007) # 10007 | proton   | Apa  | HERMES          |          |
#conf['steps'][istep]['datasets']['pidis'].append(10022) # 10022 | proton   | Apa  | SLAC(E143)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10029) # 10029 | proton   | Apa  | SLAC(E155)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10032) # 10032 | proton   | Apa  | SLACE80E130     |          |
#conf['steps'][istep]['datasets']['pidis'].append(10035) # 10035 | proton   | A1   | SMC             |          |
#conf['steps'][istep]['datasets']['pidis'].append(10036) # 10036 | proton   | A1   | SMC             |          |
#conf['steps'][istep]['datasets']['pidis'].append(10005) # 10005 | neutron  | A1   | HERMES          |          |
#conf['steps'][istep]['datasets']['pidis'].append(10001) # 10001 | deuteron | A1   | COMPASS         |          |
#conf['steps'][istep]['datasets']['pidis'].append(10006) # 10006 | deuteron | Apa  | HERMES          |          |
#conf['steps'][istep]['datasets']['pidis'].append(10021) # 10021 | deuteron | Apa  | SLAC(E143)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10027) # 10027 | deuteron | Apa  | SLAC(E155)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10033) # 10033 | deuteron | A1   | SMC             |          |
#conf['steps'][istep]['datasets']['pidis'].append(10034) # 10034 | deuteron | A1   | SMC             |          |
#conf['steps'][istep]['datasets']['pidis'].append(10016) # 10016 | deuteron | Apa  | JLabHB(EG1DVCS) |          |
#conf['steps'][istep]['datasets']['pidis'].append(10017) # 10017 | proton   | Apa  | JLabHB(EG1DVCS) |          |
#conf['steps'][istep]['datasets']['pidis'].append(10039) # 10039 | deuteron | Apa  | JLabHB(EG1b)    | E =4 GeV |
#conf['steps'][istep]['datasets']['pidis'].append(10040) # 10040 | deuteron | Apa  | JLabHB(EG1b)    | E =5 GeV |
#conf['steps'][istep]['datasets']['pidis'].append(10043) # 10043 | proton   | Apa  | JLabHB(EG1b)    | E =4 GeV |
#conf['steps'][istep]['datasets']['pidis'].append(10044) # 10044 | proton   | Apa  | JLabHB(EG1b)    | E =5 GeV |
#conf['steps'][istep]['datasets']['pidis'].append(10010) # 10010 | helium   | Apa  | JLabHA(E06-014) |          |
#conf['steps'][istep]['datasets']['pidis'].append(10014) # 10014 | helium   | Apa  | JLabHA(E99-117) |          |
#conf['steps'][istep]['datasets']['pidis'].append(10018) # 10018 | helium   | A1   | SLAC(E142)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10025) # 10025 | helium   | Apa  | SLAC(E154)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10008) # 10008 | proton   | A2   | HERMES          |          |
#conf['steps'][istep]['datasets']['pidis'].append(10023) # 10023 | proton   | Ape  | SLAC(E143)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10028) # 10028 | proton   | Ape  | SLAC(E155)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10031) # 10031 | proton   | Atpe | SLAC(E155x)     |          |
#conf['steps'][istep]['datasets']['pidis'].append(10020) # 10020 | deuteron | Ape  | SLAC(E143)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10026) # 10026 | deuteron | Ape  | SLAC(E155)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10030) # 10030 | deuteron | Atpe | SLAC(E155x)     |          |
#conf['steps'][istep]['datasets']['pidis'].append(10011) # 10011 | helium   | Ape  | JLabHA(E06-014) |          |
#conf['steps'][istep]['datasets']['pidis'].append(10015) # 10015 | helium   | Ape  | JLabHA(E99-117) |          |
#conf['steps'][istep]['datasets']['pidis'].append(10019) # 10019 | helium   | A2   | SLAC(E142)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(10024) # 10024 | helium   | Ape  | SLAC(E154)      |          |
#conf['steps'][istep]['datasets']['pidis'].append(20001) # SANE      |          |
#conf['steps'][istep]['datasets']['pidis'].append(20002) # SANE      |          |
#conf['steps'][istep]['datasets']['pidis'].append(20003) # SANE      |          |
#conf['steps'][istep]['datasets']['pidis'].append(20004) # SANE      |          |
#conf['steps'][istep]['datasets']['pjet'] = []
#conf['steps'][istep]['datasets']['pjet'].append(20001) ## STAR 2006 paper on 2003 and 2004 data
#conf['steps'][istep]['datasets']['pjet'].append(20002) ## STAR 2012 paper on 2005 data
#conf['steps'][istep]['datasets']['pjet'].append(20003) ## STAR 2012 paper on 2006 data
#conf['steps'][istep]['datasets']['pjet'].append(20004) ## STAR 2015 paper on 2009 data
#conf['steps'][istep]['datasets']['pjet'].append(20005) ## PHENIX 2011 paper on 2005 data
#conf['steps'][istep]['datasets']['pjet'].append(20006) ## STAR 2019 paper on 2012 data
#conf['steps'][istep]['datasets']['pjet'].append(20007) ## STAR 2021 paper on 2015 data
#conf['steps'][istep]['datasets']['pjet'].append(20008) ## STAR 2021 paper on 2013 data
conf['steps'][istep]['datasets']['wzrv'] = []
#conf['steps'][istep]['datasets']['wzrv'].append(1000) 
#conf['steps'][istep]['datasets']['wzrv'].append(1020) 
#conf['steps'][istep]['datasets']['wzrv'].append(1021) 

conf['steps'][istep]['datasets']['sia']=[]

conf['steps'][istep]['datasets']['sia'].append(1001)  # TASSO (12 GeV)
conf['steps'][istep]['datasets']['sia'].append(1002)  # TASSO (14 Gev)
conf['steps'][istep]['datasets']['sia'].append(1003)  # TASSO (22 GeV)
conf['steps'][istep]['datasets']['sia'].append(1004)  # TASSO (30 GeV)
conf['steps'][istep]['datasets']['sia'].append(1005)  # TASSO (34 GeV)
conf['steps'][istep]['datasets']['sia'].append(1006)  # TASSO (44 GeV)
###conf['steps'][istep]['datasets']['sia'].append(1007)  # TPC 1984
###conf['steps'][istep]['datasets']['sia'].append(1008)  # TPC 1988
conf['steps'][istep]['datasets']['sia'].append(1010)  # TPC(c)
conf['steps'][istep]['datasets']['sia'].append(1011)  # TPC(b)
conf['steps'][istep]['datasets']['sia'].append(1013)  # TOPAZ
conf['steps'][istep]['datasets']['sia'].append(1014)  # SLD
conf['steps'][istep]['datasets']['sia'].append(1015)  # SLD(uds)
conf['steps'][istep]['datasets']['sia'].append(1016)  # SLD(c)
conf['steps'][istep]['datasets']['sia'].append(1017)  # SLD(b)
conf['steps'][istep]['datasets']['sia'].append(1018)  # ALEPH
conf['steps'][istep]['datasets']['sia'].append(1019)  # OPAL
###conf['steps'][istep]['datasets']['sia'].append(1023)  # OPAL(c)
###conf['steps'][istep]['datasets']['sia'].append(1024)  # OPAL(b)
conf['steps'][istep]['datasets']['sia'].append(1025)  # DELPHI
conf['steps'][istep]['datasets']['sia'].append(1026)  # DELPHI(uds)
conf['steps'][istep]['datasets']['sia'].append(1027)  # DELPHI(b)
conf['steps'][istep]['datasets']['sia'].append(1028)  # BABAR
conf['steps'][istep]['datasets']['sia'].append(1029)  # BELL
conf['steps'][istep]['datasets']['sia'].append(1030)  # ARGUS

conf['steps'][istep]['datasets']['sia'].append(2001)  # TASSO (12 GeV)
conf['steps'][istep]['datasets']['sia'].append(2002)  # TASSO (14 GeV)
conf['steps'][istep]['datasets']['sia'].append(2003)  # TASSO (22 GeV)
conf['steps'][istep]['datasets']['sia'].append(2004)  # TASSO (30 GeV)
conf['steps'][istep]['datasets']['sia'].append(2005)  # TASSO (34 GeV)
conf['steps'][istep]['datasets']['sia'].append(2006)  # TASSO (44 GeV)
###conf['steps'][istep]['datasets']['sia'].append(2007)  # TPC 1984
conf['steps'][istep]['datasets']['sia'].append(2008)  # TPC 1988
conf['steps'][istep]['datasets']['sia'].append(2013)  # TOPAZ
conf['steps'][istep]['datasets']['sia'].append(2014)  # SLD
conf['steps'][istep]['datasets']['sia'].append(2015)  # SLD(uds)
conf['steps'][istep]['datasets']['sia'].append(2016)  # SLD(c)
conf['steps'][istep]['datasets']['sia'].append(2017)  # SLD(b)
conf['steps'][istep]['datasets']['sia'].append(2018)  # ALEPH
conf['steps'][istep]['datasets']['sia'].append(2019)  # OPAL
###conf['steps'][istep]['datasets']['sia'].append(2023)  # OPAL(c)
###conf['steps'][istep]['datasets']['sia'].append(2024)  # OPAL(b)
conf['steps'][istep]['datasets']['sia'].append(2025)  # DELPHI
conf['steps'][istep]['datasets']['sia'].append(2026)  # DELPHI(uds)
conf['steps'][istep]['datasets']['sia'].append(2027)  # DELPHI(b)
conf['steps'][istep]['datasets']['sia'].append(2028)  # BABAR
conf['steps'][istep]['datasets']['sia'].append(2029)  # BELL
conf['steps'][istep]['datasets']['sia'].append(2030)  # ARGUS

conf['steps'][istep]['datasets']['sia'].append(4000)  # ALEPH
conf['steps'][istep]['datasets']['sia'].append(4001)  # DELPHI
conf['steps'][istep]['datasets']['sia'].append(4002)  # SLD
conf['steps'][istep]['datasets']['sia'].append(4003)  # TASSO (12 GeV)
conf['steps'][istep]['datasets']['sia'].append(4004)  # TPC 1988
###conf['steps'][istep]['datasets']['sia'].append(4005)  # OPAL(b)
###conf['steps'][istep]['datasets']['sia'].append(4006)  # OPAL(c)
conf['steps'][istep]['datasets']['sia'].append(4007)  # OPAL
conf['steps'][istep]['datasets']['sia'].append(4008)  # TASSO (14 GeV)
conf['steps'][istep]['datasets']['sia'].append(4009)  # TASSO (22 GeV)
conf['steps'][istep]['datasets']['sia'].append(4010)  # TASSO (30 GeV)
###conf['steps'][istep]['datasets']['sia'].append(4011)  # TASSO (34 GeV)
conf['steps'][istep]['datasets']['sia'].append(4012)  # TASSO (44 GeV)
conf['steps'][istep]['datasets']['sia'].append(4013)  # DELPHI(b)
conf['steps'][istep]['datasets']['sia'].append(4014)  # SLD(c)
conf['steps'][istep]['datasets']['sia'].append(4015)  # SLD(b)
conf['steps'][istep]['datasets']['sia'].append(4016)  # SLD(uds)
conf['steps'][istep]['datasets']['sia'].append(4017)  # DELPHI(uds)


#conf['steps'][istep]['datasets']['psidis']=[]
#conf['steps'][istep]['datasets']['psidis'].append(20004) # 20004 | proton   | A1pi+  | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20005) # 20005 | proton   | A1pi-  | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20008) # 20008 | deuteron | A1pi+  | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20009) # 20009 | deuteron | A1pi-  | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20017) # 20017 | proton   | A1pi+  | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20018) # 20018 | proton   | A1pi-  | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20021) # 20021 | deuteron | A1pi+  | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20022) # 20022 | deuteron | A1pi-  | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20012) # 20012 | deuteron | A1K+   | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20013) # 20013 | deuteron | A1K-   | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20014) # 20014 | deuteron | A1Ksum | HERMES  
#conf['steps'][istep]['datasets']['psidis'].append(20019) # 20019 | proton   | A1K+   | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20020) # 20020 | proton   | A1K-   | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20025) # 20025 | deuteron | A1K+   | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20026) # 20026 | deuteron | A1K-   | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20000) # 20000 | proton   | A1h+   | SMC 
#conf['steps'][istep]['datasets']['psidis'].append(20001) # 20001 | proton   | A1h-   | SMC 
#conf['steps'][istep]['datasets']['psidis'].append(20002) # 20002 | deuteron | A1h+   | SMC 
#conf['steps'][istep]['datasets']['psidis'].append(20003) # 20003 | deuteron | A1h-   | SMC 
#conf['steps'][istep]['datasets']['psidis'].append(20006) # 20006 | proton   | A1h+   | HERMES 
#conf['steps'][istep]['datasets']['psidis'].append(20007) # 20007 | proton   | A1h-   | HERMES 
#conf['steps'][istep]['datasets']['psidis'].append(20010) # 20010 | deuteron | A1h+   | HERMES 
#conf['steps'][istep]['datasets']['psidis'].append(20011) # 20011 | deuteron | A1h-   | HERMES 
#conf['steps'][istep]['datasets']['psidis'].append(20015) # 20015 | helium   | A1h+   | HERMES 
#conf['steps'][istep]['datasets']['psidis'].append(20016) # 20016 | helium   | A1h-   | HERMES 
#conf['steps'][istep]['datasets']['psidis'].append(20023) # 20023 | deuteron | A1h+   | COMPASS 
#conf['steps'][istep]['datasets']['psidis'].append(20024) # 20024 | deuteron | A1h-   | COMPASS 
#conf['steps'][istep]['datasets']['sidis']=[]
#conf['steps'][istep]['datasets']['sidis'].append(1005) # deuteron , mult , pi+ , COMPASS
#conf['steps'][istep]['datasets']['sidis'].append(1006) # deuteron , mult , pi- , COMPASS
#conf['steps'][istep]['datasets']['sidis'].append(2005) # deuteron , mult , K+  , COMPASS
#conf['steps'][istep]['datasets']['sidis'].append(2006) # deuteron , mult , K-  , COMPASS
#conf['steps'][istep]['datasets']['sidis'].append(3000) # deuteron , mult , h+  , COMPASS
#conf['steps'][istep]['datasets']['sidis'].append(3001) # deuteron , mult , h-  , COMPASS
#conf['steps'][istep]['datasets']['idis']=[]
#conf['steps'][istep]['datasets']['idis'].append(10010) # proton   | F2            | SLAC
#conf['steps'][istep]['datasets']['idis'].append(10011) # deuteron | F2            | SLAC
#conf['steps'][istep]['datasets']['idis'].append(10016) # proton   | F2            | BCDMS
#conf['steps'][istep]['datasets']['idis'].append(10017) # deuteron | F2            | BCDMS
#conf['steps'][istep]['datasets']['idis'].append(10020) # proton   | F2            | NMC
#conf['steps'][istep]['datasets']['idis'].append(10021) # d/p      | F2d/F2p       | NMC
#conf['steps'][istep]['datasets']['idis'].append(10026) # proton   | sigma red     | HERA II NC e+ (1)
#conf['steps'][istep]['datasets']['idis'].append(10027) # proton   | sigma red     | HERA II NC e+ (2)
#conf['steps'][istep]['datasets']['idis'].append(10028) # proton   | sigma red     | HERA II NC e+ (3)
#conf['steps'][istep]['datasets']['idis'].append(10029) # proton   | sigma red     | HERA II NC e+ (4)
#conf['steps'][istep]['datasets']['idis'].append(10030) # proton   | sigma red     | HERA II NC e-
#conf['steps'][istep]['datasets']['idis'].append(10031) # proton   | sigma red     | HERA II CC e+
#conf['steps'][istep]['datasets']['idis'].append(10032) # proton   | sigma red     | HERA II CC e-
#conf['steps'][istep]['datasets']['idis'].append(10002) # deuteron | F2            | JLab Hall C (E00-106)
#conf['steps'][istep]['datasets']['idis'].append(10003) # proton   | F2            | JLab Hall C (E00-106)
#conf['steps'][istep]['datasets']['idis'].append(10033) # n/d      | F2n/F2d       | BONUS
#conf['steps'][istep]['datasets']['idis'].append(10041) # h/d      | F2h/F2d       | JLab Hall C (E03-103)
#conf['steps'][istep]['datasets']['idis'].append(10050) # d/p      | F2d/F2p       | MARATHON
#conf['steps'][istep]['datasets']['idis'].append(10051) # h/t      | F2h/F2t       | MARATHON
#conf['steps'][istep]['datasets']['dy']=[]
#conf['steps'][istep]['datasets']['dy'].append(10001)
#conf['steps'][istep]['datasets']['dy'].append(20001)
#conf['steps'][istep]['datasets']['dy'].append(20002)
#conf['steps'][istep]['datasets']['zrap']=[]
#conf['steps'][istep]['datasets']['zrap'].append(1000)
#conf['steps'][istep]['datasets']['zrap'].append(1001)
#conf['steps'][istep]['datasets']['wasym']=[]
#conf['steps'][istep]['datasets']['wasym'].append(1000)
#conf['steps'][istep]['datasets']['wasym'].append(1001)
#conf['steps'][istep]['datasets']['wzrv'].append(2010)
#conf['steps'][istep]['datasets']['wzrv'].append(2011)
#conf['steps'][istep]['datasets']['wzrv'].append(2012)
#conf['steps'][istep]['datasets']['wzrv'].append(2013)
#conf['steps'][istep]['datasets']['wzrv'].append(2014)
#conf['steps'][istep]['datasets']['wzrv'].append(2016)
#conf['steps'][istep]['datasets']['wzrv'].append(2017)
#conf['steps'][istep]['datasets']['wzrv'].append(2020)
#conf['steps'][istep]['datasets']['jet'] = []
#conf['steps'][istep]['datasets']['jet'].append(10001) ## D0 dataset
#conf['steps'][istep]['datasets']['jet'].append(10002) ## CDF dataset
#conf['steps'][istep]['datasets']['jet'].append(10003) ## STAR MB dataset
#conf['steps'][istep]['datasets']['jet'].append(10004) ## STAR HT dataset
#
#conf['steps'][istep]['datasets']['wc'] = []
#conf['steps'][istep]['datasets']['wc'].append(10001)  # pp | sig | ATLAS7Wm
#conf['steps'][istep]['datasets']['wc'].append(10002)  # pp | sig | ATLAS7Wp
#conf['steps'][istep]['datasets']['wc'].append(20001)  # pp | Wp+Wm | CMS7WpWm_add
#conf['steps'][istep]['datasets']['wc'].append(20002)  # pp | Wp/Wm | CMS7WpWm_ratio
#conf['steps'][istep]['datasets']['wc'].append(20003)  # pp | Wp+Wm | CMS13WpWm_add

conf['FILT'] = {_:[] for _ in ['exp','par','value']}

conf['FILT']['exp'].append(('wc',10001,10))
conf['FILT']['exp'].append(('wc',10002,10))
conf['FILT']['exp'].append(('wc',20001,10))
conf['FILT']['exp'].append(('wc',20002,10))
conf['FILT']['exp'].append(('wc',20003,10))


conf['FILT']['exp'].append(('sidis',1005,10))
conf['FILT']['exp'].append(('sidis',1006,10))
conf['FILT']['exp'].append(('sidis',2005,10))
conf['FILT']['exp'].append(('sidis',2006,10))
conf['FILT']['exp'].append(('sidis',3000,10))
conf['FILT']['exp'].append(('sidis',3001,10))
conf['FILT']['exp'].append(('idis' ,10010,10)) 
conf['FILT']['exp'].append(('idis' ,10011,10)) 
conf['FILT']['exp'].append(('idis' ,10016,10)) 
conf['FILT']['exp'].append(('idis' ,10017,10)) 
conf['FILT']['exp'].append(('idis' ,10020,10)) 
conf['FILT']['exp'].append(('idis' ,10021,10)) 
conf['FILT']['exp'].append(('idis' ,10026,10)) 
conf['FILT']['exp'].append(('idis' ,10027,10)) 
conf['FILT']['exp'].append(('idis' ,10028,10)) 
conf['FILT']['exp'].append(('idis' ,10029,10)) 
conf['FILT']['exp'].append(('idis' ,10030,10)) 
conf['FILT']['exp'].append(('idis' ,10031,10)) 
conf['FILT']['exp'].append(('idis' ,10032,10)) 
conf['FILT']['exp'].append(('idis' ,10002,10)) 
conf['FILT']['exp'].append(('idis' ,10003,10)) 
conf['FILT']['exp'].append(('idis' ,10033,10)) 
conf['FILT']['exp'].append(('idis' ,10041,10)) 
conf['FILT']['exp'].append(('idis' ,10050,10)) 
conf['FILT']['exp'].append(('idis' ,10051,10)) 
conf['FILT']['exp'].append(('dy'   ,10001,10))
conf['FILT']['exp'].append(('dy'   ,20001,10))
conf['FILT']['exp'].append(('dy'   ,20002,10))
conf['FILT']['exp'].append(('zrap' ,1000 ,10))
conf['FILT']['exp'].append(('zrap' ,1001 ,10))
conf['FILT']['exp'].append(('wasym',1000 ,10))
conf['FILT']['exp'].append(('wasym',1001 ,10))
conf['FILT']['exp'].append(('wzrv' ,2010 ,10))
conf['FILT']['exp'].append(('wzrv' ,2011 ,10))
conf['FILT']['exp'].append(('wzrv' ,2012 ,10))
conf['FILT']['exp'].append(('wzrv' ,2013 ,10))
conf['FILT']['exp'].append(('wzrv' ,2014 ,10))
conf['FILT']['exp'].append(('wzrv' ,2016 ,10))
conf['FILT']['exp'].append(('wzrv' ,2017 ,10))
conf['FILT']['exp'].append(('wzrv' ,2020 ,10))
conf['FILT']['exp'].append(('jet'  ,10001,10)) 
conf['FILT']['exp'].append(('jet'  ,10002,10)) 
conf['FILT']['exp'].append(('jet'  ,10003,10)) 
conf['FILT']['exp'].append(('jet'  ,10004,10)) 
conf['FILT']['exp'].append(('sia'  ,1001 ,20))
conf['FILT']['exp'].append(('sia'  ,1002 ,20))
conf['FILT']['exp'].append(('sia'  ,1003 ,20))
conf['FILT']['exp'].append(('sia'  ,1004 ,20))
conf['FILT']['exp'].append(('sia'  ,1005 ,20))
conf['FILT']['exp'].append(('sia'  ,1006 ,20))
conf['FILT']['exp'].append(('sia'  ,1007 ,20))
conf['FILT']['exp'].append(('sia'  ,1008 ,20))
conf['FILT']['exp'].append(('sia'  ,1012 ,20))
conf['FILT']['exp'].append(('sia'  ,1013 ,20))
conf['FILT']['exp'].append(('sia'  ,1014 ,20))
conf['FILT']['exp'].append(('sia'  ,1018 ,20))
conf['FILT']['exp'].append(('sia'  ,1019 ,20))
conf['FILT']['exp'].append(('sia'  ,1025 ,20))
conf['FILT']['exp'].append(('sia'  ,1028 ,20))
conf['FILT']['exp'].append(('sia'  ,1029 ,20))
conf['FILT']['exp'].append(('sia'  ,1030 ,20))
conf['FILT']['exp'].append(('sia'  ,1010 ,20))
conf['FILT']['exp'].append(('sia'  ,1011 ,20))
conf['FILT']['exp'].append(('sia'  ,1016 ,20))
conf['FILT']['exp'].append(('sia'  ,1017 ,20))
conf['FILT']['exp'].append(('sia'  ,1023 ,20))
conf['FILT']['exp'].append(('sia'  ,1024 ,20))
conf['FILT']['exp'].append(('sia'  ,1027 ,20))
conf['FILT']['exp'].append(('sia'  ,2001 ,20))
conf['FILT']['exp'].append(('sia'  ,2002 ,20))
conf['FILT']['exp'].append(('sia'  ,2003 ,20))
conf['FILT']['exp'].append(('sia'  ,2004 ,20))
conf['FILT']['exp'].append(('sia'  ,2005 ,20))
conf['FILT']['exp'].append(('sia'  ,2006 ,20))
conf['FILT']['exp'].append(('sia'  ,2007 ,20))
conf['FILT']['exp'].append(('sia'  ,2008 ,20))
conf['FILT']['exp'].append(('sia'  ,2012 ,20))
conf['FILT']['exp'].append(('sia'  ,2013 ,20))
conf['FILT']['exp'].append(('sia'  ,2014 ,20))
conf['FILT']['exp'].append(('sia'  ,2016 ,20))
conf['FILT']['exp'].append(('sia'  ,2017 ,20))
conf['FILT']['exp'].append(('sia'  ,2018 ,20))
conf['FILT']['exp'].append(('sia'  ,2019 ,20))
conf['FILT']['exp'].append(('sia'  ,2023 ,20))
conf['FILT']['exp'].append(('sia'  ,2024 ,20))
conf['FILT']['exp'].append(('sia'  ,2025 ,20))
conf['FILT']['exp'].append(('sia'  ,2027 ,20))
conf['FILT']['exp'].append(('sia'  ,2028 ,20))
conf['FILT']['exp'].append(('sia'  ,2029 ,20))
conf['FILT']['exp'].append(('sia'  ,2030 ,20))
conf['FILT']['exp'].append(('sia'  ,2031 ,20))
conf['FILT']['exp'].append(('sia ' ,4000 ,10))
conf['FILT']['exp'].append(('sia'  ,4001 ,10))
conf['FILT']['exp'].append(('sia'  ,4002 ,10))
conf['FILT']['exp'].append(('sia'  ,4003 ,10))
conf['FILT']['exp'].append(('sia'  ,4004 ,10))
conf['FILT']['exp'].append(('sia'  ,4005 ,10))
conf['FILT']['exp'].append(('sia'  ,4006 ,10))
conf['FILT']['exp'].append(('sia'  ,4007 ,10))
conf['FILT']['exp'].append(('sia'  ,4008 ,10))
conf['FILT']['exp'].append(('sia'  ,4009 ,10))
conf['FILT']['exp'].append(('sia'  ,4010 ,10))
conf['FILT']['exp'].append(('sia'  ,4011 ,10))
conf['FILT']['exp'].append(('sia'  ,4012 ,10))
conf['FILT']['exp'].append(('sia'  ,4013 ,10))
conf['FILT']['exp'].append(('sia'  ,4014 ,10))
conf['FILT']['exp'].append(('sia'  ,4015 ,10))
conf['FILT']['exp'].append(('sia'  ,4016 ,10))
conf['FILT']['exp'].append(('sia'  ,4017 ,10))

conf['FILT']['exp'].append(('pidis',10002,10))
conf['FILT']['exp'].append(('pidis',10003,10))
conf['FILT']['exp'].append(('pidis',10004,10))
conf['FILT']['exp'].append(('pidis',10007,10))
conf['FILT']['exp'].append(('pidis',10022,10))
conf['FILT']['exp'].append(('pidis',10029,10))
conf['FILT']['exp'].append(('pidis',10032,10))
conf['FILT']['exp'].append(('pidis',10035,10))
conf['FILT']['exp'].append(('pidis',10036,10))
conf['FILT']['exp'].append(('pidis',10005,10))
conf['FILT']['exp'].append(('pidis',10001,10))
conf['FILT']['exp'].append(('pidis',10006,10))
conf['FILT']['exp'].append(('pidis',10021,10))
conf['FILT']['exp'].append(('pidis',10027,10))
conf['FILT']['exp'].append(('pidis',10033,10))
conf['FILT']['exp'].append(('pidis',10034,10))
conf['FILT']['exp'].append(('pidis',10016,10))
conf['FILT']['exp'].append(('pidis',10017,10))
conf['FILT']['exp'].append(('pidis',10039,10))
conf['FILT']['exp'].append(('pidis',10040,10))
conf['FILT']['exp'].append(('pidis',10043,10))
conf['FILT']['exp'].append(('pidis',10044,10))
conf['FILT']['exp'].append(('pidis',10010,10))
conf['FILT']['exp'].append(('pidis',10014,10))
conf['FILT']['exp'].append(('pidis',10018,10))
conf['FILT']['exp'].append(('pidis',10025,10))
conf['FILT']['exp'].append(('pidis',10008,10))
conf['FILT']['exp'].append(('pidis',10023,10))
conf['FILT']['exp'].append(('pidis',10028,10))
conf['FILT']['exp'].append(('pidis',10031,10))
conf['FILT']['exp'].append(('pidis',10020,10))
conf['FILT']['exp'].append(('pidis',10026,10))
conf['FILT']['exp'].append(('pidis',10030,10))
conf['FILT']['exp'].append(('pidis',10011,10))
conf['FILT']['exp'].append(('pidis',10015,10))
conf['FILT']['exp'].append(('pidis',10019,10))
conf['FILT']['exp'].append(('pidis',10024,10))
conf['FILT']['exp'].append(('pidis',20001,10))
conf['FILT']['exp'].append(('pidis',20002,10))
conf['FILT']['exp'].append(('pidis',20003,10))
conf['FILT']['exp'].append(('pidis',20004,10))

conf['FILT']['exp'].append(('psidis',20012,10))
conf['FILT']['exp'].append(('psidis',20013,10))
conf['FILT']['exp'].append(('psidis',20014,10))
conf['FILT']['exp'].append(('psidis',20019,10))
conf['FILT']['exp'].append(('psidis',20020,10))
conf['FILT']['exp'].append(('psidis',20025,10))
conf['FILT']['exp'].append(('psidis',20026,10))
conf['FILT']['exp'].append(('psidis',20004,10))
conf['FILT']['exp'].append(('psidis',20005,10))
conf['FILT']['exp'].append(('psidis',20008,10))
conf['FILT']['exp'].append(('psidis',20009,10))
conf['FILT']['exp'].append(('psidis',20017,10))
conf['FILT']['exp'].append(('psidis',20018,10))
conf['FILT']['exp'].append(('psidis',20021,10))
conf['FILT']['exp'].append(('psidis',20022,10))
conf['FILT']['exp'].append(('psidis',20000,10))
conf['FILT']['exp'].append(('psidis',20001,10))
conf['FILT']['exp'].append(('psidis',20002,10))
conf['FILT']['exp'].append(('psidis',20003,10))
conf['FILT']['exp'].append(('psidis',20006,10))
conf['FILT']['exp'].append(('psidis',20007,10))
conf['FILT']['exp'].append(('psidis',20010,10))
conf['FILT']['exp'].append(('psidis',20011,10))
conf['FILT']['exp'].append(('psidis',20015,10))
conf['FILT']['exp'].append(('psidis',20016,10))
conf['FILT']['exp'].append(('psidis',20023,10))
conf['FILT']['exp'].append(('psidis',20024,10))

conf['FILT']['exp'].append(('pjet' ,20001,10)) 
conf['FILT']['exp'].append(('pjet' ,20002,10)) 
conf['FILT']['exp'].append(('pjet' ,20003,10)) 
conf['FILT']['exp'].append(('pjet' ,20004,10)) 
conf['FILT']['exp'].append(('pjet' ,20005,10)) 
conf['FILT']['exp'].append(('pjet' ,20006,10)) 
conf['FILT']['exp'].append(('pjet' ,20007,10)) 
conf['FILT']['exp'].append(('pjet' ,20008,10))

conf['FILT']['exp'].append(('wzrv' ,1000 ,10))
conf['FILT']['exp'].append(('wzrv' ,1020 ,10))
conf['FILT']['exp'].append(('wzrv' ,1021 ,10))










