#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import glob

import h5py
import numpy as np

import healpy
from scipy.stats import binned_statistic_2d

from time import time

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

from coordinate_grids import HealpixGrid, PolarGrid, CartesianGrid

parameters = {'xp': [22.674387857633945, 11.221508276482126, -14.356589454756971, -47.18817758739222, -34.59001004812521, -6.0516291196352654],
              'yp': [-20.102530541012214, -33.08246073298429, -42.24476439790574, -5.595549738219887, 42.5065445026178, 40.13395157135497],
              'min_JK': 18.8671875,
              'max_JK': 92.5}

cm_data = [[ 0.26700401,  0.00487433,  0.32941519],
       [ 0.26851048,  0.00960483,  0.33542652],
       [ 0.26994384,  0.01462494,  0.34137895],
       [ 0.27130489,  0.01994186,  0.34726862],
       [ 0.27259384,  0.02556309,  0.35309303],
       [ 0.27380934,  0.03149748,  0.35885256],
       [ 0.27495242,  0.03775181,  0.36454323],
       [ 0.27602238,  0.04416723,  0.37016418],
       [ 0.2770184 ,  0.05034437,  0.37571452],
       [ 0.27794143,  0.05632444,  0.38119074],
       [ 0.27879067,  0.06214536,  0.38659204],
       [ 0.2795655 ,  0.06783587,  0.39191723],
       [ 0.28026658,  0.07341724,  0.39716349],
       [ 0.28089358,  0.07890703,  0.40232944],
       [ 0.28144581,  0.0843197 ,  0.40741404],
       [ 0.28192358,  0.08966622,  0.41241521],
       [ 0.28232739,  0.09495545,  0.41733086],
       [ 0.28265633,  0.10019576,  0.42216032],
       [ 0.28291049,  0.10539345,  0.42690202],
       [ 0.28309095,  0.11055307,  0.43155375],
       [ 0.28319704,  0.11567966,  0.43611482],
       [ 0.28322882,  0.12077701,  0.44058404],
       [ 0.28318684,  0.12584799,  0.44496   ],
       [ 0.283072  ,  0.13089477,  0.44924127],
       [ 0.28288389,  0.13592005,  0.45342734],
       [ 0.28262297,  0.14092556,  0.45751726],
       [ 0.28229037,  0.14591233,  0.46150995],
       [ 0.28188676,  0.15088147,  0.46540474],
       [ 0.28141228,  0.15583425,  0.46920128],
       [ 0.28086773,  0.16077132,  0.47289909],
       [ 0.28025468,  0.16569272,  0.47649762],
       [ 0.27957399,  0.17059884,  0.47999675],
       [ 0.27882618,  0.1754902 ,  0.48339654],
       [ 0.27801236,  0.18036684,  0.48669702],
       [ 0.27713437,  0.18522836,  0.48989831],
       [ 0.27619376,  0.19007447,  0.49300074],
       [ 0.27519116,  0.1949054 ,  0.49600488],
       [ 0.27412802,  0.19972086,  0.49891131],
       [ 0.27300596,  0.20452049,  0.50172076],
       [ 0.27182812,  0.20930306,  0.50443413],
       [ 0.27059473,  0.21406899,  0.50705243],
       [ 0.26930756,  0.21881782,  0.50957678],
       [ 0.26796846,  0.22354911,  0.5120084 ],
       [ 0.26657984,  0.2282621 ,  0.5143487 ],
       [ 0.2651445 ,  0.23295593,  0.5165993 ],
       [ 0.2636632 ,  0.23763078,  0.51876163],
       [ 0.26213801,  0.24228619,  0.52083736],
       [ 0.26057103,  0.2469217 ,  0.52282822],
       [ 0.25896451,  0.25153685,  0.52473609],
       [ 0.25732244,  0.2561304 ,  0.52656332],
       [ 0.25564519,  0.26070284,  0.52831152],
       [ 0.25393498,  0.26525384,  0.52998273],
       [ 0.25219404,  0.26978306,  0.53157905],
       [ 0.25042462,  0.27429024,  0.53310261],
       [ 0.24862899,  0.27877509,  0.53455561],
       [ 0.2468114 ,  0.28323662,  0.53594093],
       [ 0.24497208,  0.28767547,  0.53726018],
       [ 0.24311324,  0.29209154,  0.53851561],
       [ 0.24123708,  0.29648471,  0.53970946],
       [ 0.23934575,  0.30085494,  0.54084398],
       [ 0.23744138,  0.30520222,  0.5419214 ],
       [ 0.23552606,  0.30952657,  0.54294396],
       [ 0.23360277,  0.31382773,  0.54391424],
       [ 0.2316735 ,  0.3181058 ,  0.54483444],
       [ 0.22973926,  0.32236127,  0.54570633],
       [ 0.22780192,  0.32659432,  0.546532  ],
       [ 0.2258633 ,  0.33080515,  0.54731353],
       [ 0.22392515,  0.334994  ,  0.54805291],
       [ 0.22198915,  0.33916114,  0.54875211],
       [ 0.22005691,  0.34330688,  0.54941304],
       [ 0.21812995,  0.34743154,  0.55003755],
       [ 0.21620971,  0.35153548,  0.55062743],
       [ 0.21429757,  0.35561907,  0.5511844 ],
       [ 0.21239477,  0.35968273,  0.55171011],
       [ 0.2105031 ,  0.36372671,  0.55220646],
       [ 0.20862342,  0.36775151,  0.55267486],
       [ 0.20675628,  0.37175775,  0.55311653],
       [ 0.20490257,  0.37574589,  0.55353282],
       [ 0.20306309,  0.37971644,  0.55392505],
       [ 0.20123854,  0.38366989,  0.55429441],
       [ 0.1994295 ,  0.38760678,  0.55464205],
       [ 0.1976365 ,  0.39152762,  0.55496905],
       [ 0.19585993,  0.39543297,  0.55527637],
       [ 0.19410009,  0.39932336,  0.55556494],
       [ 0.19235719,  0.40319934,  0.55583559],
       [ 0.19063135,  0.40706148,  0.55608907],
       [ 0.18892259,  0.41091033,  0.55632606],
       [ 0.18723083,  0.41474645,  0.55654717],
       [ 0.18555593,  0.4185704 ,  0.55675292],
       [ 0.18389763,  0.42238275,  0.55694377],
       [ 0.18225561,  0.42618405,  0.5571201 ],
       [ 0.18062949,  0.42997486,  0.55728221],
       [ 0.17901879,  0.43375572,  0.55743035],
       [ 0.17742298,  0.4375272 ,  0.55756466],
       [ 0.17584148,  0.44128981,  0.55768526],
       [ 0.17427363,  0.4450441 ,  0.55779216],
       [ 0.17271876,  0.4487906 ,  0.55788532],
       [ 0.17117615,  0.4525298 ,  0.55796464],
       [ 0.16964573,  0.45626209,  0.55803034],
       [ 0.16812641,  0.45998802,  0.55808199],
       [ 0.1666171 ,  0.46370813,  0.55811913],
       [ 0.16511703,  0.4674229 ,  0.55814141],
       [ 0.16362543,  0.47113278,  0.55814842],
       [ 0.16214155,  0.47483821,  0.55813967],
       [ 0.16066467,  0.47853961,  0.55811466],
       [ 0.15919413,  0.4822374 ,  0.5580728 ],
       [ 0.15772933,  0.48593197,  0.55801347],
       [ 0.15626973,  0.4896237 ,  0.557936  ],
       [ 0.15481488,  0.49331293,  0.55783967],
       [ 0.15336445,  0.49700003,  0.55772371],
       [ 0.1519182 ,  0.50068529,  0.55758733],
       [ 0.15047605,  0.50436904,  0.55742968],
       [ 0.14903918,  0.50805136,  0.5572505 ],
       [ 0.14760731,  0.51173263,  0.55704861],
       [ 0.14618026,  0.51541316,  0.55682271],
       [ 0.14475863,  0.51909319,  0.55657181],
       [ 0.14334327,  0.52277292,  0.55629491],
       [ 0.14193527,  0.52645254,  0.55599097],
       [ 0.14053599,  0.53013219,  0.55565893],
       [ 0.13914708,  0.53381201,  0.55529773],
       [ 0.13777048,  0.53749213,  0.55490625],
       [ 0.1364085 ,  0.54117264,  0.55448339],
       [ 0.13506561,  0.54485335,  0.55402906],
       [ 0.13374299,  0.54853458,  0.55354108],
       [ 0.13244401,  0.55221637,  0.55301828],
       [ 0.13117249,  0.55589872,  0.55245948],
       [ 0.1299327 ,  0.55958162,  0.55186354],
       [ 0.12872938,  0.56326503,  0.55122927],
       [ 0.12756771,  0.56694891,  0.55055551],
       [ 0.12645338,  0.57063316,  0.5498411 ],
       [ 0.12539383,  0.57431754,  0.54908564],
       [ 0.12439474,  0.57800205,  0.5482874 ],
       [ 0.12346281,  0.58168661,  0.54744498],
       [ 0.12260562,  0.58537105,  0.54655722],
       [ 0.12183122,  0.58905521,  0.54562298],
       [ 0.12114807,  0.59273889,  0.54464114],
       [ 0.12056501,  0.59642187,  0.54361058],
       [ 0.12009154,  0.60010387,  0.54253043],
       [ 0.11973756,  0.60378459,  0.54139999],
       [ 0.11951163,  0.60746388,  0.54021751],
       [ 0.11942341,  0.61114146,  0.53898192],
       [ 0.11948255,  0.61481702,  0.53769219],
       [ 0.11969858,  0.61849025,  0.53634733],
       [ 0.12008079,  0.62216081,  0.53494633],
       [ 0.12063824,  0.62582833,  0.53348834],
       [ 0.12137972,  0.62949242,  0.53197275],
       [ 0.12231244,  0.63315277,  0.53039808],
       [ 0.12344358,  0.63680899,  0.52876343],
       [ 0.12477953,  0.64046069,  0.52706792],
       [ 0.12632581,  0.64410744,  0.52531069],
       [ 0.12808703,  0.64774881,  0.52349092],
       [ 0.13006688,  0.65138436,  0.52160791],
       [ 0.13226797,  0.65501363,  0.51966086],
       [ 0.13469183,  0.65863619,  0.5176488 ],
       [ 0.13733921,  0.66225157,  0.51557101],
       [ 0.14020991,  0.66585927,  0.5134268 ],
       [ 0.14330291,  0.66945881,  0.51121549],
       [ 0.1466164 ,  0.67304968,  0.50893644],
       [ 0.15014782,  0.67663139,  0.5065889 ],
       [ 0.15389405,  0.68020343,  0.50417217],
       [ 0.15785146,  0.68376525,  0.50168574],
       [ 0.16201598,  0.68731632,  0.49912906],
       [ 0.1663832 ,  0.69085611,  0.49650163],
       [ 0.1709484 ,  0.69438405,  0.49380294],
       [ 0.17570671,  0.6978996 ,  0.49103252],
       [ 0.18065314,  0.70140222,  0.48818938],
       [ 0.18578266,  0.70489133,  0.48527326],
       [ 0.19109018,  0.70836635,  0.48228395],
       [ 0.19657063,  0.71182668,  0.47922108],
       [ 0.20221902,  0.71527175,  0.47608431],
       [ 0.20803045,  0.71870095,  0.4728733 ],
       [ 0.21400015,  0.72211371,  0.46958774],
       [ 0.22012381,  0.72550945,  0.46622638],
       [ 0.2263969 ,  0.72888753,  0.46278934],
       [ 0.23281498,  0.73224735,  0.45927675],
       [ 0.2393739 ,  0.73558828,  0.45568838],
       [ 0.24606968,  0.73890972,  0.45202405],
       [ 0.25289851,  0.74221104,  0.44828355],
       [ 0.25985676,  0.74549162,  0.44446673],
       [ 0.26694127,  0.74875084,  0.44057284],
       [ 0.27414922,  0.75198807,  0.4366009 ],
       [ 0.28147681,  0.75520266,  0.43255207],
       [ 0.28892102,  0.75839399,  0.42842626],
       [ 0.29647899,  0.76156142,  0.42422341],
       [ 0.30414796,  0.76470433,  0.41994346],
       [ 0.31192534,  0.76782207,  0.41558638],
       [ 0.3198086 ,  0.77091403,  0.41115215],
       [ 0.3277958 ,  0.77397953,  0.40664011],
       [ 0.33588539,  0.7770179 ,  0.40204917],
       [ 0.34407411,  0.78002855,  0.39738103],
       [ 0.35235985,  0.78301086,  0.39263579],
       [ 0.36074053,  0.78596419,  0.38781353],
       [ 0.3692142 ,  0.78888793,  0.38291438],
       [ 0.37777892,  0.79178146,  0.3779385 ],
       [ 0.38643282,  0.79464415,  0.37288606],
       [ 0.39517408,  0.79747541,  0.36775726],
       [ 0.40400101,  0.80027461,  0.36255223],
       [ 0.4129135 ,  0.80304099,  0.35726893],
       [ 0.42190813,  0.80577412,  0.35191009],
       [ 0.43098317,  0.80847343,  0.34647607],
       [ 0.44013691,  0.81113836,  0.3409673 ],
       [ 0.44936763,  0.81376835,  0.33538426],
       [ 0.45867362,  0.81636288,  0.32972749],
       [ 0.46805314,  0.81892143,  0.32399761],
       [ 0.47750446,  0.82144351,  0.31819529],
       [ 0.4870258 ,  0.82392862,  0.31232133],
       [ 0.49661536,  0.82637633,  0.30637661],
       [ 0.5062713 ,  0.82878621,  0.30036211],
       [ 0.51599182,  0.83115784,  0.29427888],
       [ 0.52577622,  0.83349064,  0.2881265 ],
       [ 0.5356211 ,  0.83578452,  0.28190832],
       [ 0.5455244 ,  0.83803918,  0.27562602],
       [ 0.55548397,  0.84025437,  0.26928147],
       [ 0.5654976 ,  0.8424299 ,  0.26287683],
       [ 0.57556297,  0.84456561,  0.25641457],
       [ 0.58567772,  0.84666139,  0.24989748],
       [ 0.59583934,  0.84871722,  0.24332878],
       [ 0.60604528,  0.8507331 ,  0.23671214],
       [ 0.61629283,  0.85270912,  0.23005179],
       [ 0.62657923,  0.85464543,  0.22335258],
       [ 0.63690157,  0.85654226,  0.21662012],
       [ 0.64725685,  0.85839991,  0.20986086],
       [ 0.65764197,  0.86021878,  0.20308229],
       [ 0.66805369,  0.86199932,  0.19629307],
       [ 0.67848868,  0.86374211,  0.18950326],
       [ 0.68894351,  0.86544779,  0.18272455],
       [ 0.69941463,  0.86711711,  0.17597055],
       [ 0.70989842,  0.86875092,  0.16925712],
       [ 0.72039115,  0.87035015,  0.16260273],
       [ 0.73088902,  0.87191584,  0.15602894],
       [ 0.74138803,  0.87344918,  0.14956101],
       [ 0.75188414,  0.87495143,  0.14322828],
       [ 0.76237342,  0.87642392,  0.13706449],
       [ 0.77285183,  0.87786808,  0.13110864],
       [ 0.78331535,  0.87928545,  0.12540538],
       [ 0.79375994,  0.88067763,  0.12000532],
       [ 0.80418159,  0.88204632,  0.11496505],
       [ 0.81457634,  0.88339329,  0.11034678],
       [ 0.82494028,  0.88472036,  0.10621724],
       [ 0.83526959,  0.88602943,  0.1026459 ],
       [ 0.84556056,  0.88732243,  0.09970219],
       [ 0.8558096 ,  0.88860134,  0.09745186],
       [ 0.86601325,  0.88986815,  0.09595277],
       [ 0.87616824,  0.89112487,  0.09525046],
       [ 0.88627146,  0.89237353,  0.09537439],
       [ 0.89632002,  0.89361614,  0.09633538],
       [ 0.90631121,  0.89485467,  0.09812496],
       [ 0.91624212,  0.89609127,  0.1007168 ],
       [ 0.92610579,  0.89732977,  0.10407067],
       [ 0.93590444,  0.8985704 ,  0.10813094],
       [ 0.94563626,  0.899815  ,  0.11283773],
       [ 0.95529972,  0.90106534,  0.11812832],
       [ 0.96489353,  0.90232311,  0.12394051],
       [ 0.97441665,  0.90358991,  0.13021494],
       [ 0.98386829,  0.90486726,  0.13689671],
       [ 0.99324789,  0.90615657,  0.1439362 ]]

test_cm = LinearSegmentedColormap.from_list(__file__, cm_data)

def make_idx(variable, range, nbins):
    
    bins = np.linspace(range[0], range[1], nbins+1)
    idx = np.searchsorted(bins, variable)
    
    if np.any(idx == 0) | np.any(idx == len(bins)):
        print 'Warning: there where values out of range.'
        exit()
        
    idx -= 1
    idx = idx.astype('int')
    
    offset = np.amin(idx)
    length = np.ptp(idx) + 1
    
    return idx, length, offset

def sigma_clip(values, niter=5):
    
    select = range(len(values))
    for i in range(niter):
        mean = np.mean(values[select])
        std = np.std(values[select])
        select = np.argwhere(np.abs(values-mean) < 2*std)
        
    return mean

def index_statistic(indices, values, statistic='mean', keeplength=True):
    
    if np.isscalar(indices):
        indices = np.array([indices])
        values = np.array([values])
    
    indices = indices.astype('int')
    
    known_statistics = ['mean', 'std', 'count', 'sum', 'median']
    if statistic not in known_statistics:
        raise ValueError('invalid statistic %r' % (statistic,))
    
    flatcount = np.bincount(indices, None)
    a = flatcount.nonzero()
    
    if statistic == 'mean':
        flatsum = np.bincount(indices, values)
        result = flatsum[a] / flatcount[a]
        
    elif statistic == 'std':
        flatsum = np.bincount(indices, values)
        flatsum2 = np.bincount(indices, values ** 2)
        result = np.sqrt(flatsum2[a] / flatcount[a] - (flatsum[a] / flatcount[a]) ** 2)
        
    elif statistic == 'count':
        result = flatcount[a]
        
    elif statistic == 'sum':
        flatsum = np.bincount(indices, values)
        result = flatsum[a]
        
    elif statistic == 'median':
        result = np.zeros(len(flatcount))
        for i in np.unique(indices):
            result[i] = np.median(values[indices == i])
        #result = result[a]
        
    if keeplength:
        return np.repeat(result, flatcount[a])
    else:
        return result

def sysrem(ind1, ind2, values, errors, a2=None, maxiter=50, eps=1e-3, idx=None):
    
    npoints = len(values)
    npars1 = len(np.unique(ind1))
    npars2 = len(np.unique(ind2))
    npars = npars1 + npars2
    
    if a2 is None:
        a2 = np.ones(np.amax(ind2)+1) # WRONG!!!
    else:
        fi = np.copy(a2)
    
    s = values/errors**2
    r = 1./errors**2
    t = np.bincount(ind2)
    
    niter = 0
    end_crit = True
    while (niter < maxiter) & (end_crit):
        
        print niter
        
        a1 = np.bincount(ind1, s*a2[ind2])/np.bincount(ind1, r*(a2**2)[ind2])
        a2 = (np.bincount(ind2, s*a1[ind1])+t/fi)/(np.bincount(ind2, r*(a1**2)[ind1])+t/fi**2)
        
        if niter == 0:
            end_crit = True
        else:
            end_crit = (np.nanmax(np.abs((a1o-a1)/a1o)) > eps) | (np.nanmax(np.abs((a2o-a2)/a2o)) > eps)
            print np.nanmax(np.abs((a1o-a1)/a1o)), np.nanmax(np.abs((a2o-a2)/a2o))
        
        a1o = np.copy(a1)
        a2o = np.copy(a2)
            
        niter += 1
    
    chisq = np.sum(r*(values-a1[ind1]*a2[ind2])**2)/(npoints-npars)
    chisq_pbin1 = np.bincount(ind1, r*(values-a1[ind1]*a2[ind2])**2)
    chisq_pbin2 = np.bincount(ind2, r*(values-a1[ind1]*a2[ind2])**2)
    
    return a1, a2, niter, chisq, chisq_pbin1, chisq_pbin2, npoints, npars
    
    
def transmission(filename):

    with h5py.File(filename) as f:
        
        hdr = f['table_header']
        ascc = hdr['ascc'].value
        vmag = hdr['vmag'].value
        ra = hdr['ra'].value
        dec = hdr['dec'].value
        Nobs = hdr['nobs'].value.astype('int')
        
        select = np.append(0, np.cumsum(Nobs))
    
        lst = np.zeros(np.sum(Nobs))
        flux0 = np.zeros(np.sum(Nobs))
        eflux0 = np.zeros(np.sum(Nobs))
        sky = np.zeros(np.sum(Nobs))
        flags = np.zeros(np.sum(Nobs))
        x = np.zeros(np.sum(Nobs))
        y = np.zeros(np.sum(Nobs))
        
        data = f['data']
        for i in range(len(ascc)):
            lst[select[i]:select[i+1]] = data[ascc[i]]['lst']
            flux0[select[i]:select[i+1]] = data[ascc[i]]['flux0']
            #eflux0[select[i]:select[i+1]] = data[ascc[i]]['eflux0']
            eflux0[select[i]:select[i+1]] = index_statistic(data[ascc[i]]['lstidx']//50, data[ascc[i]]['flux0'], statistic='std')
            sky[select[i]:select[i+1]] = data[ascc[i]]['sky']
            flags[select[i]:select[i+1]] = data[ascc[i]]['flag']
            x[select[i]:select[i+1]] = data[ascc[i]]['x']
            y[select[i]:select[i+1]] = data[ascc[i]]['y']
    
    dec_1 = np.copy(dec)
    
    print 'Done reading.'

    # Sky coordinates of data.
    ha = np.mod(lst*15.-np.repeat(ra,Nobs), 360.)
    dec = np.repeat(dec, Nobs)

    # Create the indices.
    star_id = np.repeat(np.arange(len(ascc)), Nobs)
    
    # Remove bad data.
    here, = np.where((flags < 1)&(flux0>0)&(sky>0)&(eflux0>0))
    flux0 = flux0[here]
    eflux0 = eflux0[here]
    ha = ha[here]
    dec = dec[here]
    x = x[here]
    y = y[here]
    star_id = star_id[here]

    hg = CartesianGrid(800, 600, margin=50)
    bins, binnum = hg.find_gridpoint(x, y, compact=True)
    
    count = index_statistic(binnum, flux0, statistic='count', keeplength=False)
    print np.percentile(count, 5), np.percentile(count, 95)
    
    # Compute the transmission using sysrem.
    idx = np.digitize(dec_1, np.linspace(-90,90,721))
    a2, a1, niter, chisq, chisq_pbin2, chisq_pbin1, npoints, npars = sysrem(binnum, star_id, flux0, eflux0, a2=1e7*10**(vmag/-2.5))
    
    trans = hg.put_values_on_grid(a2, ind=bins, fill_value=np.nan)
    count = hg.put_values_on_grid(count, ind=bins)
    chisq_map = hg.put_values_on_grid(chisq_pbin2, ind=bins, fill_value=np.nan)
    stars = np.unique(star_id)
    
    plt.subplot(211)
    trans = trans/np.percentile(trans[1:-1,1:-1], 99.9)
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', origin='lower', vmin=0, vmax=1, extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label('Transmission')
    plt.xlim(0, 4008)
    plt.ylim(0, 2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(212)
    arr = a1[stars]/(1e7*10**(vmag[stars]/-2.5))
    plt.plot(dec_1[stars], arr, '.', alpha=.1)
    plt.xlim(np.amin(dec_1), np.amax(dec_1))
    plt.ylim(.8, 1.2)
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    plt.show()
    
    exit()
    median = index_statistic(idx[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), statistic='median', keeplength=False)
    
    plt.subplot(211)
    plt.semilogy(dec_1[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), '.', alpha=.1)
    plt.xlim(np.amin(dec), np.amax(dec))
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    a1[stars] = a1[stars]/median[idx[stars]]
    
    plt.subplot(212)
    plt.semilogy(dec_1[stars], a1[stars]/(1e7*10**(vmag[stars]/-2.5)), '.', alpha=.1)
    plt.xlim(np.amin(dec), np.amax(dec))
    plt.xlabel('Dec [deg]')
    plt.ylabel('F/V')
    
    plt.tight_layout()
    plt.show()
    plt.close()
    
    trans[:,np.unique(idx[stars])] = trans[:,np.unique(idx[stars])]*median[np.unique(idx[stars])]
    
    ax = plt.subplot(221)
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', vmin=0, vmax=1.5, extent=(0,360,-90,90), cmap=test_cm)
    plt.colorbar().set_label('Transmission')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(count[1:-1,1:-1].T, interpolation='None', aspect='auto', origin='lower', extent=(0,360,-90,90), cmap=test_cm)
    plt.colorbar().set_label('# points')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.subplot(223, sharex=ax, sharey=ax)
    arr = chisq_map[1:-1,1:-1].T/count[1:-1,1:-1].T
    plt.imshow(arr, interpolation='None', aspect='auto', origin='lower', vmin=0, vmax=np.percentile(arr[np.isfinite(arr)], 98), extent=(0,360,-90,90), cmap=test_cm)
    plt.colorbar().set_label(r'$\chi^2/N$')
    plt.xlim(np.amin(ha), np.amax(ha))
    plt.ylim(np.amin(dec), np.amax(dec))
    plt.xlabel('HA [deg]')
    plt.ylabel('Dec [deg]')
    
    plt.subplot(224)
    arr = chisq_pbin1[stars]/np.bincount(star_id)[stars]
    plt.plot(vmag[stars], arr, '.', alpha=.1)
    plt.xlim(8.4, 2)
    plt.ylim(0, np.percentile(arr[np.isfinite(arr)], 98))
    plt.xlabel('V [mag]')
    plt.ylabel(r'$\chi^2/N$')
    
    plt.tight_layout()
    #plt.savefig('0203LPE_maps.png')
    plt.show()
    plt.close()
    
    cg = CartesianGrid(200, 150, margin=50)
    binnum = cg.find_gridpoint(x, y)
    a2 = np.bincount(binnum, flux0*a1[star_id]/eflux0**2)/np.bincount(binnum, (a1**2)[star_id]/eflux0**2)
    trans = cg.put_values_on_grid(a2)
    
    count = np.bincount(binnum)
    count = cg.put_values_on_grid(count)
    
    chisq_map = np.bincount(binnum, (flux0-a1[star_id]*a2[binnum])**2/eflux0**2)
    chisq_map = cg.put_values_on_grid(chisq_map)
    
    chisq_stars = np.bincount(binnum, (flux0-a1[star_id]*a2[binnum])**2/eflux0**2)
    
    ax = plt.subplot(221)
    trans = trans/np.percentile(trans[1:-1,1:-1], 99)
    plt.imshow(trans[1:-1,1:-1].T, interpolation='None', origin='lower', aspect='equal', vmin=0, vmax=1, extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label('Transmission')
    plt.xlim(0,4008)
    plt.ylim(0,2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(222, sharex=ax, sharey=ax)
    plt.imshow(count[1:-1,1:-1].T, interpolation='None', origin='lower', aspect='equal', extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label('# points')
    plt.xlim(0,4008)
    plt.ylim(0,2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(223, sharex=ax, sharey=ax)
    arr = chisq_map[1:-1,1:-1].T/count[1:-1,1:-1].T
    plt.imshow(arr, interpolation='None', origin='lower', aspect='equal', vmin=0, vmax=np.percentile(arr[np.isfinite(arr)], 98), extent=(50,4008-50,50,2672-50), cmap=test_cm)
    plt.colorbar().set_label(r'$\chi^2/N$')
    plt.xlim(0,4008)
    plt.ylim(0,2672)
    plt.xlabel('x [pix]')
    plt.ylabel('y [pix]')
    
    plt.subplot(224)
    arr = chisq_stars[stars]/np.bincount(star_id)[stars]
    plt.plot(vmag[stars], arr, '.', alpha=.1)
    plt.xlim(8.4, 2)
    plt.ylim(0, np.percentile(arr[np.isfinite(arr)], 98))
    plt.xlabel('V [mag]')
    plt.ylabel(r'$\chi^2/N$')
    
    plt.tight_layout()
    #plt.savefig('0203LPE_maps.png')
    plt.show()
    
    
    
    
    return 0

if __name__ == '__main__':
    #import argparse
    
    #parser = argparse.ArgumentParser(description='Compute the transmission map for a particular night and camera.')
    
    #parser.add_argument('path', type=str, help='the global path for search')
    #parser.add_argument('-n', '--night', default='', type=str, help='the night to reduce')
    #parser.add_argument('-c', '--camera',default='', type=str, help='the camera to reduce')
    
    #args = parser.parse_args()
    
    filelist = glob.glob('/data2/mascara/LaPalma/20150203LPE/fLC/fLC_*.hdf5')
    #filelist = np.sort(filelist)

    #filelist = ['/data2/talens/fLC_20150203LPC.hdf5']

    for filename in filelist:
        print 'Data:', filename
        transmission(filename)
