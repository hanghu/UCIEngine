"""
   A Module of utility functions to analyze 
   analyaze Configuration Interaction Strings 
"""

import numpy as np
from scipy.special import comb

def compute_NConf_GHF_CAS(no,ne):
    return int(comb(no,ne))

def compute_NConf_RHF_CAS(no,ne):
    return int(comb(no,int(ne/2)))**2

def compute_NConf_GHF_GAS_Cat(nos,nes):
    assert len(nes) == len(nos)
    return np.prod([int(comb(nos[i], nes[i])) for i in range(len(nes))])

def compute_NConf_GHF_RAS(mh, me, nos, ne):
    assert len(nos) == 3
    NConfs = 0
    for ih in range(0,mh+1):
        for ie in range(0,me+1):
            NConfs += compute_NConf_GHF_GAS_Cat(nos, [nos[0]-ih,ne-nos[0]+ih-ie,ie])
    return NConfs

def compute_NConf_RHF_RAS(mh, me, nos, ne_a, ne_b):
    assert len(nos) == 3
    NConfs = 0
    for ih_a in range(0,mh+1):
        for ie_a in range(0,me+1):
            for ih_b in range(0,mh+1-ih_a):
                for ie_b in range(0,me+1-ie_a):
                    NConfs += (compute_NConf_GHF_GAS_Cat(nos, [nos[0]-ih_a,ne_a-nos[0]+ih_a-ie_a,ie_a])*  
                               compute_NConf_GHF_GAS_Cat(nos, [nos[0]-ih_b,ne_b-nos[0]+ih_b-ie_b,ie_b]))
   
    return NConfs

