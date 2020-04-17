"""
  A module for CI String Addressing
"""

from scipy.special import comb
import numpy as np
import re

def addressing_array(ne,no):
    """Generate and return the addrssing array for given number of
    orbitals no and number of electron ne

    Args:
        ne: int, number of orbitals
        no: int, number of electrons

    Returns:
        Z: numpy 2d array with shape (ne,no), addressing array.
    """
    assert no >= ne and ne >= 0

    if ne == 0: return None

    Z = np.zeros(shape=(ne,no), dtype=np.int)
    for k in range(1,ne):
        for l in range(k,no-ne+k+1):
            for m in range(no-l+1, no-k+1):
                #Z[k-1][l-1] += binom(m, ne-k) - binom(m-1, ne-k-1)
                Z[k-1][l-1] += comb(m-1, ne-k)
    for l in range(ne,no+1):
        Z[ne-1][l-1] = l - ne

    return Z

def addressing_single_graph(config, Z):
    """
       addressing strings within single CAS-like graph
       (with all possible combinations of electron occupying orbitals)
    """
    assert isinstance(config,str)
    assert re.search(r'[^01]',config) is None
    assert len(config) == Z.shape[1]

    addr = 1
    ie = 0
    for io in range(len(config)):
        ie += int(config[io])
        if config[io] == '1':
            addr += Z[ie-1][io]

    return addr

def de_addressing_array(Z):
    """Generate and return the deaddrssing array 
       for an addressing array 
    """
    assert len(Z.shape) == 2
    
    if(Z.shape[0] == 1): return Z.copy()
    
    Zd = np.zeros(Z.shape, dtype=np.int)
    Zd[-1] = Z[-1]
    for i in range(-2,-Z.shape[0]-1,-1):
        Zd[i] = np.roll(Zd[i+1],-1)+Z[i]

    return Zd     

def de_addressing_single_graph(addr, Z, Zd=None):
    """
        deaddressing address within single CAS-like graph
        (with all possible combinations of electron occupying orbitals)
    """
    assert isinstance(addr,int) or isinstance(addr,np.int64)
    assert addr > 0 
    if Zd is None:
        Zd = de_addressing_array(Z)
    else:
        assert Z.shape == Zd.shape

    config = ['0']*Z.shape[1]
    
    io_search = [Z.shape[1]-Z.shape[0], -1] # great to small
    addr_p  = addr - 1

    for ie in range(Z.shape[0]):
        for io in range(*io_search,-1):
            if addr_p >= Zd[ie][io]:
                config[io] = '1'
                addr_p -= Z[ie][io] 
                io_search[0] +=1
                io_search[1] = io
                break
    
    return "".join(config)























