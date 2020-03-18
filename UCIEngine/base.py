"""
  A python module to implement vectorized universal CI Enegine
  based on generalized active spce (GAS) concept
"""
import numpy as np
from utils import compute_NConf_GHF_GAS_Cat

class UCIEngine(object):
    """base class of Universal Configuration Interaction Engine
       Note this base class is GHF-like implementation. 
    
    Attributes:
        NOs : list of int; number of orbitals in each space
        NO  : int; total number of orbials in each space
        NE  : int; total number of electrons in each space
        NCnf: int; total number of configurations in the engine
        NGAS: int, number of generalized active spaces
        NCat: int, number of categories
        
        Cats: dict; information of categorires
              keys:
                   NEs: list of list of int; number of electrons in each space 
                                             in each category
                   NCnfs: list of int; number of configurations in each category
                   PSt : list of int; starting pointer of each category 
                   graph_1e: adjacency list; graph representation of 1e interactions  

        Z_cache: dict with keys as (int,int) and values as np.matrix; 
                 Addressing array cache
                 key: (number of electron, number of orbitals)
                 values: addressing array
        Zd_cache: dict with keys as (int,int) and values as np.matrix; 
                  Addressing array cache
                  key: (number of electron, number of orbitals)
                  values: deaddressing array

    """
    def __init__(self, NOs=None, NE=None, init_by=None,
                 *args, **kwargs):
        """
        Args:    
            NOs: list of int, number of orbtals in each space

            NE: int, total number of electron in space

            init_by: None,
                     Categories, 
                     GAS Specfication,
                     Configuraitons,

            None: init an empty class and 
                  latter assign catergories/configuraitons by hand
                  Additional Args: None

            Categories: occupations of each categroies
                  Additional Args: 
                     Cats_occ: list of list of int, electrons occupancy in each space

            GAS Specificaiton: [TODO]
                  Additional Args: default is equivalent to FCI in defined space
                             Cat_ref:  list of int, electron occupancy in reference categries, 
                                       default as HF
                      GAS_excitations: String,  
                                       'F' (default) -- all possible excitations
                                       'S' -- all possible single excitations
                                       'D' -- all possible double excitations
                                       'T' -- all possible triple excitations
                                       'Q' -- all possible quadruple excitation
                                       
                                       or list of tuples, explicit excitations between space

                      Occ_restriction: list of tuples, min and max occupation in each space
                                       default as all possible occupation

            Configurations: explicit configurations, automatic vectroiziation [TODO]
                  Additional Args: 
                     Configs: list of strings as  explicit configurations    
                   
        """
        if init_by is None: return

        if init_by in ['Categories','categroies','Cat','cat']:
            
            NEs = kwargs.get('Cats_occ')
            assert NEs is not None 
            assert NOs is not None
            if NE is None: NE = sum(NEs[0])
            
            self.NOs  = NOs
            self.NO   = sum(NOs)
            self.NE   = NE
            self.NGAS = len(NOs)

            assert self.NO > self.NE
            
            for i in range(len(NEs)):    
                assert len(NEs[i]) == len(NOs)
                assert sum(NEs[i]) == NE
        
        else:
            raise TypeError('Wrong initialization method')
        
        
        # init categories

        self._init_categories(NEs)

        return 
        
    def _init_categories(self,NEs):
        """ calculate information of catergories, 
            Note that assuming there is no redundency in NEs"""
        self.Cats = {}
        #TODO: add some sorting for categories
        self.Cats['NEs']      = []
        self.Cats['NCnfs']    = []
        self.Cats['PSt']      = []
        
        PSt = 0
        for i in range(len(NEs)):
            if NEs[i] in self.Cats['NEs']:
                print('Warning! The %s_th cat is redundant' %i)
                continue

            for j in range(self.NGAS):
                if NEs[i][j] > self.NOs[j]:
                    print('Warning! The %s_th cat is not valid' %i)
                break

            self.Cats['NEs'].append(NEs[i])
            self.Cats['NCnfs'].append(compute_NConf_GHF_GAS_Cat(self.NOs,NEs[i]))
            self.Cats['PSt'].append(PSt)
            PSt += self.Cats['NCnfs'][-1]
        
        self.NCnf = sum(self.Cats['NCnfs'])
        self.NCat = len(self.Cats['NCnfs'])

        # generate 1e interaction graph
        graph_1e = {}
        valid_NEs = np.array(self.Cats['NEs'])
        for i in range(self.NCat):
            graph_1e[i] = {i:[]}
            for x in range(self.NGAS):
                if(valid_NEs[i][x] != 0): 
                    graph_1e[i][i].append((x,x))
        
        for i in range(self.NCat):     
            for j in range(i+1,self.NCat):
                Ex_level = valid_NEs[j] - valid_NEs[i]
                if sum(abs(Ex_level)) == 2:
                    x_gain = np.argwhere(Ex_level== 1)[0][0] 
                    x_lose = np.argwhere(Ex_level== -1)[0][0]

                    graph_1e[i][j] = (x_lose, x_gain)
                    graph_1e[j][i] = (x_gain, x_lose)

        self.Cats['graph_1e'] = graph_1e

        return 












        

