from scipy.special import comb,binom
import numpy as np
from addressing import *



class RASAddrEngine():
    """Compute the address of given electronic configuration or
    the electron configuration of given address

    Attributes:
        NOrb: int, number of total orbitals
        NElec: int, number of total electrons
        GHF: boolean, generalized Hatree-Fock or restricted Hatree-Fock

        AddrArray: list of list: defining addressing array for each
                   catogories, weight by holes (RAS1) or electrons (RAS2/3)

        CatOffs: list of list, contains informations of configuation
                  catogories
            [Cat_i][X]
                X:  0 ... RAS1 Offset
                    1 ... RAS2 Offset
                    2 ... RAS3 Offset
                    3 ... Total Configs in Cat_i
                    4 ... Cat_i Offset, if equals -1 -> not valid
    """

    def __init__(self, NOrb, NElec,
                 MxHole, NORAS1, MxElec, NORAS3,
                 GHF=True, **kwargs):
        """Constructor of the Engine
        """
        assert MxElec <= NElec
        assert NORAS1 >= MxHole
        assert NORAS3 >= MxElec
        assert NORAS1 <= NElec
        assert NORAS1 + NORAS3 <= NOrb

        self.GHF    = GHF
        self.NOrb   = NOrb
        self.NElec  = NElec
        self.MxHole = MxHole
        self.MxElec = MxElec
        self.NORAS  = [NORAS1, NOrb-NORAS1-NORAS3, NORAS3]
        self.NCat   = (MxHole+1)*(MxElec+1)

        # calculate addressing arrays that might be used
        self._initAddrArray()
        self._initCategory()

        return

    def _initCategory(self):
        """
        """
        if self.GHF:
            self.CatOffs = np.zeros(shape=(self.NCat,5),dtype=np.int)
            self.NConfigs = 0

            i_cat = 0
            for i_e in range(self.MxElec+1):
                for i_h in range(self.MxHole+1):
                    NERAS2 = self.NElec - self.NORAS[0] + i_h - i_e
                    if NERAS2 >= 0 and NERAS2 <= self.NORAS[1]:
                        self.CatOffs[i_cat][0]= 1
                        self.CatOffs[i_cat][1]= binom(self.NORAS[0],i_h)
                        self.CatOffs[i_cat][2]= binom(self.NORAS[1],NERAS2)*self.CatOffs[i_cat][1]
                        self.CatOffs[i_cat][3]= binom(self.NORAS[2],i_e)*self.CatOffs[i_cat][2]

                    if self.CatOffs[i_cat][3] != 0:
                        self.CatOffs[i_cat][4] = self.NConfigs
                        self.NConfigs += self.CatOffs[i_cat][3]
                    else:
                        self.CatOffs[i_cat][4] = -1

                    i_cat +=1

        else:
            raise Exception('No implementation!')
        return

    def _initAddrArray(self):
        """ initialize the addressing arrays"""
        self.AddrArray = []
        self.deAddrArray = []

        if self.GHF:
            # RAS1
            addrarry_ras = []
            deaddrarry_ras = []
            for i in range(self.MxHole+1):
                addrarry_ras.append(addressing_array(i,self.NORAS[0]))
#                deaddrarry_ras.append(de_addressing_array(addrarry_ras[-1]))
            self.AddrArray.append(addrarry_ras)
            self.deAddrArray.append(deaddrarry_ras)

            # RAS2
            addrarry_ras = []
            deaddrarry_ras = []
            for i in range(self.NORAS[1]):
                addrarry_ras.append(addressing_array(i,self.NORAS[1]))
#                deaddrarry_ras.append(de_addressing_array(addrarry_ras[-1]))
            self.AddrArray.append(addrarry_ras)
            self.deAddrArray.append(deaddrarry_ras)

            # RAS3
            addrarry_ras = []
            deaddrarry_ras = []
            for i in range(self.MxElec+1):
                addrarry_ras.append(addressing_array(i,self.NORAS[2]))
#                deaddrarry_ras.append(de_addressing_array(addrarry_ras[-1]))
            self.AddrArray.append(addrarry_ras)
            self.deAddrArray.append(deaddrarry_ras)

        else:
            raise Exception('No implementation!')

        return

    def addressing(self,elec_config):
        """Compute the address of electron configurations

        Args:
            elec_config: strings of 0 and 1 with length of NOrb

        Returns:
            config_addr: configuraiton address, weight by electron

        """

        assert set(elec_config) == {'0', '1'}
        assert len(elec_config) == self.NOrb
        assert elec_config.count('1') == self.NElec


        ras_config = [elec_config[:self.NORAS[0]],
                      elec_config[self.NORAS[0]:self.NORAS[0]+self.NORAS[1]],
                      elec_config[self.NORAS[0]+self.NORAS[1]:sum(self.NORAS)]]
        ras_occ    = [ras_config[0].count('0'),
                      ras_config[1].count('1'),
                      ras_config[2].count('1'),]

        #complimentary RAS1 configs
        ras1_config = []
        for i in range(self.NORAS[0]):
            ras1_config.append(1 - int(ras_config[0][i]))
        ras_config[0] = ''.join(list(map(str,ras1_config)))

        #print(ras_occ)

        assert ras_occ[0] <= self.MxHole
        assert ras_occ[2] <= self.MxElec

        ras_addr = np.zeros(3, dtype=np.int)

        for X in range(0,3):
            if ras_occ[X] == 0 or ras_occ[X] == self.NORAS[X]:
                ras_addr[X] = 1
            else:
                ras_addr[X] = addressing_single_graph(
                    ras_config[X], self.AddrArray[X][ras_occ[X]])

        i_cat =  ras_occ[0] + (self.MxHole+1)*ras_occ[2]

        #complimentary RAS1 address
        ras_addr[0] = self.CatOffs[i_cat][1] - ras_addr[0] + 1

        #print(ras_addr)
        config_addr = np.dot(ras_addr-1, self.CatOffs[i_cat][:3]) \
            + 1 + self.CatOffs[i_cat][4]

        return config_addr

    def de_addressing(self,config_addr, join_ras=False):
        """Compute the address of electron configurations

        Args:
            config_addr: configuraiton address, weight by electron

        Returns:
            elec_config: strings of 0 and 1 with length of NOrb
        """

        assert isinstance(config_addr, int)
        assert config_addr > 0 and config_addr <= self.NConfigs

        i_cat_prev = 0
        for i_cat in range(len(self.CatOffs)):
            if self.CatOffs[i_cat][4] == -1: continue

            if self.CatOffs[i_cat][4] >= config_addr:
                i_cat = i_cat_prev
                break

            i_cat_prev = i_cat

        ras_occ = np.zeros(3, dtype=np.int)
        ras_occ[2] = int(i_cat/(self.MxHole+1))
        ras_occ[0] = i_cat - ras_occ[2]*(self.MxHole+1)
        ras_occ[1] = self.NElec - self.NORAS[0] + ras_occ[0] - ras_occ[2]

        #print(ras_occ)

        cat_addr = config_addr - self.CatOffs[i_cat][4]

        ras_addr = np.zeros(3, dtype=int)
        cat_addr -= 1

        for X in range(2, -1, -1):
            # print(self.CatOffs[i_cat])
            ras_addr[X] = int(cat_addr / self.CatOffs[i_cat][X])
            cat_addr -= ras_addr[X]*self.CatOffs[i_cat][X]
            ras_addr[X] +=1

        #complimentary RAS1 address
        ras_addr[0] = self.CatOffs[i_cat][1] - ras_addr[0] + 1
        #print(ras_addr)

        ras_config = []
        ras_occ_type = [0, 1, 1]
        for X in range(0,3):
            if ras_occ[X] == 0:
                ras_config.append(str(1-ras_occ_type[X])*self.NORAS[X])
            elif ras_occ[X] == self.NORAS[X]:
                ras_config.append(str(ras_occ_type[X])*self.NORAS[X])
            else:
                #print(ras_addr[X], type(ras_addr[X]))
                ras_config.append(de_addressing_single_graph(
                    ras_addr[X], self.AddrArray[X][ras_occ[X]]))#, 
                    #self.deAddrArray[X][ras_occ[X]]))

        #complimentary RAS1 configs
        if(self.CatOffs[i_cat][1] != 1):
            ras1_config = []
            for i in range(self.NORAS[0]):
                ras1_config.append(1 - int(ras_config[0][i]))
            ras_config[0] = ''.join(list(map(str,ras1_config)))

        if not join_ras:
            return ras_config
        else:
            return ''.join(ras_config)
