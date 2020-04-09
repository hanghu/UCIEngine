import re
from gLog_tools import match_partial_line
from RAS_addressing import RASAddrEngine

class RASCIStates:
    """
    Attributes:
        Energies: list
        Coeffs: list of dictionary
        CAS_spec: specification of CAS
        RAS_spec: specificaiton of RAS
    """
    def __init__(self, filename, maxNCoeff=50):

        self._readStates(filename, maxNCoeff)

        return

    def _readStates(self,filename, maxNCoeff):
        """
        Read RAS configuration and State info
        """
        f = open(filename,'r')
        iterCoeff = -1
        i_na = 0
        self.coeffs = []
        self.energies = []
        self.CAS_spec = None
        self.RAS_spec = None
        self.complex = True

        for line in f:
            if match_partial_line(' RAS(',line): 
                info = list(filter(None, re.split('[(),\s]', line)))[1:]
                self.RAS_spec = list(map(int, info))

            elif match_partial_line(' CAS(',line):
                info = list(filter(None, re.split('[(),\s]', line)))[1:]
                self.CAS_spec = list(map(int, info))

            elif match_partial_line(' State:',line):
                self.energies.append(float(line.split()[-1]))
                iterCoeff = 0
                coeffs_i = {}

            elif iterCoeff >= 0:
                coeffs = list(filter(None, re.split('[(),\s]',line)))
                
                if(iterCoeff == 0):
                    try:
                        int(coeffs[3])
                    except:
                        self.complex = False
                
                if(self.complex):
                    n_coeffs = int(len(coeffs)/3)
                else:
                    n_coeffs = int(len(coeffs)/2)

                n_left_to_record = maxNCoeff - iterCoeff

                if n_left_to_record >= n_coeffs:
                    n_record = n_coeffs
                else:
                    n_record = n_left_to_record

                for j in range(n_record):
                    jj = j*3 if self.complex else j*2

                    if coeffs[jj][0] == '*':
                        addr = 'na' + str(i_na)
                        i_na += 1
                    else:
                        addr = int(coeffs[jj])

                    coeffs_i[addr] = float(coeffs[jj+1])
                    coeffs_i[addr] += 1j*float(coeffs[jj+2]) if self.complex else 0

                iterCoeff += n_record

                if len(coeffs) == 0 or iterCoeff == maxNCoeff:
                    iterCoeff = -1
                    self.coeffs.append(coeffs_i)
            else:
                pass

        f.close()
        return

    def gen_state_config_report(self, threshold = 0.1, orb_assign = None):

        if orb_assign is not None: assert isinstance(orb_assign, dict)
        if self.RAS_spec is None: self.RAS_spec = [0, 0, 0, 0]

        print('N_TOTAL_ACTIVE_ORBITALS = %i, N_TOTAL_ACTIVE_ELECTRON = %i'
              %(self.CAS_spec[1], self.CAS_spec[0]))
        print('N_MAX_HOLE = %i, N_RAS1 = %i, N_MAX_ELECTRON = %i, N_RAS3 = %i' % tuple(self.RAS_spec))
        print(' ')
        print('Analyzing leading determinants with squared oefficient absolute value larger than %f' %threshold)
        print(' ')

        ENGINE = RASAddrEngine(NOrb=self.CAS_spec[1], NElec=self.CAS_spec[0],
                                MxHole=self.RAS_spec[0], NORAS1=self.RAS_spec[1],
                                MxElec=self.RAS_spec[2], NORAS3=self.RAS_spec[3])

        for i in range(len(self.coeffs)):
            print('State %i, Energy: %20.10f' %(i+1, self.energies[i]))
            for addr, coeff in self.coeffs[i].items():
                norm_sqaure = abs(coeff**2)
                if(norm_sqaure < threshold): continue

                if isinstance(addr, int):

                    elec_config = ENGINE.de_addressing(addr, join_ras=True)

                    if orb_assign is None:
                        print_elec_config = str(elec_config)
                    else:
                        print_elec_config = ''
                        for orb_name, orb_range in orb_assign.items():
                            elec_occ = sum([int(elec_config[x-1]) for x in orb_range])
                            print_elec_config += ' %s%2i' %(orb_name, elec_occ)

                    print('  |%8i>: %8.5f+%8.5fi; |C|^2: %.5E; ' %
                          (addr, coeff.real, coeff.imag, norm_sqaure) + print_elec_config)

                else:
                    print('  |%8i>: %8.5f+%8.5fi; |C|^2: %.5E ' %
                          (addr, coeff.real, coeff.imag, norm_sqaure))


#---------------------------------------------
# Running as a Script
#---------------------------------------------

if __name__ == '__main__':
    filename = 'rasci_example.log'
    rasci = RASCIStates('rasci_example.log')
    list_5f =  list(range(28,41))
    list_5f.append(7)

    rasci.gen_state_config_report(threshold=0.01,
        orb_assign={'U-F Bonding':range(1,7),
                    'U(5f)': list_5f,
                    'U(7s8s)': range(8,12),
                    'U(6d)': range(12,22),
                    'U(7p)': range(22,28),
                    'U-F Anti-bonding':range(41,47)})
