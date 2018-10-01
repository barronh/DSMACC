"""
This file is an example of how to use dsmacc in python
to modify reaction rates by applying a scalar.

This only works with kpp_lsode, but the default integrator
is rosenbrock for cri and kpp_lsode for geoschem

In this simple test, the JPL (15-10) uncertainty equation (Eq 1) is used
to modify the "O3 + NO --> NO2 + O2" rxn

f(T) = f(298K)exp|g(1/T - 1/298)|       Eq 1

"""
from dsmacc import model
from scipy.constants import Avogadro
import numpy as np
class jpluncertainty(model):
    def __init__(self, *args, rxnname = None, f298K = 1., g = 0., pow = 0, **kwds):
        """
        Description
        -----------
        Updates rate constant (RCONST) values by applying JPL
        uncertainty
        Parameters
        ----------
        args and kwds : positional arguments and keywords are passed 
                        to the standard model (see model.__doc__)
        
        rxnname : name excluding preceding or trailing spaces for reaction
        f298K   : JPL uncertainty estimate at 298K
        g       : JPL temperature dependent uncertainty factor
        pow     : power to apply to unceratinty

        Notes
        -----
        fT = f298Kexp|g(1/T - 1/298)|       Eq 1
        
        RCONST[rxnidx] *= fT**pow
        """
        super(jpluncertainty, self).__init__(*args, **kwds)
        self.rxnname = rxnname
        self.rxnidx = self.eqn_names.index(rxnname)
        self.f298K = f298K
        self.g = g
        self.pow = pow
     
    def custom_after_rconst(self):
        pyglob = self.pyglob
        scale = self.f298K * np.exp(self.g * np.abs(1 / pyglob.TEMP - 1 / 298.))
        scale = scale**(self.pow)
        self.pyglob.rconst[self.rxnidx] *= scale

from dsmacc import model, dynenv
from scipy.constants import Avogadro as N_A, R

# Choose model (cri or geoschem)
modelname = 'geoschem'

# Set some global variables to intialize
globvar = {'LAT': 30.0, 'LON': -90.0,
           'TEMP': 288.0, 'PRESS': 101325.0,
           'JDAY_GMT': 2006186.0}

# Add water in #/cm3
exec('H2O = 9.50e-9 * N_A * PRESS / TEMP / R', globals(), globvar)

# Set the starting date
startdate = globvar['JDAY_GMT']

# Initial conditions are taken from ../cri/cri_init.dat
initcond = {
            'O3': 41.3038, 'CO': 142.2, 'SO2': 13.63,
            'NO': 1.7879, 'NO2': 5.9328, 'NO3': 0.0002,
            'N2O5': 0.0028, 'HNO3': 0.754, 'PAN': 0.6496,
            'PPN': 0.0825,
            'HCHO': 2.7141, 'CH4': 1879., 'CH3OH': 6.1246,
            'HCOOH': 0.1881, 'C2H4': 1.458, 'C2H6': 8.658,
            'C3H8': 6.59, 'CH3CHO': 1.8779,
            'CH3COCH3': 2.3931, 'NC4H10': 7.171, 'MEK': 0.7013,
            }

if modelname == 'geoschem':
    """
    requires linked FAST-JX
    """
    cri2gc = {'CH3CHO': 'ALD2',
              'CH3COCH3': 'ACET',
              'NC4H10': 'ALK4',
              'HCHO': 'CH2O',
              'CH3OH': 'MOH'}
    initcond = {cri2gc.get(k, k): v for k, v in initcond.items()}
    initcond['O2'] = 0.2095e9
    initcond['N2'] = 0.7905e9
    del initcond['C2H4']

for pow in [-2, -1, 0, 1, 2]:
    prefix = ('O3plNO_NO2plO2_k_pos%ssigma' % str(pow)).replace('pos-', 'neg')
    print(prefix)
    mod = jpluncertainty(modelname = modelname,
                    outconcpath = prefix + 'conc.dat',
                    outratepath = prefix + 'rate.dat',
                    rxnname = 'NO + O3 --> NO2 + O2', f298K = 1.1, g = 200, pow = pow)

    mod.run(startdate, 24, 180, conc_ppb = initcond, globvar = globvar)
    
    del mod

import pandas as pd
mn2 = pd.read_csv('O3plNO_NO2plO2_k_neg2sigmaconc.dat', index_col = 'JDAY_GMT')
mn1 = pd.read_csv('O3plNO_NO2plO2_k_neg1sigmaconc.dat', index_col = 'JDAY_GMT')
orig = pd.read_csv('O3plNO_NO2plO2_k_pos0sigmaconc.dat', index_col = 'JDAY_GMT')
pl1 = pd.read_csv('O3plNO_NO2plO2_k_pos1sigmaconc.dat', index_col = 'JDAY_GMT')
pl2 = pd.read_csv('O3plNO_NO2plO2_k_pos2sigmaconc.dat', index_col = 'JDAY_GMT')
paired = orig.join(pl1, lsuffix = '', rsuffix = '+1s')\
              .join(pl2, rsuffix = '+2s')\
              .join(mn1, rsuffix = '-1s')\
              .join(mn2, rsuffix = '-2s')
ratios = (pl1 / orig)
paired.plot(y = ['O3-2s', 'O3-1s', 'O3', 'O3+1s', 'O3+2s'], color = ['darkblue', 'lightblue', 'k', 'pink', 'red']).figure.savefig('ratesens.png')
print(ratios[::60].filter(['O3']).to_csv(float_format = '%.6f'))
