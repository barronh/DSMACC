"""
This file is an example of how to use dsmacc in python
"""
from dsmacc import base_model, dynenv_model, gasplusiso_model
from scipy.constants import Avogadro as N_A, R

# Choose model (cri or geoschem)
modelname = 'cri'

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
    'N2O5': 0.0028, 'HNO3': 0.754, 'PAN': 0.6496, 'PPN': 0.0825,
    'CH4': 1879., 'HCHO': 2.7141, 'CH3OH': 6.1246, 'HCOOH': 0.1881,
    'C2H4': 1.458, 'C2H6': 8.658, 'C3H8': 6.59, 'NC4H10': 7.171,
    'CH3CHO': 1.8779, 'CH3COCH3': 2.3931, 'MEK': 0.7013,
}

if modelname == 'geoschem':
    """
    requires linked FAST-JX
    """
    cri2gc = {
        'CH3CHO': 'ALD2', 'CH3COCH3': 'ACET', 'NC4H10': 'ALK4',
        'HCHO': 'CH2O', 'CH3OH': 'MOH'
    }
    initcond = {cri2gc.get(k, k): v for k, v in initcond.items()}
    initcond['O2'] = 0.2095e9
    initcond['N2'] = 0.7905e9
    del initcond['C2H4']

#*** This will only work if CRI is compiled with kpp_lsode
#*** other solvers update RCONST during integration steps,
#*** which cannot be superseded by updateenv. Test with and
#*** without RCONST override. If the answer is identical, then
#*** the solver must not be kpp_lsode
class override_na(gasplusiso_model):
    def custom_after_rconst(self):
        rxnids = self.find_eqn(prd='NA', rct='HNO3', return_index=True, return_name=False)
        pyglob = self.pyglob
        for rxnid in rxnids:
            pyglob.rconst[rxnid] *= 0

mod = override_na(
    modelname=modelname, outconcpath='boxisoconc.dat', outratepath='boxisorate.dat',
    mech2iso={'NA': 'NO3', 'HNO3': 'NO3', 'SA': 'SO4'},
    iso2mech={'NO3': 'HNO3', 'SO4': 'SA'},
    isoconst={  # currently const are provided in moles/m3, which is isorropia's unit
        'NH4': 8.462914138362976e-04,  # CRI does not have NH3 or NH4, so providing as a constant
        'NA': 8.62068966e-10,  # CRI does not have sodium, so providing as a constant
        'CL': 8.62068966e-10  # CRI does not have chloride, so providing as a constant
    },
    #verbose=0
)

mod.run(startdate, 24, 180, conc_ppb = initcond, globvar = globvar)
