"""
This file is an example of how to use dsmacc in python
"""
from dsmacc import model, dynenv
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

mod = model(modelname = modelname, outconcpath = 'boxconc.dat', outratepath = 'boxrate.dat')

mod.run(startdate, 24, 180, conc_ppb = initcond, globvar = globvar)
