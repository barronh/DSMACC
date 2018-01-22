"""
This file is an example of how to use dsmacc in python
"""
from dsmacc import model, dynenv
from scipy.constants import Avogadro
globvar = {'LAT': 30.0, 'LON': -90.0, 'TEMP': 288.0, 'PRESS': 101325.0, 'JDAY_GMT': 2006186.0, 'H2O': 9500000.0 / 1e9 * Avogadro * 1e-6 * 101325. / 298./ 8.314 }
startdate = 2006186.0
initcond = {'HNO3': 0.754, 'MEK': 0.7013, 'CH4': 1879., 'NO': 1.7879, 'NC4H10': 7.171, 'NO3': 0.0002, 'C2H6': 8.658, 'PPN': 0.0825, 'NO2': 5.9328, 'N2O5': 0.0028, 'CH3CHO': 1.8779, 'CH3OH': 6.1246, 'HCOOH': 0.1881, 'SO2': 13.63, 'CH3COCH3': 2.3931, 'C3H8': 6.59, 'O3': 41.3038, 'HCHO': 2.7141, 'CO': 142.2, 'C2H4': 1.458, 'PAN': 0.6496}

mod = model(modelname = 'cri', outconcpath = 'boxconc.dat', outratepath = 'boxrate.dat')

mod.run(startdate, 24, 180, conc_ppb = initcond, globvar = globvar)
