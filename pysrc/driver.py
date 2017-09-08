from scipy.constants import Avogadro, R, centi, nano
import numpy as np
import kpp

def h2o_from_rh_and_temp(RH, TEMP):
    """
    Return H2O in molecules/cm**3 from RH (0-100) and
    TEMP in K
    """
    TC = TEMP - 273.15
    frh = RH / 100.
    svp_millibar = 6.11 * 10**((7.5 * TC)/(TC+237.3))
    svp_pa = svp_millibar * 100
    vp_pa = svp_pa * frh
    molecule_per_cubic_m = vp_pa * Avogadro / R / TEMP
    molecule_per_cubic_cm = molecule_per_cubic_m * centi**3
    #print RH, TEMP, molecule_per_cubic_cm
    return molecule_per_cubic_cm

#Prepare kpp objects for use
pyglob = kpp.pyglob
pyrate = kpp.pyrate
pymon = kpp.pymon

# Get spc_names as a string
spc_names = [n.decode('ASCII') for n in np.char.strip(kpp.pymon.getnames().copy().view('S24')[:, 0]).tolist()]

# Prepare integrator
integrate =kpp.pyint.integrate
RSTATE = np.zeros(30, dtype = 'd')
ERROR = np.zeros(1, dtype = 'd')
ICNTRL_U = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], dtype = 'i')
pyglob.atol[:] = 1e-3
pyglob.rtol[:] = 1e-5

# Set day in YYYYJJJ format
base_jday = 2014182
# Set time for start
pyglob.time = 8.5 * 3600
dt=180

# Integrate for 8 hours
tend = 8*3600 + pyglob.time

def updateenv(LON_degE = 0., LAT_degN = 30., TEMP_K = 298.15, PRESS_Pa = 101325, O2_vmr = 0.21, H2_vmr = 1e-6, RH_pct = 20):
    """
    Default environment updater
    - sets LON, LAT, TEMP, PRESS, M, O2, N2, H2, H2O in molecules/cm3
    - sets CFACTOR to convert from ppb to molecules/cm3
    - updates rate constants
    
    Can be edited called with arguments
    """
    t = pyglob.time
    fday = pyglob.time / 24. / 3600
    pyglob.JDAY_GMT = base_jday + fday
    # Set longitude to GMT and latitude to midlatitude
    pyglob.LON=LON_degE
    pyglob.LAT=LAT_degN

    # Set basic temperature and pressure
    pyglob.TEMP = TEMP_K
    pyglob.PRESS = PRESS_Pa

    # Set some globals
    pyglob.M=pyglob.PRESS * Avogadro / R / pyglob.TEMP * centi **3
    pyglob.CFACTOR = pyglob.M * nano#{ppb-to-molecules/cm3}
    pyglob.O2=O2_vmr*pyglob.M
    pyglob.N2=pyglob.M-pyglob.O2
    pyglob.H2=H2_vmr*pyglob.M
    pyglob.H2O=h2o_from_rh_and_temp(RH_pct, pyglob.TEMP)
    # Calculate the CFACTOR (c_air [=] molecules/cm3
    pyrate.update_rconst()

# Use default environment
updateenv()

def initialize(default = 1e-32, **conc_ppb):
    """
    set concentrations to initial values
    """
    CFACTOR = pyglob.CFACTOR
    
    # Create default concentration and set for all species
    pyglob.c[:] = default*CFACTOR;

    # Set initial values for ozone (O3), isoprene (C5H8) and nitric oxide (NO)
    for k, v in conc_ppb.items():
        pyglob.c[spc_names.index(k)] = v * CFACTOR

initialize(O3 = 30, C5H8 = 850., NO = 128)

# Globals to write out.
# Write out initial results
out = '\t'.join(['time', 'cfactor'] + spc_names)
out += '\n%d\t%.8e\t' % (pyglob.time, pyglob.CFACTOR) + '\t'.join(['%.8e' % (ci/pyglob.CFACTOR) for ci in pyglob.c[:]])

print('JDAY_LST', pyglob.THETA)
# Loop through time until at end
while pyglob.time < tend:
    tout = pyglob.time+dt
    while pyglob.time < tout:
        updateenv()
        istatus, rstatus, ierr = integrate( tin = pyglob.time, tout = tout, icntrl_u = ICNTRL_U)
        if ierr != 1:
            raise ValueError('Integration failed at ' + str(t))
        pyglob.time = rstatus[0]
    # show Local time for clarity
    print(pyglob.JDAY_GMT+pyglob.LON/15./24, pyglob.THETA)
    # Write out tout results
    out += '\n%d\t%.8e\t' % (pyglob.time, pyglob.CFACTOR) + '\t'.join(['%.8e' % (ci/pyglob.CFACTOR) for ci in pyglob.c[:]])

    
# Archive results in a file
outfile = open('cri_isop.tsv', 'w')
print(out, file = outfile)
