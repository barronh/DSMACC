from __future__ import print_function
from warnings import warn
from scipy.constants import Avogadro, R, centi, nano
import numpy as np
from . import kpp
#Prepare kpp objects for easy access
pyglob = kpp.pyglob
pyrate = kpp.pyrate
pymon = kpp.pymon

# Get spc_names as a string
spc_names = [n.decode('ASCII') for n in np.char.strip(kpp.pymon.getnames().copy().view('S24')[:, 0]).tolist()]


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


class model(object):
    def __init__(self, outpath = None, delimiter = ',', globalkeys = ['time', 'LAT', 'LON', 'PRESS', 'TEMP', 'THETA', 'H2O', 'CFACTOR']):
        """
            outpath - string for path of output
            globalkeys - keys to output with concentrations
        """
        self.outpath = outpath
        self.delimiter = delimiter
        self.globalkeys = globalkeys
        
    def custom_updateenv(self):
        """
        Called by updateenv before update_rconst. By default,
        this does nothing. Overwrite in subclass with custom 
        updates to pyglob.
        """
        return
    
    def updateenv(self, O2_vmr = 0.21, H2_vmr = 1e-6, CH4_vmr = 1.85e-6, **kwds):
        """
        Arguments:
            O3_vmr - default molecular oxygen volume mixing ratio
            H2_vmr - deafault molecular hydrogen volume mixing ratio
            CH4_vmr - default methane volume mixing ratio
        Actions:
        - uses BASE_JDAY_GMT to set JDAY_GMT
        - sets LON (degE), LAT (degN), TEMP (K), PRESS (Pa)
        - sets O2, N2, H2, H2O in molecules/cm3
        - sets CFACTOR to convert from ppb to molecules/cm3
        - calls custom_updateenv, which should be overwritten
        - updates rate constants
    
        Can be edited called with arguments
        """
        t = pyglob.time
        fday = pyglob.time / 24. / 3600
        pyglob.JDAY_GMT = pyglob.BASE_JDAY_GMT + fday

        # Calculate the CFACTOR (c_air [=] molecules/cm3
        for k, v in kwds.items():
            if hasattr(pyglob, k):
                setattr(pyglob, k, v)
            elif 'RH_pct' == k:
                pyglob.H2O=h2o_from_rh_and_temp(v, pyglob.TEMP)
            else:
                warn('%s not set; not a global variable' % (k,))

        # Set some globals
        pyglob.M=pyglob.PRESS * Avogadro / R / pyglob.TEMP * centi **3
        pyglob.CFACTOR = pyglob.M * nano#{ppb-to-molecules/cm3}
        pyglob.O2=O2_vmr*pyglob.M
        pyglob.N2=pyglob.M-pyglob.O2
        pyglob.H2=H2_vmr*pyglob.M
        pyglob.CH4=CH4_vmr*pyglob.M

        
        self.custom_updateenv()
        pyrate.update_rconst()

    def initialize(self, JDAY_GMT, conc_ppb = {}, globvar = {}, default = 1e-32):
        """
        Arguments:
            JDAY_GMT - YYYYJJJ.FFFFFF where FFFFFF is a fraction of a day
            conc_ppb - dictionary-like of initial concentrations in ppb for any species
            globvar - dictionary-like of global variables
            default - default concentration for any species not specified in conc_ppb
        
        Actions:
            - set BASE_JDAY_GMT and initial time from JDAY_GMT
            - update the global environment
            - set default concentrations
            - set specified concentrations.       
        """
        # Set the base JDAY to the integer portion of the input
        pyglob.BASE_JDAY_GMT = (JDAY_GMT // 1)
        # Set the initial time to the fraction portion of the 
        # input jday converted to seconds
        pyglob.time = (JDAY_GMT % 1) * 24 * 3600.
        
        # Use any global environment variables provided to
        # update the environment
        self.updateenv(**globvar)
    
        # Set default concentration for all species
        CFACTOR =pyglob.CFACTOR
        print(CFACTOR)
        pyglob.c[:] = default*CFACTOR;

        # Set initial values for any species in conc_ppb
        for k, v in conc_ppb.items():
            pyglob.c[spc_names.index(k)] = v * CFACTOR

    def output(self, globalkeys = None, restart = False):
        """
        Arguments:
            globalkeys - list of keys to print from global environment
            restart - boolean indicating to restart the output
        
        Actions:
            saves current state to self.out
        """
        if globalkeys is None:
            globalkeys = self.globalkeys
        vals = tuple([getattr(pyglob, gk) for gk in globalkeys] + [(ci/pyglob.CFACTOR) for ci in pyglob.c[:]])
        if restart:
            self.out = self.delimiter.join(globalkeys + spc_names)
            self.fmt = self.delimiter.join(['%d'] + ['%.8e'] * (len(vals) - 1))

        # Write out tout results
        out = self.fmt % vals
        self.out += '\n' + out

    def save(self, path = None, clear = True):
        """
        Arguments:
            path - string path for output to be saved to
            clear - boolean to erase self.out after savign
            
        Actions
        """
        if path is None:
            path = self.outpath
        if path is None:
            path = 'output.dat'
        # Archive results in a file
        outfile = open(path, 'w')
        outfile.write(self.out)
        outfile.close()
        if clear: self.out = ""

    def run(self, jday_gmt, run_hours, dt, conc_ppb, globvar, initkwds = {}, atol = 1e-3, rtol = 1e-5):
        """
        - jday_gmt, conc_ppb, globvar and initkwds are used to initialize the model (See initialize)
        - run_hours and dt are used to configure integration and steps
        - atol and rtol are used to configure the solver
        
        basically:
            - configure integrate inputs
                - RSTATE (30 double zeros)
                - ERROR = (1 double zero)
                - ICNTRL_U = (20 integer zeros)
            - set global variables atol and rtol that integrator uses
            - call initialize 
            - output initial values
            
                
            
        """
        self.dt = dt
        # Prepare integrator
        integrate =kpp.pyint.integrate
        RSTATE = np.zeros(30, dtype = 'd')
        ERROR = np.zeros(1, dtype = 'd')
        ICNTRL_U = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], dtype = 'i')
        pyglob.atol[:] = atol
        pyglob.rtol[:] = rtol

        # Globals to write out.
        # Write out initial results
        self.initialize(jday_gmt, conc_ppb = conc_ppb, globvar = globvar, **initkwds)
        tend = pyglob.time + run_hours * 3600
        updateenv = self.updateenv
        output = self.output
        output(restart = True)
        
        # Loop through time until at end
        while pyglob.time < tend:
            tout = pyglob.time+dt
            while pyglob.time < tout:
                updateenv()
                istatus, rstatus, ierr = integrate( tin = pyglob.time, tout = tout, icntrl_u = ICNTRL_U)
                if ierr != 1:
                    raise ValueError('Integration failed at ' + str(pyglob.time))
                pyglob.time = rstatus[0]
            # show Local time for clarity
            #print(pyglob.JDAY_GMT+pyglob.LON/15./24, pyglob.THETA)
            # Write out tout results
            output(restart = False)

    
        self.save()

class dynenv(model):
    def __init__(self, outpath = None, delimiter = ',', envdata = None, emissdata = None, bkgdata = None, globalkeys = ['time', 'LAT', 'LON', 'PRESS', 'TEMP', 'THETA', 'H2O', 'CFACTOR', 'PBL']):
        """
        Arguments:
            outpath - path for output to be saved
            envdata - dictionary of vectors (must contain JDAY_GMT) and optionally other
                      global variables that influence model. PBL can also be supplied in meters.
            emissdata - dictionary of vectors (must contain JDAY_GMT) and optionally species
                        data in moles/area/s where area is cm2
            bkgdata - dictionary of scalars optionally provides time-constant species
                        data in ppb
            globalkeys - keys for outputing global variables
        """
        import pandas as pd
        self.globalkeys = globalkeys
        self.outpath = outpath
        self.delimiter = delimiter
        self.envdata = envdata
        self.emissdata = emissdata
        self.updatebkg = not bkgdata is None
        if self.updatebkg:
            bkgc_ppb = self._bkgc_ppb = np.zeros_like(pyglob.c[:])
            for k, v in bkgdata.items():
                ki = spc_names.index(k)
                bkgc_ppb[ki] = v
                
    def updateenv(self, **kwds):
        """
        updateenv is overwritten to first calculate interpolated
        values and then call the original model updateenv
        
        uses optional keyword PBL to entrain air
        """
        globvar = kwds.copy()
        t = pyglob.BASE_JDAY_GMT + pyglob.time / 3600. / 24.
        xt = self.envdata['JDAY_GMT']
        for k, yv in self.envdata.items():
            if k != 'JDAY_GMT':
                v = np.interp(t, xt, yv)
                globvar[k] = v
        emisvar = {}
        for k, yv in self.emissdata.items():
            if k != 'JDAY_GMT':
                v = np.interp(t, xt, yv)
                emisvar[k] = v
        
        if 'PBL' in globvar:
            pyglob.PBL = newpbl = globvar.pop('PBL')
            if self.updatebkg:
                if hasattr(self, 'oldpbl'):
                    dpbl = newpbl - self.oldpbl
                    if dpbl > 0:
                        fnew = dpbl / self.oldpbl
                        fold = 1 - fnew
                        pyglob.c[:] = pyglob.c[:] * fold + self._bkgc_ppb[:] * pyglob.CFACTOR * fnew
                
            newpbl_cm = newpbl * 100
            for ek, ev in emisvar.items():
                ki = spc_names.index(ek)
                pyglob.c[:] += self.dt * ev * Avogadro / newpbl_cm
            self.oldpbl = newpbl
        
        model.updateenv(self, **globvar)
