from __future__ import print_function
from collections import defaultdict

__all__ = ['model', 'dynenv', 'envutil']

from warnings import warn
from scipy.constants import Avogadro, R, centi, nano
import numpy as np
from . import envutil


# You can only have 1 instance at a time.
_instance_count = defaultdict(lambda: 0)

class model(object):
    def __del__(self):
        global _instance_count
        _instance_count[self.modelname] -= 1
    
    def __init__(self, outconcpath = None, outratepath = None, delimiter = ',', globalkeys = ['time', 'JDAY_GMT', 'LAT', 'LON', 'PRESS', 'TEMP', 'THETA', 'H2O', 'CFACTOR'], modelname = 'cri'):
        """
            outconcpath - path for output concentrations to be saved (ppb)
            outratepath - path for output rates to be saved (1/s, cm3/molecules/s)
            globalkeys - keys to output with concentrations
        """
        global _instance_count
        self.modelname = modelname
        if _instance_count[self.modelname] > 0:
            raise ValueError('You may only have one DSMACC instance open per session; currently %d' % _instance_count[self.modelname])
        else:
            _instance_count[self.modelname] += 1
            _instance_count

        from importlib import import_module
        #Prepare kpp objects for easy access
        kpp = self.kpp = import_module('.'+self.modelname, package = 'dsmacc')
        self.pyint = kpp.pyint
        self.pyglob = kpp.pyglob
        self.pyrate = kpp.pyrate
        self.pymon = kpp.pymon
        # Get spc_names as a string
        self.spc_names = [n.decode('ASCII') for n in np.char.strip(self.pymon.getnames().copy().view('S24')[:, 0]).tolist()]

        # Get spc_names as a string
        self.eqn_names = [n.decode('ASCII') for n in np.char.strip(self.pymon.geteqnnames().copy().view('S100')[:, 0]).tolist()]


        self.outconcpath = outconcpath
        self.outratepath = outratepath
        self.delimiter = delimiter
        self.globalkeys = globalkeys
        
    def custom_before_rconst(self):
        """
        Called by updateenv before update_rconst. By default,
        this does nothing. Overwrite in subclass with custom 
        updates to pyglob.
        """
        return
    
    def custom_after_rconst(self):
        """
        Called by updateenv after update_rconst. By default,
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
        - calls custom_before_rconst, which should be overwritten
        - updates rate constants
        - calls custom_after_rconst, which should be overwritten
    
        Can be edited called with arguments
        """
        pyglob = self.pyglob
        t = pyglob.time
        fday = pyglob.time / 24. / 3600
        pyglob.JDAY_GMT = pyglob.BASE_JDAY_GMT + fday

        # Calculate the CFACTOR (c_air [=] molecules/cm3
        for k, v in kwds.items():
            if hasattr(pyglob, k):
                globv = getattr(pyglob, k)
                if globv.size > 1:
                    globv[:] = np.ones_like(globv[:]) * v
                else:
                    setattr(pyglob, k, v)
            elif 'RH_pct' == k:
                pyglob.H2O=envutil.h2o_from_rh_and_temp(v, pyglob.TEMP)
            else:
                warn('%s not set; not a global variable' % (k,))

        # Set some globals
        pyglob.M=pyglob.PRESS * Avogadro / R / pyglob.TEMP * centi **3
        pyglob.CFACTOR = pyglob.M * nano#{ppb-to-molecules/cm3}
        pyglob.O2=O2_vmr*pyglob.M
        pyglob.N2=pyglob.M-pyglob.O2
        pyglob.H2=H2_vmr*pyglob.M
        pyglob.CH4=CH4_vmr*pyglob.M

        
        self.custom_before_rconst()
        self.pyrate.update_rconst()
        self.custom_after_rconst()

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
        pyglob = self.pyglob
        spc_names = self.spc_names
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
        pyglob.c[:] = default*CFACTOR;

        # Set initial values for any species in conc_ppb
        for k, v in conc_ppb.items():
            if k in spc_names:
                pyglob.c[spc_names.index(k)] = v * CFACTOR
            else:
                warn('%s not in mechanism' % k)

    def output(self, globalkeys = None, restart = False):
        """
        Arguments:
            globalkeys - list of keys to print from global environment
            restart - boolean indicating to restart the output
        
        Actions:
            saves current state to self.out
        """
        spc_names = self.spc_names
        eqn_names = self.eqn_names
        pyglob = self.pyglob
        if globalkeys is None:
            globalkeys = self.globalkeys
        outconcvals = tuple([getattr(pyglob, gk, np.nan) for gk in globalkeys] + [(ci/pyglob.CFACTOR) for ci in pyglob.c[:]])
        outratevals = tuple([getattr(pyglob, gk, np.nan) for gk in globalkeys] + [(ri) for ri in pyglob.rconst[:]])
        if restart:
            conckeys = globalkeys + spc_names
            ratekeys = globalkeys + eqn_names
            self.outconc = self.delimiter.join(conckeys)
            self.outrate = self.delimiter.join(ratekeys)
            concfmts = ['%.8e'] * len(conckeys)
            ratefmts = ['%.8e'] * len(ratekeys)
            if 'JDAY_GMT' in conckeys:
                concfmts[conckeys.index('JDAY_GMT')] = '%.5f'
            if 'JDAY_GMT' in ratekeys:
                ratefmts[ratekeys.index('JDAY_GMT')] = '%.5f'

            self.fmtoutconc = self.delimiter.join(concfmts)
            self.fmtoutrate = self.delimiter.join(ratefmts)

        # Write out tout results
        outconc = self.fmtoutconc % outconcvals
        self.outconc += '\n' + outconc
        outrate = self.fmtoutrate % outratevals
        self.outrate += '\n' + outrate

    def save(self, concpath = None, ratepath = None, clear = True):
        """
        Arguments:
            concpath - string path for output concentrations to be saved to
            ratepath - string path for output rates to be saved to
            clear - boolean to erase self.out after savign
            
        Actions
        """
        if concpath is None:
            concpath = self.outconcpath
        if concpath is None:
            concpath = 'output.dat'
        if ratepath is None:
            ratepath = self.outratepath
        if ratepath is None:
            ratepath = 'rate.dat'
        # Archive results in a file
        outfile = open(concpath, 'w')
        outfile.write(self.outconc)
        outfile.close()
        outfile = open(ratepath, 'w')
        outfile.write(self.outrate)
        outfile.close()
        if clear: self.outconc = ""
        if clear: self.outrate = ""

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
        pyglob = self.pyglob
        self.dt = dt
        # Prepare integrator
        integrate =self.pyint.integrate
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
        ierr = 1
        while pyglob.time < tend and ierr == 1:
            tout = pyglob.time+dt
            while pyglob.time < tout and ierr == 1:
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
    def __init__(self, outconcpath = None, outratepath = None, delimiter = ',', envdata = None, emissdata = None, bkgdata = None, globalkeys = ['time', 'JDAY_GMT', 'LAT', 'LON', 'PRESS', 'TEMP', 'THETA', 'H2O', 'CFACTOR', 'PBL'], modelname = 'cri'):
        """
        Arguments:
            outconcpath - path for output concentrations to be saved (ppb)
            outratepath - path for output rates to be saved (1/s, cm3/molecules/s)
            envdata - dictionary of vectors (must contain JDAY_GMT) and optionally other
                      global variables that influence model. PBL can also be supplied in meters.
            emissdata - dictionary of vectors (must contain JDAY_GMT) and optionally species
                        data in moles/area/s where area is cm2
            bkgdata - dictionary of scalars optionally provides time-constant species
                        data in ppb
            globalkeys - keys for outputing global variables
        """
        import pandas as pd
        super(dynenv, self).__init__(outconcpath = outconcpath, outratepath = outratepath, delimiter = delimiter, globalkeys = globalkeys, modelname = modelname)
        spc_names = self.spc_names
        self.envdata = envdata
        self.emissdata = emissdata
        self.updatebkg = not bkgdata is None
        if self.updatebkg:
            bkgc_ppb = self._bkgc_ppb = np.zeros_like(self.pyglob.c[:])
            for k, v in bkgdata.items():
                ki = spc_names.index(k)
                bkgc_ppb[ki] = v
                
    def updateenv(self, **kwds):
        """
        updateenv is overwritten to first calculate interpolated
        values and then call the original model updateenv
        
        uses optional keyword PBL to entrain air
        """
        pyglob = self.pyglob
        globvar = kwds.copy()
        spc_names = self.spc_names
        pyglob = self.pyglob
        t = pyglob.BASE_JDAY_GMT + pyglob.time / 3600. / 24.
        if not self.envdata is None:
            xt = self.envdata['JDAY_GMT']
            for k, yv in self.envdata.items():
                if k != 'JDAY_GMT':
                    v = np.interp(t, xt, yv)
                    globvar[k] = v
        emisvar = {}
        if not self.emissdata is None:
            xt = self.emissdata['JDAY_GMT']
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
                if ek in spc_names:
                    ki = spc_names.index(ek)
                    pyglob.c[:] += self.dt * ev * Avogadro / newpbl_cm
                else:
                    warn(ek + ' in emissions, but not mechanism')
            self.oldpbl = newpbl
        
        super(dynenv, self).updateenv(**globvar)
