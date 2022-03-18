__all__ = ['base_model', 'dynenv_model', 'gasplusiso_model']

from collections import defaultdict
from warnings import warn
from scipy.constants import Avogadro, R, centi, nano
import numpy as np

_minval = 1e-36

# 345678901234567890123456789012345678901234567890123456789012345678901234567890
# You can only have 1 instance at a time.
_instance_count = defaultdict(lambda: 0)


class base_model(object):
    def __del__(self):
        global _instance_count
        _instance_count[self.modelname] -= 1

    def __init__(
        self, outconcpath=None, outratepath=None, delimiter=',',
        globalkeys=(
            'time', 'JDAY_GMT', 'LAT', 'LON', 'PRESS', 'TEMP', 'THETA',
            'CFACTOR', 'H2O', 'N2', 'O2', 'M', 'H2', 'CH4'
        ), modelname='cri', verbose=0
    ):
        """
        Arguments
        ---------
        outconcpath: str
            path for saved output concentrations (ppb)
        outratepath: str
            path for saved output rates (1/s, cm3/molecules/s)
        delimeter: str
            Output delimiter
        globalkeys: iterable
            keys to output with concentrations
        modelname: str
            Name of model in DSMACC.
        """
        global _instance_count
        self.verbose = verbose
        self.modelname = modelname
        if _instance_count[self.modelname] > 0:
            raise ValueError(
                'You may only have one DSMACC instance open per session;'
                + 'currently %d' % _instance_count[self.modelname]
            )
        else:
            _instance_count[self.modelname] += 1
            _instance_count

        from importlib import import_module
        # Prepare kpp objects for easy access
        kpp = self.kpp = import_module(
            '.'+self.modelname, package='dsmacc'
        )
        self.pyint = kpp.pyint
        self.pyglob = kpp.pyglob
        self.pyrate = kpp.pyrate
        self.pymon = kpp.pymon
        # Get spc_names as a string
        self.spc_names = [
            n.decode('ASCII')
            for n in np.char.strip(
                self.pymon.getnames().copy().view('S24')[:, 0]
            ).tolist()
        ]

        # Get spc_names as a string
        self.eqn_names = [
            n.decode('ASCII')
            for n in np.char.strip(
                self.pymon.geteqnnames().copy().view('S100')[:, 0]
            ).tolist()
        ]

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

    def updateenv(
        self, O2_vmr=0.21, H2_vmr=550e-9, CH4_vmr=1.85e-6, **kwds
    ):
        """
        Arguments:
            O2_vmr - default molecular oxygen volume mixing ratio
            H2_vmr - deafault molecular hydrogen volume mixing ratio
            CH4_vmr - default methane volume mixing ratio
        Actions:
        - uses BASE_JDAY_GMT to set JDAY_GMT
        - sets LON (degE), LAT (degN), TEMP (K), PRESS (Pa)
        - sets O2, N2, H2, H2O, and CH4 in molecules/cm3 (from C if available)
        - sets CFACTOR to convert from ppb to molecules/cm3
        - calls custom_before_rconst, which should be overwritten
        - updates rate constants
        - calls custom_after_rconst, which should be overwritten

        Can be edited called with arguments
        """
        from . import envutil
        pyglob = self.pyglob
        t = pyglob.time
        spc_names = self.spc_names
        fday = t / 24. / 3600
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
                pyglob.H2O = envutil.h2o_from_rh_and_temp(v, pyglob.TEMP)
            else:
                warn('%s not set; not a global variable' % (k,))

        # Set some globals
        pyglob.M = pyglob.PRESS * Avogadro / R / pyglob.TEMP * centi**3
        pyglob.CFACTOR = pyglob.M * nano  # {ppb-to-molecules/cm3}

        def setfix(pyglob, key, defval, kwds):
            if key in spc_names:
                idx = spc_names.index(key)
                newval = pyglob.c[idx]
            elif key in kwds:
                newval = kwds[key]
            else:
                newval = defval
            setattr(pyglob, key, newval)

        setfix(pyglob, 'O2', O2_vmr * pyglob.M, kwds)
        setfix(pyglob, 'H2', H2_vmr * pyglob.M, kwds)
        setfix(pyglob, 'CH4', CH4_vmr * pyglob.M, kwds)
        setfix(
            pyglob, 'N2', pyglob.M - pyglob.O2 - pyglob.H2 - pyglob.CH4, kwds
        )

        self.custom_before_rconst()
        self.pyrate.update_rconst()
        self.custom_after_rconst()

    def initialize(self, JDAY_GMT, conc_ppb=None, globvar=None, default=1e-32):
        """
        Arguments
        ---------
        JDAY_GMT: float
            YYYYJJJ.FFFFFF where FFFFFF is a fraction of a day
        conc_ppb: dict
            dictionary-like of initial concentrations in ppb for any species
        globvar: dict
            dictionary-like of global variables
        default: dict
            default concentration for any species not specified in conc_ppb

        Notes
        -----
        Actions:
        - set BASE_JDAY_GMT and initial time from JDAY_GMT
        - update the global environment
        - set default concentrations
        - set specified concentrations.
        """
        if conc_ppb is None:
            conc_ppb = {}
        if globvar is None:
            globvar = {}
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
        CFACTOR = pyglob.CFACTOR
        pyglob.c[:] = default * CFACTOR

        # Set initial values for any species in conc_ppb
        for k, v in conc_ppb.items():
            if k in spc_names:
                pyglob.c[spc_names.index(k)] = v * CFACTOR
            else:
                warn('%s not in mechanism' % k)

    def output(self, globalkeys=None, restart=False):
        """
        Arguments
        ---------
        globalkeys: list
            list of keys to print from global environment
        restart: bool
            boolean indicating to restart the output

        Actions:
            saves current state to self.out
        """
        spc_names = self.spc_names
        eqn_names = self.eqn_names
        pyglob = self.pyglob
        if globalkeys is None:
            globalkeys = self.globalkeys
        outconcvals = tuple(
            [getattr(pyglob, gk, np.nan) for gk in globalkeys]
            + [(ci/pyglob.CFACTOR) for ci in pyglob.c[:]]
        )
        outratevals = tuple(
            [getattr(pyglob, gk, np.nan) for gk in globalkeys]
            + [(ri) for ri in pyglob.rconst[:]]
        )
        if restart:
            conckeys = list(globalkeys) + list(spc_names)
            ratekeys = list(globalkeys) + list(eqn_names)
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

    def save(self, concpath=None, ratepath=None, clear=True):
        """
        Arguments
        ---------
        concpath: str
            string path for output concentrations to be saved to
        ratepath: str
            string path for output rates to be saved to
        clear: bool
            If true, erase self.out after saving

        Actions
        -------
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
        if clear:
            self.outconc = ""
        if clear:
            self.outrate = ""

    def run(
        self, jday_gmt, run_hours, dt, conc_ppb, globvar, initkwds=None,
        atol=1e-3, rtol=1e-5
    ):
        """
        Overview
        --------
        basically:
        - configure integrate inputs
            - RSTATE (30 double zeros)
            - ERROR = (1 double zero)
            - ICNTRL_U = (20 integer zeros)
        - set global variables atol and rtol that integrator uses
        - call initialize
        - output initial values

        Arguments
        ---------
        jday_gmt: float
            YYYYJJJ.SSSSS
        conc_ppb: dict
            concentrations in ppb
        globvar: dict
            Global variables
        initkwds: dict
            key words used to initialize the model (See initialize)
        run_hours: scalar
            duration in hours for the integration
        dt: scalar
            integrater time steps in seconds
        atol, rtol: arrays
            absolute and relative tolerance scalars or arrays
        """
        if initkwds is None:
            initkwds = {}
        pyglob = self.pyglob
        self.dt = dt
        # Prepare integrator
        integrate = self.pyint.integrate
        RSTATE = np.zeros(30, dtype='d')
        ERROR = np.zeros(1, dtype='d')
        ICNTRL_U = np.array([0] * 20, dtype='i')
        pyglob.atol[:] = atol
        pyglob.rtol[:] = rtol

        # Globals to write out.
        # Write out initial results
        self.initialize(
            jday_gmt, conc_ppb=conc_ppb, globvar=globvar, **initkwds
        )
        tend = pyglob.time + run_hours * 3600
        updateenv = self.updateenv
        output = self.output
        output(restart=True)

        # Loop through time until at end
        ierr = 1
        while pyglob.time < tend and ierr == 1:
            tout = pyglob.time+dt
            while pyglob.time < tout and ierr == 1:
                updateenv()
                istatus, rstatus, ierr = integrate(
                    tin=pyglob.time, tout=tout, icntrl_u=ICNTRL_U
                )
                if ierr != 1:
                    self.save()
                    raise ValueError(
                        'Integration failed at ' + str(pyglob.time)
                        + '; saved partial run'
                    )
                pyglob.time = rstatus[0]
            # show Local time for clarity
            # print(pyglob.JDAY_GMT+pyglob.LON/15./24, pyglob.THETA)
            # Write out tout results
            output(restart=False)

        self.save()

    def find_eqn(self, spc=None, rct=None, prd=None, return_index=False, return_name=True):
        import re
        eqns = [(eqi, eq) for eqi, eq in enumerate(self.eqn_names)]
        if spc is not None:
            tmpre = re.compile(fr'\b{spc}\b')
            eqns = [(eqi, eq) for eqi, eq in eqns if tmpre.search(eq) is not None]
        if prd is not None:
            tmpre = re.compile(fr'\b{prd}\b')
            eqns = [(eqi, eq) for eqi, eq in eqns if tmpre.search(eq.split('-->')[1]) is not None]
        if rct is not None:
            tmpre = re.compile(fr'\b{rct}\b')
            eqns = [(eqi, eq) for eqi, eq in eqns if tmpre.search(eq.split('-->')[0]) is not None]
        if return_index and return_name:
            return eqns
        if return_index:
            return [eqi for eqi, eq in eqns]
        
        return [eq for eqi, eq in eqns]

    def find_spc(self, spc):
        if spc in self.spc_names:
            return self.spc_names.index(spc)
        else:
            return None

    def get_conc(self, spc):
        spci = self.find_spc(spc)
        if spci is None:
            return None

        return self.pyglob.c[spci]

    def aero_dict(self, outunits=True):
        from . import isoropia
        showvals = self.aerosol_phase
        if outunits:
            showvals = showvals / self.pyglob.CFACTOR
        return dict(zip(isoropia.wSPCS, showvals))

class dynenv_model(base_model):
    def __init__(
        self, outconcpath=None, outratepath=None, delimiter=',',
        envdata=None, emissdata=None, bkgdata=None,
        globalkeys=(
            'time', 'JDAY_GMT', 'LAT', 'LON', 'PRESS', 'TEMP',
            'THETA', 'H2O', 'CFACTOR', 'PBL'
        ), modelname='cri', verbose=0
    ):
        """
        The dynenv_model allows for emissions and global variables (key in
        globalkeys) to evolve over time.

        Arguments
        ---------
        outconcpath: str
            path for output concentrations to be saved (ppb)
        outratepath: str
            path for output rates to be saved (1/s, cm3/molecules/s)
        envdata: dict
            dictionary of vectors (must contain JDAY_GMT) and optionally other
            global variables that  influence model. PBL can also be supplied
            in meters.
        emissdata: dict
            dictionary of vectors (must contain JDAY_GMT) and optionally
            species data in moles/area/s where area is cm2
        bkgdata: dict
            dictionary of scalars optionally provides time-constant species
            data in ppb
        globalkeys: iterable
            names of global variables to output
        modelname: str
        """
        # import pandas as pd
        super(dynenv_model, self).__init__(
            outconcpath=outconcpath, outratepath=outratepath,
            delimiter=delimiter, globalkeys=globalkeys,
            modelname=modelname, verbose=verbose
        )
        spc_names = self.spc_names
        self.envdata = envdata
        self.emissdata = emissdata
        self.updatebkg = bkgdata is not None
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
        if self.envdata is not None:
            xt = self.envdata['JDAY_GMT']
            for k, yv in self.envdata.items():
                if k != 'JDAY_GMT':
                    v = np.interp(t, xt, yv)
                    globvar[k] = v
        emisvar = {}
        if self.emissdata is not None:
            xt = self.emissdata['JDAY_GMT']
            for k, yv in self.emissdata.items():
                if k != 'JDAY_GMT':
                    v = np.interp(t, xt, yv)
                    emisvar[k] = v

        if 'PBL' in globvar:
            pyglob.PBL = newpbl = np.array(globvar['PBL'])
            if self.updatebkg:
                if hasattr(self, 'oldpbl'):
                    dpbl = newpbl - self.oldpbl
                    if dpbl > 0:
                        fnew = dpbl / self.oldpbl
                        fold = 1 - fnew
                        pyglob.c[:] = (
                            pyglob.c[:] * fold + self._bkgc_ppb[:]
                            * pyglob.CFACTOR * fnew
                        )

            newpbl_cm = newpbl * 100
            for ek, ev in emisvar.items():
                if ek in spc_names:
                    ki = spc_names.index(ek)
                    pyglob.c[ki] += self.dt * ev * Avogadro / newpbl_cm
                else:
                    warn(ek + ' in emissions, but not mechanism')
            self.oldpbl = newpbl

        super(dynenv_model, self).updateenv(**globvar)


class gasplusiso_model(dynenv_model):
    def __init__(self, *args, mech2iso=None, iso2mech=None, isoconst=None, **kwds):
        """
        See dynenv_model. This is a thin wrapper that adds a place holder for
        call_isorropia.

        mech2iso: dict
            mechspc, isospc key/value pairs
        iso2mech: dict
            isospc, mechspc key/value pairs
        """
        dynenv_model.__init__(self, *args, **kwds)
        self._mech2iso = None
        self._iso2mech = None
        self._isoconst = None
        self.mech2iso = mech2iso
        self.iso2mech = iso2mech
        self.isoconst = isoconst
        self.aerosol_phase = np.zeros(8, dtype='f')

    @property
    def mech2iso(self):
        return self._mech2iso

    @mech2iso.setter
    def mech2iso(self, mech2isomap):
        from . import isoropia
        ispci = isoropia.wSPCS.index
        spc_names = self.spc_names
        mspci = spc_names.index
        self._mech2iso = []
        if mech2isomap is None:
            mech2isomap = {}
            if 'SO4' in spc_names:
                mech2isomap['SO4'] = 'SO4'
            if 'NH3' in spc_names:
                mech2isomap['NH3'] = 'NH4'
            if 'NH4' in spc_names:
                mech2isomap['NH3'] = 'NH4'
            if 'HNO3' in spc_names:
                mech2isomap['HNO3'] = 'NO3'
            if 'NA' in spc_names:
                mech2isomap['NA'] = 'NO3'
            if 'CL' in spc_names:
                mech2isomap['CL'] = 'CL'
            if 'CA' in spc_names:
                mech2isomap['CA'] = 'CA'
            if 'K' in spc_names:
                mech2isomap['K'] = 'K'
            if 'MG' in spc_names:
                mech2isomap['MG'] = 'MG'
            print('Default mech2iso')
            print(mech2isomap)

        for srckey, destkey in mech2isomap.items():
            self._mech2iso.append((ispci(destkey), mspci(srckey)))

    @property
    def iso2mech(self):
        return self._iso2mech

    @iso2mech.setter
    def iso2mech(self, iso2mechmap):
        from . import isoropia
        ispci = isoropia.wSPCS.index
        spc_names = self.spc_names
        mspci = spc_names.index
        self._iso2mech = []
        if iso2mechmap is None:
            iso2mechmap = {}
            if 'SA' in spc_names:
                iso2mechmap['SO4'] = 'SA'
            if 'NH3' in spc_names:
                iso2mechmap['NH4'] = 'NH3'
            if 'HNO3' in spc_names:
                iso2mechmap['NO3'] = 'HNO3'
            if 'CL' in spc_names:
                iso2mechmap['CL'] = 'CL'
            if 'CA' in spc_names:
                iso2mechmap['CA'] = 'CA'
            if 'K' in spc_names:
                iso2mechmap['K'] = 'K'
            if 'MG' in spc_names:
                iso2mechmap['MG'] = 'MG'
            print('Default iso2mech')
            print(iso2mechmap)

        for srckey, destkey in iso2mechmap.items():
            self._iso2mech.append((mspci(destkey), ispci(srckey)))

    @property
    def isoconst(self):
        return self._isoconst

    @isoconst.setter
    def isoconst(self, isoconstmap):
        from . import isoropia
        ispci = isoropia.wSPCS.index
        mspci = self.spc_names.index
        self._isoconst = []
        if isoconstmap is None:
            isoconstmap = {}
            if 'NA' not in self.mech2iso:
                isoconstmap['NA'] = 8.62068966e-10
            if 'CL' not in self.mech2iso:
                isoconstmap['CL'] = 8.62068966e-10
            print('Default isoconst')
            print(isoconstmap)

        for destkey, val in isoconstmap.items():
            self._isoconst.append((ispci(destkey), val))

    def gas_to_iso(self):
        """
        Returns
        -------
        out : np.array
            output correctly ordred for ISOROPIA
        """
        from . import isoropia
        from scipy.constants import Avogadro
        i = isoropia
        pyglob = self.pyglob
        WI = np.zeros(8, dtype='f')
        spc_names = self.spc_names
        spci = self.spc_names.index
        # i.w_NA is sodium
        # WI[i.w_NA] = pyglob.c[pyglob.ind_NA] * Avogadro
        # The gas phase mechanism must implement the rest of these species
        for desti, srci in self.mech2iso:
            # Convert from molecules/cm3 to moles/m3
            WI[desti] = pyglob.c[srci] / Avogadro * 1e6
            pyglob.c[srci] = 0

        # previous time step aerosols
        WI[:] += (self.aerosol_phase / Avogadro * 1e6)

        for desti, val in self._isoconst:
            WI[desti] = val

        return np.ma.maximum(WI, _minval).filled(_minval)

    def process_isoout(self, iout):
        """
        Parameters
        ----------
        iout : dictionary of arrays
            keys are TOT, GAS, AERLIQ, AERSLD and correspond to isoropia
            outputs

        Returns
        -------
        None
        """
        from . import isoropia
        i = isoropia
        pyglob = self.pyglob
        # separate gas and total aerosol phase
        ga = isoropia.total_gasaero(iout)
        gases = ga['GAS']
        spc_names = self.spc_names
        spci = spc_names.index
        # map isorropia resultant gases into the gas vector
        # pyglob.c[pyglob.ind_NA] = gases[i.w_NA] * Avogadro
        for dsti, srci in self.iso2mech:
            # convert moles/m3 to molecules/cm3
            pyglob.c[dsti] = gases[srci] * Avogadro / 1e6

        # store aerosol phase sum separately
        # convert moles/m3 to molecules/cm3
        self.aerosol_phase[:] = np.ma.maximum(ga['AERO'] * Avogadro / 1e6, _minval).filled(_minval)

    def isorropia(self):
        """
        get indices and stuff you'll need create input condition arrays
        """
        from . import isoropia
        RHI = 0.5  # just a place holder.
        TEMPI = float(self.pyglob.TEMP)
        WI = self.gas_to_iso()
        if self.verbose > 1:
            print('ISO IN')
            print(list(zip(isoropia.wSPCS, WI)))
        # run isorropia
        isoresult = isoropia.isoropia(WI, RHI, TEMPI, METASTABLE=True)
        if self.verbose > 1:
            print('ISO OUT')
            print(list(zip(isoropia.wSPCS, isoresult['WT'])))
        self.process_isoout(isoresult)

    def updateenv(self, **kwds):
        """
        Thin wrapper that calls isoropia before updateenv
        """
        self.isorropia()
        return super(gasplusiso_model, self).updateenv(**kwds)
