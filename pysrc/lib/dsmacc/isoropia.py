import numpy as np

# 345678901234567890123456789012345678901234567890123456789012345678901234567890
g_NA = np.array([0, 0, 0])
g_SO4 = np.array([0, 0, 0])
g_NH4 = np.array([1, 0, 0])
g_NO3 = np.array([0, 1, 0])
g_CL = np.array([0, 0, 1])
g_CA = np.array([0, 0, 0])
g_K = np.array([0, 0, 0])
g_MG = np.array([0, 0, 0])

al_NA = np.array([0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
al_SO4 = np.array([0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
al_NH4 = np.array([0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])
al_NO3 = np.array([0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0])
al_CL = np.array([0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
al_CA = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0])
al_K = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0])
al_MG = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1])

as_NA = np.array([1, 0, 1, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
as_SO4 = np.array([0, 0, 0, 0, 1, 1, 1, 1, 2, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0])
as_NH4 = np.array([0, 1, 0, 1, 0, 2, 0, 1, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
as_NO3 = np.array([1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2, 0])
as_CL = np.array([0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 2])
as_CA = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
as_K = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 1, 0, 0, 0])
as_MG = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1])

w_NA = 0
w_SO4 = 1
w_NH4 = 2
w_NO3 = 3
w_CL = 4
w_CA = 5
w_K = 6
w_MG = 7

wSPCS = 'NA SO4 NH4 NO3 CL CA K MG'.split()


def isoropia(WI, RHI, TEMPI, REVERSE=False, METASTABLE=False):
    """
    Parameters
    ----------
    WI        : double array (moles/m3) containing total sodium, sulfate,
                ammonium, nitrate, chloride, calcium, potassium, magnesium
    RHI       : relative humidity (0, 1)
    TEMPI     : temperature K
    REVERSE   : default False, run in reverse mode
    METASTABLE: default False, liquide only

    Returns
    -------
    out       : dictionary WT, GAS, AERLIQ, AERSLD, SCASI, OTHER

    Notes
    -----

    Descriptions of options copied from isoropiaIIcode.F
    (1-based indices in fortran should be decremented)

      INPUT:
      1. [WI]
         DOUBLE PRECISION array of length [8].
         Concentrations, expressed in moles/m3. Depending on the type of
         problem solved (specified in CNTRL(1)), WI contains either
         GAS+AEROSOL or AEROSOL only concentratios.
         WI(1) - sodium
         WI(2) - sulfate
         WI(3) - ammonium
         WI(4) - nitrate
         WI(5) - chloride
         WI(6) - calcium
         WI(7) - potassium
         WI(8) - magnesium

      2. [RHI]
         DOUBLE PRECISION variable.
         Ambient relative humidity expressed on a (0,1) scale.

      3. [TEMPI]
         DOUBLE PRECISION variable.
         Ambient temperature expressed in Kelvins.

      4. [CNTRL]
         DOUBLE PRECISION array of length [2].
         Parameters that control the type of problem solved.

         CNTRL(1): Defines the type of problem solved.
         0 - Forward problem is solved. In this case, array WI contains
             GAS and AEROSOL concentrations together.
         1 - Reverse problem is solved. In this case, array WI contains
             AEROSOL concentrations only.

         CNTRL(2): Defines the state of the aerosol
         0 - The aerosol can have both solid+liquid phases (deliquescent)
         1 - The aerosol is in only liquid state (metastable aerosol)

      OUTPUT:
      1. [WT]
         DOUBLE PRECISION array of length [8].
         Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3.
         If the foreward probelm is solved (CNTRL(1)=0), array WT is
         identical to array WI.
         WT(1) - total sodium
         WT(2) - total sulfate
         WT(3) - total ammonium
         WT(4) - total nitrate
         WT(5) - total chloride
         WT(6) - total calcium
         WT(7) - total potassium
         WT(8) - total magnesium

      2. [GAS]
         DOUBLE PRECISION array of length [03].
         Gaseous species concentrations, expressed in moles/m3.
         GAS(1) - NH3
         GAS(2) - HNO3
         GAS(3) - HCl

      3. [AERLIQ]
         DOUBLE PRECISION array of length [15].
         Liquid aerosol species concentrations, expressed in moles/m3.
         AERLIQ(01) - H+(aq)
         AERLIQ(02) - Na+(aq)
         AERLIQ(03) - NH4+(aq)
         AERLIQ(04) - Cl-(aq)
         AERLIQ(05) - SO4--(aq)
         AERLIQ(06) - HSO4-(aq)
         AERLIQ(07) - NO3-(aq)
         AERLIQ(08) - H2O
         AERLIQ(09) - NH3(aq) (undissociated)
         AERLIQ(10) - HNCl(aq) (undissociated)
         AERLIQ(11) - HNO3(aq) (undissociated)
         AERLIQ(12) - OH-(aq)
         AERLIQ(13) - Ca2+(aq)
         AERLIQ(14) - K+(aq)
         AERLIQ(15) - Mg2+(aq)

      4. [AERSLD]
         DOUBLE PRECISION array of length [19].
         Solid aerosol species concentrations, expressed in moles/m3.
         AERSLD(01) - NaNO3(s)
         AERSLD(02) - NH4NO3(s)
         AERSLD(03) - NaCl(s)
         AERSLD(04) - NH4Cl(s)
         AERSLD(05) - Na2SO4(s)
         AERSLD(06) - (NH4)2SO4(s)
         AERSLD(07) - NaHSO4(s)
         AERSLD(08) - NH4HSO4(s)
         AERSLD(09) - (NH4)4H(SO4)2(s)
         AERSLD(10) - CaSO4(s)
         AERSLD(11) - Ca(NO3)2(s)
         AERSLD(12) - CaCl2(s)
         AERSLD(13) - K2SO4(s)
         AERSLD(14) - KHSO4(s)
         AERSLD(15) - KNO3(s)
         AERSLD(16) - KCl(s)
         AERSLD(17) - MgSO4(s)
         AERSLD(18) - Mg(NO3)2(s)
         AERSLD(19) - MgCl2(s)

      5. [SCASI]
         CHARACTER*15 variable.
         Returns the subcase which the input corresponds to.

      6. [OTHER]
         DOUBLE PRECISION array of length [9].
         Returns solution information.

         OTHER(1): Shows if aerosol water exists.
         0 - Aerosol is WET
         1 - Aerosol is DRY

         OTHER(2): Aerosol Sulfate ratio, defined as (in moles/m3) :
                   (total ammonia + total Na) / (total sulfate)

         OTHER(3): Sulfate ratio based on aerosol properties that defines
                   a sulfate poor system:
                   (aerosol ammonia + aerosol Na) / (aerosol sulfate)

         OTHER(4): Aerosol sodium ratio, defined as (in moles/m3) :
                   (total Na) / (total sulfate)

         OTHER(5): Ionic strength of the aqueous aerosol (if it exists).

         OTHER(6): Total number of calls to the activity coefficient
                   calculation subroutine.

         OTHER(7): Sulfate ratio w/ crustal species, defined as (in moles/m3) :
                   (total ammonia + total crustal species + total Na) /
                   (total sulfate)

         OTHER(8): Crustal species + sodium ratio, defined as (in moles/m3) :
                   (total crustal species + total Na) / (total sulfate)

         OTHER(9): Crustal species ratio, defined as (in moles/m3) :
                   (total crustal species) / (total sulfate)

     *** COPYRIGHT 1996-2006, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY,
     *** GEORGIA INSTITUTE OF TECHNOLOGY
     *** WRITTEN BY ATHANASIOS NENES
     *** UPDATED BY CHRISTOS FOUNTOUKIS

    """
    from dsmacc import _isoropia
    WI = np.asarray(WI, dtype='d')
    RHI = np.asarray(RHI, dtype='d')
    TEMPI = np.asarray(TEMPI, dtype='d')
    CNTRL = np.zeros(2, dtype='d')
    if REVERSE:
        CNTRL[0] = 1
    if METASTABLE:
        CNTRL[1] = 1
    out = _isoropia.dsmacc_isoropia(WI, RHI, TEMPI, CNTRL)
    outd = dict(
        WT=out[0], GAS=out[1], AERLIQ=out[2],
        AERSLD=out[3], SCASI=out[4], OTHER=out[5]
    )
    return outd


def total_gasaero(isoout):
    """
    Parameters
    ----------
    isoout : output dictionary from isoropia

    Returns
    -------
    gasaer : dictionary with GAS and AERO items.
             GAS and AERO are in moles/m3 and use
             indices w_*
    """
    WT = isoout['WT']
    GAS = isoout['GAS']
    AERLIQ = isoout['AERLIQ']
    AERSLD = isoout['AERSLD']
    TA = np.zeros_like(WT)
    TA[w_NA] = (AERLIQ * al_NA).sum() + (AERSLD * as_NA).sum()
    TA[w_SO4] = (AERLIQ * al_SO4).sum() + (AERSLD * as_SO4).sum()
    TA[w_NH4] = (AERLIQ * al_NH4).sum() + (AERSLD * as_NH4).sum()
    TA[w_NO3] = (AERLIQ * al_NO3).sum() + (AERSLD * as_NO3).sum()
    TA[w_CL] = (AERLIQ * al_CL).sum() + (AERSLD * as_CL).sum()
    TA[w_CA] = (AERLIQ * al_CA).sum() + (AERSLD * as_CA).sum()
    TA[w_K] = (AERLIQ * al_K).sum() + (AERSLD * as_K).sum()
    TA[w_MG] = (AERLIQ * al_MG).sum() + (AERSLD * as_MG).sum()
    GA = np.zeros_like(WT)
    GA[w_NA] = (GAS * g_NA).sum()
    GA[w_SO4] = (GAS * g_SO4).sum()
    GA[w_NH4] = (GAS * g_NH4).sum()
    GA[w_NO3] = (GAS * g_NO3).sum()
    GA[w_CL] = (GAS * g_CL).sum()
    GA[w_CA] = (GAS * g_CA).sum()
    GA[w_K] = (GAS * g_K).sum()
    GA[w_MG] = (GAS * g_MG).sum()
    return dict(GAS=GA, AERO=TA)


def test_conserve(isoin, isoout):
    import numpy as np
    WT = isoout['WT']
    # AERLIQ = isoout['AERLIQ']
    # AERSLD = isoout['AERSLD']
    out = total_gasaero(isoout)
    GAS = out['GAS']
    AERO = out['AERO']
    RC = GAS + AERO
    print('WI == WT?', np.array_equal(isoin, WT))
    print('WI == Reconstructed?')
    print(' - SPCN, Equal, WT, RC, RC/WT-1')
    for si, spcn in enumerate(wSPCS):
        print(
            ' -', spcn, np.isclose(WT[si], RC[si]),
            WT[si], RC[si], RC[si] / WI[si] - 1
        )


if __name__ == '__main__':
    #     WI taken from ...
    #     from PseudoNetCDF import pncopen
    #     f = pncopen(
    #         'CAMx.v6.40.midwest.36.12.noMPI.20020603.avrg.grd02',
    #         format='uamiv'
    #     ).sliceDimensions(TSTEP=18, LAY=0, ROW=45, COL=45)
    #     nacl = 0.05 # CAMx chmdat.f
    #     WI[w_NA] = nacl / 1e6 / 58. # sodium
    #     WI[w_SO4] = f.variables['PSO4'] / 1e6 / 96. # sulfate
    #     WI[w_NH4] = f.variables['PNH4'] / 1e6 / 18. # ammonium
    #     WI[w_NH4] += f.variables['NH3'] * PRESS / R / TEMPI
    #     WI[w_NO3] = f.variables['PNO3'] / 1e6 / 62.  # nitrate
    #     WI[w_NO3] += f.variables['HNO3'] * PRESS / R / TEMPI
    #     WI[w_CL] = nacl / 1e6 / 58. # chloride
    #     WI[w_CA] = 0. # calcium
    #     WI[w_K] = 0. # potassium
    #     WI[w_MG] = 0. # magnesium
    WI = np.array([
        8.62068966e-10, 1.41172123e-08, 1.23921670e-02,
        1.77862160e-02, 8.62068966e-10, 0., 0., 0.
    ], dtype='d')  # input moles/m3
    RHI = np.array(.5, dtype='d')  # (0, 1)
    TEMPI = np.array(295., dtype='d')  # temperature K
    PRESS = np.array(101325., dtype='d')  # pressure Pa
    # dim(2); forward = 0, backward = 1; solid/liquid = 0, solidonly = 1
    CNTRL = np.array([0, 0], dtype='d')
    print('TEST with CAMx = 1e-30')
    out = isoropia(WI, .5, 295.)
    aout = total_gasaero(out)
    print(wSPCS)
    print('INP', WI)
    print('GAS', aout['GAS'])
    print('AER', aout['AERO'])
    test_conserve(WI, out)
    WI = np.maximum(WI, 1e-30)
    print('TEST with minvals = 1e-30')
    out = isoropia(WI, .5, 295.)
    test_conserve(WI, out)
    WI = np.maximum(WI, 1e-20)
    print('TEST with minvals = 1e-20')
    out = isoropia(WI, .5, 295.)
    test_conserve(WI, out)
