__all__ = ['h2o_from_rh_and_temp']


def h2o_from_rh_and_temp(RH, TEMP):
    """
    Return H2O in molecules/cm**3 from RH (0-100) and
    TEMP in K
    """
    from scipy.constants import Avogadro, R, centi
    TC = TEMP - 273.15
    frh = RH / 100.
    svp_millibar = 6.11 * 10**((7.5 * TC)/(TC+237.3))
    svp_pa = svp_millibar * 100
    vp_pa = svp_pa * frh
    molecule_per_cubic_m = vp_pa * Avogadro / R / TEMP
    molecule_per_cubic_cm = molecule_per_cubic_m * centi**3
    # print RH, TEMP, molecule_per_cubic_cm
    return molecule_per_cubic_cm
