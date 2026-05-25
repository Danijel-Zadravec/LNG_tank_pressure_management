# -*- coding: utf-8 -*-
"""

"""
import CoolProp.CoolProp as CP
import scipy.optimize as opt
fluid = 'INCOMP::MEG-50%'



def density(T):
    '''
    Proracun gustoce kapljevine LNG u zavisnosti o temperaturi.
    Tlak zasicenja za temperaturu.
    Coolprop - linearna regresija

    Parameters
    ----------
    T : float
        Temperatura [K]

    Returns
    -------
    ro : float
        gustoca [kg/m^3]

    '''
    ro = CP.PropsSI('D', 'T', T, 'Q', 0, fluid)
    return ro


def viscosity(T):
    return 0.000001



def cp(T, p):
    kapacitet = CP.PropsSI('C', 'T', T, 'P', 100000.0, fluid) #J/kK
    return kapacitet


def enthalpy(T):
    entalpija = CP.PropsSI('H', 'T', T, 'P', 100000.0, fluid) #J/kK
    return entalpija



density(300)
