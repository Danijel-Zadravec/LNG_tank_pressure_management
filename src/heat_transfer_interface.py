# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import scipy.optimize as opt
#import src.properties.liq_props as liq
from src.properties.liq_props_lowlev import MetLiq
from src.properties.vap_props_lowlev import MetVap

import matplotlib.pyplot as plt

g = 9.81
met_liq = MetLiq()
met_vap = MetVap()

def interf_ht_fun(T_int, Tl, Tv, L_ekv, p, area):
    alpha_liq, _, _, _ = liq_alpha_interf(Tl, T_int, L_ekv)
    alpha_vap, _, _, _ = vap_alpha_interf(Tv, T_int, p, L_ekv)
    alpha_vap = max(alpha_vap, 0.01)
    resistance = (1./alpha_liq + 1./alpha_vap) / area
    heat_flow = (Tv - Tl) / resistance
    heat_flow_mid = (Tv - T_int) * alpha_vap * area
    error = (heat_flow-heat_flow_mid)**2
    return error


def liq_alpha_interf(Tl, T_int, L_ekv):
    T = (Tl + T_int) / 2.
    dT = Tl - T_int
    met_liq.set_T(T)
    beta = met_liq.get_beta()
    a = met_liq.get_thermal_diff()
    nu = met_liq.get_kinematic_viscosity()
    conductivity = met_liq.get_lambda()
    surface = 'down' #orijentacija ploce
    alpha, Nu, Pr, Ra = hor_heat_transfer_coef(dT, surface,  beta, L_ekv, a, nu, conductivity)
    return (alpha, Nu, Pr, Ra)


def vap_alpha_interf(Tv, T_int, p, L_ekv):
    T = (Tv + T_int) / 2.
    dT = Tv-T_int
    met_vap.set_pT_robust(p, Tv)
    beta =  met_vap.get_beta() #0.007526359799240458
    a = met_vap.get_thermal_diff() #4.925508533067638e-06 #
    nu =  met_vap.get_kinematic_viscosity() #3.925878367930164e-06 #
    conductivity =  met_vap.get_lambda() #0.014714359414435208 #
    surface = 'up' #orijentacija ploce
    alpha, Nu, Pr, Ra =hor_heat_transfer_coef(dT, surface,  beta, L_ekv, a, nu, conductivity)
    return (alpha, Nu, Pr, Ra)




def hor_heat_transfer_coef(dT, surface, beta, L_ekv, a, nu, conductivity):
    (Nu, Pr, Ra) = nusselt(dT, surface, beta, L_ekv, a, nu)
    alpha = Nu * conductivity / L_ekv
    return (alpha, Nu, Pr, Ra)

def rayleigh(dT, beta, L, a, nu):
    Ra = g * beta * np.absolute(dT) * L**3. / (a * nu)
    if Ra < 1.0e4:
        print('dt: ' + str(dT) + ' a: ' + str(a) + ' nu: ' + str(nu) + ' b: ' + str(beta) + ' Ra: ' + str(Ra))
    return Ra


#@jit(nopython=True, cache=True)
def unstable(Pr, Ra):

    if Ra > 1.0e4 and Ra < 1.0e7:
        assert Pr > 0.7, 'Prandtl < 0.7P'
        Nu = 0.54 * Ra**(1./4.)
    elif Ra  >= 1.0e7 and Ra < 50e11:
        Nu = 0.15 * Ra**(1./3.)
    else:
        Nu = 0.0
        assert False, "unstable Ra outside range"
    return Nu

#@jit(nopython=True, cache=True)
def stable(Pr, Ra):
    Nu = 0.52 * Ra**(1./5.)
    Ra = min(Ra, 30e11)
    #Ra = max(Ra, 1.0e4+1)

    #print(Ra)
    assert Ra > 1.0e4, "stable Ra outside range"
    assert Ra < 50e11, "stable Ra outside range"
    return Nu

def nusselt(dT, surface, beta, L, a, nu):
    Pr = nu/a
    #if Pr > 30.0:
    #    print('a: ' + str(a) + ' nu: ' + str(nu))
    Ra = rayleigh(dT, beta, L, a, nu)

    if surface == 'up':
        if dT > 0.0: #plate is cold, dT = T_fluid - T_wall
            Nu = stable(Pr, Ra)
        else:
            Nu = unstable(Pr, Ra)
    elif surface == 'down':
        if dT > 0.0: #cold plate
            Nu = unstable(Pr, Ra) #lower side of cold plate
        else:
            Nu = stable(Pr, Ra)
    else:
        Nu = 0.
        assert False, 'krivi surface'
    return (Nu, Pr, Ra)



if __name__ == "__main__":
    a_liq = liq_alpha_interf(112, 114, 1.0)
    a_vap = vap_alpha_interf(140, 114, 200000, 1.0)
    kk = 1/(1./a_liq + 1./a_vap)
    p=120000
    area_int = 0.4989130422727751
    L_ekv = 0.1711115281835091
    interf_ht_fun(113, 112, 130, L_ekv, p, area_int)

