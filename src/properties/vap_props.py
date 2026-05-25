import CoolProp.CoolProp as CP
import numpy as np
M = 16.04 #kg/kmol
fluid = 'Methane'


M = np.array(16.04) #kg/kmol - molarna masa metana
OMEGA = 0.001142 #acentric factor
P_CRIT = 4.599e6 #critical pressure, Pa
T_CRIT = 190.564 #critical temperature, K
Rm = 8.3145 #J/molK
B = 0.0778*Rm*T_CRIT/P_CRIT
f_omega = 0.37464 + 1.54226*OMEGA - 0.26992*OMEGA**2
A = 0.45724 * Rm**2 * T_CRIT**2 / P_CRIT


def molar_volume(T,p):
    '''
    Proracun molarnog volumena iz P-R EOS

    Parameters
    ----------
    T : float
        temperatura, K.
    p : float
        tlak, Pa.

    Returns
    -------
    v : float
        molarni volumen, m^3 / kmol.

    '''
    z_vap = cardano(T, p) # -
    v = z_vap * Rm * T / p # m^3/mol
    v_kmol  = v * 1000 # m^3 / kmol
    return v_kmol


def expansion_coef(T, p):
    dT = 0.1
    v2 = molar_volume(T+dT/2., p)
    v1 = molar_volume(T-dT/2, p)
    v = molar_volume(T, p)
    beta = (v2-v1)/(dT)/v
    return beta


def cardano(t, p):
    '''
    Proracun kompresibilnosti Z Cardanovom formulom iz Peng-Robinsonove EOS

    Parameters
    ----------
    t : float
        temperatura, K.
    p : float
        tlak, Pa.

    Returns
    -------
    z : float
        compressibility factor, -.

    '''
    T_r = t / T_CRIT
    a = A * (1. + f_omega * (1. - T_r**(0.5)))**2.
    c1 = B * p / (Rm * t)
    c2 = a * p / (Rm**2 * t**2)
    U = c1 - 1.
    S = c2 - 3*c1**2 - 2*c1
    T = c1**3 + c1**2 - c2*c1
    P = (3*S - U**2) / 3
    Q = 2*U**3 / 27 - U*S/3 + T
    D = (P/3)**3 + (Q/2)**2
    #assert D < 0, 'D < 0 - 3 realna rjesenja'
    if D < 0.0:
        z = car1(P, Q, U)
    else:
        z = car2(P, Q, U, D)
    return z




#@jit(nopython=True, cache=True)
def car1(P, Q, U):
    Theta = (-P**3./27.)**0.5
    Fi = np.arccos(-Q/(2*Theta))
    z1 = 2 * Theta**(1/3.) * np.cos(Fi/3.) - U/3.0
    z2 = 2 * Theta**(1/3.) * np.cos(Fi/3. + 2*np.pi/3.) - U/3.0
    z3 = 2 * Theta**(1/3.) * np.cos(Fi/3. + 4*np.pi/3.) - U/3.0
    z0 = np.maximum(z1, z2)
    z = np.maximum(z0, z3)
    return z

#@jit(nopython=True, cache=True)
def car2(P, Q, U, D):
    term1 = (D**(0.5)-Q/2.)**(1./3.)
    term2 = P/(3.*term1)
    z = term1 - term2 - U/3.
    return z





#@jit(nopython=True, cache=True)
def viscosity(T, p):
    mi = CP.PropsSI('V', 'T', T, 'P', p, fluid)
    return mi


#@jit(nopython=True, cache=True)
def kinematic_viscosity(T, p):
    mi = viscosity(T, p)
    ro = density(T, p)
    nu = mi/ro
    return nu






def thermal_diffusivity(T, p):
    lmbda = thermal_conductivity(T, p)
    ro = density(T, p)
    cp = heat_capacity(T, p)
    a = lmbda/(ro*cp)
    return a

def thermal_conductivity(T, p):
    lmbda = CP.PropsSI('L', 'T', T, 'P', p, fluid) #kJ/kmolK
    return lmbda

def heat_capacity(T, p):
    kapacitet = CP.PropsSI('C', 'T', T, 'P', p, fluid) #kJ/kmolK
    return kapacitet











def density(T, p):
    '''

    Proracun gustoce pare LNG u zavisnosti o temperaturi.
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
    if CP.PropsSI('T', 'P', p, 'Q', 1.0, fluid) < T:
        ro = CP.PropsSI('D', 'T', T, 'P', p, fluid)
    else:
        print('Vap state density!!')
        ro = CP.PropsSI('D', 'T', T, 'Q', 1.0, fluid)
    return ro


def viscosity(T, p):
    '''

    Proracun gustoce pare LNG u zavisnosti o temperaturi.
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
    if CP.PropsSI('T', 'P', p, 'Q', 1.0, fluid) < T:
        visc = CP.PropsSI('V', 'T', T, 'P', p, fluid)
    else:
        print('Vap state viscosity!!')
        visc = CP.PropsSI('V', 'T', T, 'Q', 1.0, fluid)
    return visc #Pas


def enthalpy(T, p):
    '''

    Parameters
    ----------
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    enthalpy_molar : TYPE
        DESCRIPTION.

    '''
    if CP.PropsSI('T', 'P', p, 'Q', 1.0, fluid) < T:
        h_mol = CP.PropsSI('HMOLAR', 'T', T, 'P', p, fluid) #kJ/kmol
    else:
        print('Vap state enthalpy!!')
        print(T)
        print(p)
        h_mol = CP.PropsSI('HMOLAR', 'T', T, 'Q', 1.0, fluid)  #kJ/kmol
    return h_mol


def enthalpy_saturation(p):
    '''

    Parameters
    ----------
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    enthalpy_molar : TYPE
        DESCRIPTION.

    '''
    h_mol = CP.PropsSI('HMOLAR', 'P', p, 'Q', 1.0, fluid)

    return h_mol




def sat_vap_enthalpy(T):
    h_mol = CP.PropsSI('HMOLAR', 'T', T, 'Q', 1.0, fluid)
    return h_mol


def states_ent_dens(molar_ent, molar_dens):
    '''
    TO DO!!!!!!!!!!!!!!

    Parameters
    ----------
    ent : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    temperature = CP.PropsSI('T', 'HMOLAR', molar_ent, 'DMOLAR', molar_dens, fluid)
    pressure = CP.PropsSI('P', 'HMOLAR', molar_ent, 'DMOLAR', molar_dens, fluid)
    return (temperature, pressure)



def Cmp(T, p):
    kapacitet = CP.PropsSI('CPMOLAR', 'T', T, 'P', p, fluid) #kJ/kmolK
    return kapacitet
