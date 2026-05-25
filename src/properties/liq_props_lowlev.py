import numpy as np
import CoolProp
import CoolProp.CoolProp as CP
import scipy.optimize as opt
M = 16.04
fluid = "Methane"
class MetLiq():
    def __init__(self):
        self.states=CoolProp.AbstractState("HEOS", fluid)
        #self.states=CoolProp.AbstractState("BICUBIC&HEOS", fluid)

    def set_T(self, T):
        self.states.update(CoolProp.QT_INPUTS, 0.0, T)
        #K

    def set_p(self,p):
        self.states.update(CoolProp.PQ_INPUTS, p, 0.0)
        #Pa

    def set_molH(self, h_molar):
        rez = opt.minimize_scalar(self._ent_fun_opt, bounds=(91., 189.0), args=(h_molar), method='bounded')
        T=rez.x
        self.set_T(T)
        #kJ/kmol
    def set_pT(self, p, T):
        self.states.update(CoolProp.PT_INPUTS, p, T)
        ph = self.get_phase()
        assert ph == 0, 'phase: ' + str(ph) + ' p: ' + str(p) + ' T: ' + str(T)
    def set_pT_all(self, p, T):
        self.states.update(CoolProp.PT_INPUTS, p, T)

    def set_pT_robust(self, p, T):
        self.states.update(CoolProp.QT_INPUTS, 0.0, T)
        T_sat = self.get_T()
        if T < T_sat:
            self.states.update(CoolProp.PT_INPUTS, p, T)




    def set_pSmol(self, p, S_molar):
        self.states.update(CoolProp.PSmolar_INPUTS, p, S_molar)
        ph = self.get_phase()
        assert ph == 0, 'phase: ' + str(ph) + ' p: ' + str(p) + ' S: ' + str(S_molar)

    def set_pHmol(self, p, H_molar):
        self.states.update(CoolProp.HmolarP_INPUTS, H_molar, p)
        ph = self.get_phase()
        assert ph == 0, 'phase: ' + str(ph) + ' p: ' + str(p) + ' H: ' + str(H_molar)
    def set_pHmol_wet(self, p, H_molar):
        self.states.update(CoolProp.HmolarP_INPUTS, H_molar, p)


    def _ent_fun_opt(self, T, h_molar):
        self.set_T(T)
        entalpija_dobivena = self.get_mol_enthalpy()
        error = (h_molar - entalpija_dobivena) ** 2
        return error
    def get_phase(self):
        #5 - gas
        #2 -  supercritical gas
        #0 - kapljevina
        #6 - 2 faze
        return self.states.phase()

    def get_lambda(self):
        return self.states.conductivity() #W/mK
    def get_density(self):
        return self.states.rhomass() #kg/m^3
    def get_mol_enthalpy(self):
        return self.states.hmolar() #J/mol ili kJ/kmol
    def get_mas_enthalpy(self):
        return self.states.hmass()
    def get_mol_entropy(self):
        return self.states.smolar() #J/molK ili kJ/kmolK

    def get_T(self):
        return self.states.T() #J/mol ili kJ/kmol
    def get_p(self):
        return self.states.p()
    def get_beta(self):
        return isobaric_exp_coef(self.get_T())
    def get_cp(self):
        return self.states.cpmass()
    def get_thermal_diff(self):
        lmb = self.get_lambda()
        ro = self.get_density()
        cp = self.get_cp()
        a = lmb/(ro*cp)
        return a
    def get_dynamic_viscosity(self):
        return self.states.viscosity()
    def get_kinematic_viscosity(self):
        mu = self.get_dynamic_viscosity()
        ro = self.get_density()
        nu = mu/ro
        return nu
    def get_Pr(self):
        mu = self.get_dynamic_viscosity()
        cp = self.get_cp()
        l = self.get_lambda()
        Pr =mu*cp/l
        return Pr



#vrela kapljevina!
C_DENSITY = np.array([381.8810787978164, 4.909851329148779, -0.07737476436109958,
                      0.00043813256292211077, -9.841388792746806e-07])
C_CP = np.array([-5848., 231.4, -1.965, 0.005753])
C_MI = np.array([0.001331, -2.231e-5, 1.332e-7, -2.743e-10])
C_LAMBDA = np.array([0.2109, 1.546e-3, -2.186e-5, 5.268e-8])



def density_old(T):
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
    ro = CP.PropsSI('D', 'T', T, 'Q', 0.0, fluid)
    return ro

def density(T):
    '''
    Proracun gustoce kapljevine metana u zavisnosti o temperaturi.
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
    ro = C_DENSITY[0] + C_DENSITY[1]*T + C_DENSITY[2]*T**2 + C_DENSITY[3]*T**3 + C_DENSITY[4]*T**4
    return ro


#@jit(nopython=True, cache=True)
def heat_capacity(T):
    topl_kapacitet = C_CP[0] + C_CP[1]*T + C_CP[2]*T**2 + C_CP[3]*T**3
    return topl_kapacitet

#@jit(nopython=True, cache=True)
def viscosity(T):
    mi = C_MI[0] + C_MI[1]*T + C_MI[2]*T**2 + C_MI[3]*T**3
    return mi

#@jit(nopython=True, cache=True)
def thermal_conductivity(T):
    lmbd = C_LAMBDA[0] + C_LAMBDA[1]*T + C_LAMBDA[2]*T**2 + C_LAMBDA[3]*T**3
    return lmbd


#@jit(nopython=True, cache=True)
def isobaric_exp_coef(T):
    dro_dt =C_DENSITY[1] +  2.*C_DENSITY[2]*T + 3.*C_DENSITY[3]*T**2 + 4.*C_DENSITY[4]*T**3
    ro = density(T)
    beta = -dro_dt / ro
    return beta



def isobaric_exp_coef_old(T):
    beta = CP.PropsSI('ISOBARIC_EXPANSION_COEFFICIENT', 'T', T, 'Q', 0.0, fluid) #Pas
    return beta #1/K

def thermal_diffusivity(T):
    lmbda =  thermal_conductivity(T)
    ro =  density(T)
    cp =  CP.PropsSI('C', 'T', T, 'Q', 0.0, fluid)
    a = lmbda/(ro*cp)
    return a
def kinematic_viscosity(T):
    mi = viscosity(T)
    ro = density(T)
    nu = mi/ro
    return nu

def thermal_conductivity_old(T):
    lmbda =  CP.PropsSI('CONDUCTIVITY', 'T', T, 'Q', 0.0, fluid)
    return lmbda

def viscosity_old(T):
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
    visc = CP.PropsSI('V', 'T', T, 'Q', 0.0, fluid) #Pas
    return visc


def enthalpy_sat(T):
    '''

    ----------
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    enthalpy_molar : TYPE
        kJ/kmol.

    '''
    enthalpy_molar = CP.PropsSI('HMOLAR', 'T', T, 'Q', 0.0, fluid)
    return enthalpy_molar


def enthalpy(T, p):
    '''

    ----------
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    enthalpy_molar : TYPE
        kJ/kmol.

    '''
    Tsat = CP.PropsSI('T', 'P', p, 'Q', 0.0, fluid)
    #assert T < Tsat, 'temperatura previsoka!'
    if T < (Tsat-0.005):
        enthalpy_molar = CP.PropsSI('HMOLAR', 'T', T, 'P', p, fluid)
    else:
        enthalpy_molar = CP.PropsSI('HMOLAR', 'T', T, 'Q', 0.0, fluid)

    return enthalpy_molar

def tempearture_sat_H(H):
    temperature = CP.PropsSI('T', 'HMOLAR', H, 'Q', 0.0, fluid) #kJ/kmol
    return temperature


def entropy(T, p):
    '''

    ----------
    T : TYPE
        DESCRIPTION.

    Returns
    -------
    enthalpy_molar : TYPE
        kJ/kmol.

    '''
    Tsat = CP.PropsSI('T', 'P', p, 'Q', 0.0, fluid)
    #assert T < Tsat, 'temperatura previsoka!'
    if T < (Tsat-0.005):
        entropy_molar = CP.PropsSI('SMOLAR', 'T', T, 'P', p, fluid) #kJ/kmolK
    else:
        entropy_molar = CP.PropsSI('SMOLAR', 'T', T, 'Q', 0.0, fluid) #kJ/kmolK
    return entropy_molar


def isentropic_H(p, S):
    enthalpy = CP.PropsSI('HMOLAR', 'SMOLAR', S, 'P', p, fluid) #kJ/kmol
    return enthalpy


def tempearture_Hp(H, p):
    temperature = CP.PropsSI('T', 'HMOLAR', H, 'P', p, fluid) #kJ/kmol
    return temperature

def temp_ent_fun_opt(T, ent):
    entalpija_dobivena = CP.PropsSI('HMOLAR', 'T', T, 'Q', 0.0, fluid)
    error = (ent - entalpija_dobivena) ** 2
    return error
def temp_ent(ent):
    '''
    # TO DO!!!!!!!!!!!!!!!
    #optim, bolje napraviti neku interpolaciju
    Parameters

    Parameters
    ----------
    ent : TYPE
        DESCRIPTION.

    Returns
    -------
    temperature, K

    '''
    #temperature = CP.PropsSI('T', 'HMOLAR', ent, 'Q', 0.0, fluid)
    rez = opt.minimize_scalar(temp_ent_fun_opt, bounds=(91., 189.0), args=(ent), method='bounded')
    temperature = rez.x
    return temperature


def saturation_temperature(p):
    '''
    Temperatura zasićenja ovisno o tlaku
    '''
    Tsat = CP.PropsSI('T', 'P', p, 'Q', 0.0, fluid)
    return Tsat

def Cmp(T, p):
    kapacitet = CP.PropsSI('CPMOLAR', 'T', T, 'Q', 0.0, fluid) #kJ/kmolK
    return kapacitet
