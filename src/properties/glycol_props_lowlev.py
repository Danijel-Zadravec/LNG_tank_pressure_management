# -*- coding: utf-8 -*-
"""

"""
import CoolProp.CoolProp as CP
import numpy as np
import CoolProp
fluid = 'INCOMP::MEG-50%'



class GlycolClass():
    def __init__(self, mas_udio_gly):
        self.states=CoolProp.AbstractState("INCOMP", 'MEG')
        self.states.set_mass_fractions([mas_udio_gly])

        #self.states=CoolProp.AbstractState("BICUBIC&HEOS", fluid)

    def set_T(self, T):
        self.states.update(CoolProp.PT_INPUTS, 300000.0, T)
        #ph = self.get_phase()
        #assert ((ph == 5) or (ph == 2)), 'phase: ' + str(ph) + ' p: ' + str(p) + ' T: ' + str(T)


    def get_density(self):
        return self.states.rhomass()
    def get_lambda(self):
        return self.states.conductivity() #W/mK
    def get_dynamic_viscosity(self):
        return self.states.viscosity()
    def get_kinematic_viscosity(self):
        mu = self.get_dynamic_viscosity()
        ro = self.get_density()
        nu = mu/ro
        return nu
    def get_cp(self):
        return self.states.cpmass()
    def get_mass_enthalpy(self):
        return self.states.hmass()





def enthalpy(T):
    entalpija = CP.PropsSI('H', 'T', T, 'P', 100000.0, fluid) #J/kK
    return entalpija

if __name__ == "__main__":
    glik = GlycolClass(0.5)
    glik.set_T(300)
    glik.get_density()
    glik.get_lambda()
    glik.get_dynamic_viscosity()
    glik.get_kinematic_viscosity()
    glik.get_cp()
    glik.get_mass_enthalpy()


    T_1 = 350.0
    T_2 = 300.0
    glik.set_T(T_1)
    h_1 = glik.get_mass_enthalpy()
    glik.set_T(T_2)
    h_2 = glik.get_mass_enthalpy()
    dh = h_2-h_1

    T_m = 0.5*(T_1+T_2)
    glik.set_T(T_m)
    cm = glik.get_cp()
    cdT =cm*(T_2-T_1)
    razlika =cdT-dh
    rel_razl =razlika/dh*100 #%


