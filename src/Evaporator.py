# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import scipy.optimize as opt

import src.properties.glycol_props as glyc
from src.properties.liq_props_lowlev import MetLiq
from src.properties.vap_props_lowlev import MetVap, MetVap_slow

met_liq = MetLiq()
met_vap = MetVap()


class EVparams:
    def __init__(self, areas, k):
        '''
        Heat exchanger parameters

        Parameters
        ----------
        area : TYPE
            DESCRIPTION.
        k : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        '''
        self.As = areas
        self.k = k


class EV_glycol:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, T_out, mas_flow):
        self.T_in = T_in #K
        self.T_out = T_out #K  ######################
        self.mas_flow = mas_flow #kg/s #######################
        self.h_mas_in = glyc.enthalpy(T_in) #J/kg #######################
        self.h_mas_out =  glyc.enthalpy(T_out) #J/kg ###################


class EV_lng:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, p, mol_flow):
        self.T_in = T_in #K
        self.p = p #Pa
        met_liq.set_p(p)
        T_sat = met_liq.get_T() + 0.05
        self.T_sat = T_sat
        self.mol_flow = mol_flow #kmol/s
        met_liq.set_pT(p, T_in)
        self.H_molar_in = met_liq.get_mol_enthalpy() #kJ/kmol
        met_vap.set_p_sat(p)
        self.H_molar_out =  met_vap.get_mol_enthalpy() #kJ/kmol

class StatesCommonEV:
    def __init__(self):
        self.pi1=0.0 #############
        self.pi2=0.0 ##############
        self.A = 0.0 ###########
        self.heat_flow = 0.0


class EV_States:
    def __init__(self, lng, glycol, common):
        self.lng = lng
        self.glycol = glycol
        self.common = common
        self.common.heat_flow = self.calc_heat_flow()
    def calc_heat_flow(self):
        enthalpy_in = self.lng.H_molar_in
        enthalpy_out = self.lng.H_molar_out
        flow_l = self.lng.mol_flow
        heat_flow = (enthalpy_out-enthalpy_in) * flow_l*1000 #W
        return heat_flow

class SaveEvap:
    '''
    Saving heat exchanger states
    '''
    def __init__(self, states):
        '''
        https://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
        '''
        #atributi u listu
        self.lng_st_list = [a for a in dir(states.lng) if not a.startswith('__') and not callable(getattr(states.lng, a))]
        #atributi u listu
        self.glycol_st_list = [a for a in dir(states.glycol) if not a.startswith('__') and not callable(getattr(states.glycol, a))]
        #atributi zajednoicki
        self.common_st_list = [a for a in dir(states.common) if not a.startswith('__') and not callable(getattr(states.common, a))]

        self.lng = {i : [] for i in self.lng_st_list}
        self.glycol = {i : [] for i in self.glycol_st_list}
        self.common ={i : [] for i in self.common_st_list}
    def save_states(self,states):
        #iteriraj po svim atributima
        for lng_st in self.lng_st_list:
            #trenutna vrijednost atributa
            current_state = getattr(states.lng, lng_st)
            # spremi u dictionary - dodaj zadnju vrijednost liste
            self.lng[lng_st].append(current_state)
        for glycol_st in self.glycol_st_list:
            current_state = getattr(states.glycol, glycol_st)
            self.glycol[glycol_st].append(current_state)
        for common_st in self.common_st_list:
            current_state = getattr(states.common, common_st)
            self.common[common_st].append(current_state)

class Evaporator:
    '''
    klasa izmjenjivač topline
    '''
    def __init__(self, ev_params, states):
        self.params = ev_params
        self.states = states
        self.save = SaveEvap(self.states) #inicijalizacija spremanja rezultata
        self.save.save_states(self.states) #spremanje pocetnih stanja


    def update_states(self, lng_flow, lng_T_in, p, glycol_T_in):
        self.states.lng.T_in = lng_T_in
        self.states.lng.p = p
        met_liq.set_p(p)
        T_sat = met_liq.get_T() + 0.05
        self.states.lng.T_sat = T_sat
        met_liq.set_pT(p, lng_T_in)
        self.states.lng.H_molar_in = met_liq.get_mol_enthalpy() #kJ/kmol
        met_vap.set_p_sat(p)
        self.states.lng.H_molar_out =  met_vap.get_mol_enthalpy() #kJ/kmol
        self.states.lng.mol_flow = lng_flow
        self.states.common.heat_flow = self.states.calc_heat_flow() #W

        self.states.glycol.T_in = glycol_T_in
        self.states.glycol.h_mas_in = glyc.enthalpy(glycol_T_in)
        self.calc_glycol()
        self.save.save_states(self.states) #spremanje pocetnih stanja

    def calc_glycol(self):
        T_sat_lng = self.states.lng.T_sat
        T_in_gl =self.states.glycol.T_in
        h_gl_in = self.states.glycol.h_mas_in
        flow_l = self.states.lng.mol_flow
        heat_flow = self.states.common.heat_flow
        #print(heat_flow)
        k = self.params.k
        As=self.params.As
        #print(As)
        success = False
        A_res = None
        for A in As:
            T_glyc_out_init = 250.0
            (_, _, _, Fi_glyc_init) = evap_calc(
                T_glyc_out_init, T_in_gl, T_sat_lng, A, k)
            if Fi_glyc_init > heat_flow:
                sucess = False
                h_gl_out = 0.0
                glyc_flow = 0.0
                #print('a prevelika: ' + str(A))
            else:
                if isinstance(A_res, type(None)):
                    A_res = A
                opt_rez = opt.minimize_scalar(T_func, bounds=(250.0, T_in_gl),
                                              args=(T_in_gl, T_sat_lng, A, k, heat_flow), method='bounded')
                T_glyc_out = opt_rez.x
                (pi1, pi2, _, Fi_glyc) = evap_calc(
                    T_glyc_out, T_in_gl, T_sat_lng, A, k)
                if abs(Fi_glyc - heat_flow) < 0.01:
                    sucess = True
                    self.states.glycol.T_out = T_glyc_out
                    h_gl_out =  glyc.enthalpy(T_glyc_out) #J/kg
                    glyc_flow = heat_flow/(h_gl_in - h_gl_out) #kg/s
                    self.states.glycol.h_mas_out = h_gl_out
                    self.states.glycol.mas_flow = glyc_flow
                    #print('success!  ' + str(T_glyc_out))
                    self.states.common.A = A
                    break
                # print('nema rjesenja: ' + str(A))
                # print('T: ' + str(T_glyc_out))
                # print('heat flow: ' + str(Fi_glyc))
                # print(" ")
                sucess = False
                h_gl_out = 0.0
                glyc_flow = 0.0
        T_in_lng = self.states.lng.T_in
        assert sucess == True, str(A_res) + ' ' + str(flow_l) + ' ' +str(heat_flow) + ' ' +  str(T_in_lng) + ' ' + str(T_sat_lng)
        self.states.common.pi1 = pi1
        self.states.common.pi2 = pi2

def evap_calc(T_out_gl, T_in_gl, T_sat, A, k):
    pi1 =(T_in_gl-T_out_gl) /(T_in_gl-T_sat)
    pi2 = - np.log(1-pi1)
    C1 = A*k/pi2
    Fi = C1 * (T_in_gl - T_out_gl)
    return (pi1, pi2, C1, Fi)

def T_func(T_out_gl, T_in_gl, T_sat, A, k, heat_flow):
    (_, _, _, Fi) = evap_calc(T_out_gl, T_in_gl, T_sat, A, k)
    error = (heat_flow-Fi)**2
    return error

if __name__ == "__main__":
    #primjer koristenja:
    areas_lst = [5., 4., 3.5, 3., 2.5, 2.0, 1.5, 1.25, 1.0, 0.8, 0.6] #m2
    k_val = 1200. #W/m2K
    evaporator_params = EVparams(areas_lst, k_val)
    pressure = 400000. #Pa
    T_lng_in = 113.0 #K
    lng_mol_flow = 0.09 #kmol/s
    evaporator_lng =EV_lng(T_lng_in, pressure, lng_mol_flow)

    T_glycol_in = 318.0
    T_glycol_out = 318.0
    glycol_flow = 0.0

    evaporator_glycol =EV_glycol(T_glycol_in, T_glycol_out, glycol_flow)
    st_common = StatesCommonEV()
    evaporator_states = EV_States(evaporator_lng, evaporator_glycol, st_common)

    ev = Evaporator(evaporator_params, evaporator_states)
    met_liq.set_p(pressure)
    T_in_lng = met_liq.get_T()-2.0
    ev.update_states(0.015, T_in_lng, pressure, T_glycol_in) #maksimalan protok 0.09 kmol/s