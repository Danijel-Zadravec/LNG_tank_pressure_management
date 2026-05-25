# -*- coding: utf-8 -*-

import numpy as np
import scipy.optimize as opt
#import src.properties.vap_props as vap
#import src.properties.liq_props as liq
import src.properties.glycol_props as glyc
from src.properties.vap_props_lowlev import MetVap, MetVap_slow

met_vap = MetVap()



class HEparams:
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

class LngStream:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, T_out, mol_flow, p):

        self.T_in = T_in #K
        self.T_out = T_out #K
        self.p = p
        self.mol_flow = mol_flow #kmol/s
        met_vap.set_pT(p, T_in)
        self.H_molar_in = met_vap.get_mol_enthalpy()*1000.0 #J/kmol
        met_vap.set_pT(p, T_out)
        self.H_molar_out =  met_vap.get_mol_enthalpy()*1000.0 #J/kmol

class GlycolStream:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, T_out, mas_flow):

        self.T_in = T_in #K
        self.T_out = T_out #K
        self.mas_flow = mas_flow #kg/s
        self.h_mas_in =  glyc.enthalpy(T_in) #J/kg
        self.h_mas_out =  glyc.enthalpy(T_out) #J/kg


class StatesCommonHE:
    def __init__(self):
        self.pi1=0.0
        self.pi2=0.0
        self.pi3=0.0
        self.A = 0.0
        self.heat_flow = 0.0

class StatesHE:
    def __init__(self, lng, glycol, common):
        self.lng = lng
        self.glycol = glycol
        self.common = common
        self.common.heat_flow = self.calc_heat_flow()
    def calc_heat_flow(self):
        H_molar_in = self.lng.H_molar_in
        H_molar_out = self.lng.H_molar_out
        flow_l = self.lng.mol_flow
        heat_flow = (H_molar_out - H_molar_in) * flow_l
        # print(Cmp_out)
        # print(Cmp_in)
        # print(flow_l)
        # print('HF:  ' + str(heat_flow))
        # print('')
        return heat_flow

class SaveHE:
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



class HeatExchanger:
    '''
    klasa izmjenjivač topline
    '''
    def __init__(self, he_params, states):
        self.params = he_params
        self.states = states
        self.save = SaveHE(self.states) #inicijalizacija spremanja rezultata
        self.save.save_states(self.states) #spremanje pocetnih stanja


    def update_states(self, lng_flow, lng_T_in, lng_T_out, glycol_T_in, p):
        self.states.lng.T_in = lng_T_in
        self.states.lng.p = p
        met_vap.set_pT(p, lng_T_in)
        self.states.lng.H_molar_in =met_vap.get_mol_enthalpy()*1000.0
        self.states.lng.T_out = lng_T_out
        met_vap.set_pT(p, lng_T_out)
        self.states.lng.H_molar_out = met_vap.get_mol_enthalpy()*1000.0
        self.states.lng.mol_flow = lng_flow
        self.states.common.heat_flow = self.states.calc_heat_flow()

        self.states.glycol.T_in = glycol_T_in
        self.states.glycol.h_mas_in =  glyc.enthalpy(glycol_T_in)
        self.calc_glycol()
        self.save.save_states(self.states) #spremanje pocetnih stanja



    def calc_glycol(self):
        T_out_l = self.states.lng.T_out
        T_in_l = self.states.lng.T_in
        T_in_gl =self.states.glycol.T_in
        #flow_l = self.lng.mol_flow
        heat_flow = self.states.common.heat_flow
        #print('T_out: ' + str(T_out_l))
        #print('T_in: ' + str(T_in_l))
        #print('')
        #print('')

        C_lng = heat_flow / (T_out_l - T_in_l) #toplinski kapacitet glikola, kJ/K
        k = self.params.k
        As=self.params.As
        #print(As)
        pi1 = (T_in_l - T_out_l) / (T_in_l - T_in_gl) # pi 1 parametar
        # print(" ")
        # print("pi1: "+str(pi1))
        success = False
        for A in As:
            pi2 = k * A / C_lng # pi 2 parametar
            #print('pi2: ' + str(pi2))
            #print(" ")
            #print("pi2: "+str(pi2))
            opt_rez = opt.minimize_scalar(pi_3_func, bounds=(0, 1), args=(pi1, pi2), method='bounded')
            pi3 = opt_rez.x # pi 3 parametar
            if pi3 < 1.0:
                pass
            else:
                print("pi3>1!!!")
            C_glyc = C_lng / pi3 # toplinski kapacitet glikola
            T_out_gl = T_in_gl - heat_flow / C_glyc
            error = pi_3_func(pi3, pi1, pi2)
            # print('error: ' + str(error))
            # print('A: '+ str(A))
            # print("pi3: "+str(pi3))
            # print(T_out_gl)
            if error < 0.00001 and T_out_gl > 250:
                success = True
                self.states.common.A = A

                break
        assert success, 'problem pregrijač površina'


        self.states.glycol.T_out = T_out_gl
        h_glyc_out = glyc.enthalpy(T_out_gl)
        self.states.glycol.h_mas_out = h_glyc_out
        h_glyc_in =self.states.glycol.h_mas_in
        flow_glyc = heat_flow / (h_glyc_in - h_glyc_out) #kg/s
        self.states.glycol.mas_flow = flow_glyc
        self.states.common.pi1=pi1
        self.states.common.pi2=pi2
        self.states.common.pi3=pi3

def pi_3_func(pi3, pi1, pi2):
    lhs = pi1
    exponent = np.exp(-(1-pi3)*pi2)
    rhs = (1-exponent) / (1 - pi3*exponent)
    error = (rhs-lhs)**2
    return error



