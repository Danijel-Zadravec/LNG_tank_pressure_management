# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import scipy.optimize as opt
import time

import src.properties.glycol_props as glyc
from src.properties.liq_props_lowlev import MetLiq
from src.properties.vap_props_lowlev import MetVap

M=16.04
met_liq = MetLiq()
met_vap = MetVap()


class PipeParams:
    def __init__(self, d, length, roughness, medium, k, T_env):
        self.d = d
        self.length = length
        self.roughness = roughness
        self.medium = medium #'lng_liq, lng_vap, glycol
        self.k = k
        self.T_env = T_env
        self.A_cross = d ** 2.0 * 3.1415 / 4.0 #površina poprečnog presjeka
        self.A_wall = d*3.1415 * length

class PipeStates:
    def __init__(self, T_in, T_out, p_in, p_out, flow):
        self.T_in = T_in
        self.T_out = T_out
        self.p_in = p_in
        self.p_out = p_out
        self.flow = flow #kmol/s
        self.h_in = 0.0
        self.h_out = 0.0
        self.Fi = 0.0
        self.dT = 0.0
        self.M = 0.0
        self.ro = 0.0
        self.mu = 0.0
        self.dp = 0
        self.speed = 0
        self.friction_coeff = 0.0


class PipeSave:
    '''
    Saving pipe states
    '''
    def __init__(self, states):

        #atributi u listu
        self.state_list = [a for a in dir(states) if not a.startswith('__') and not callable(getattr(states, a))]

        self.states = {i : [] for i in self.state_list}

    def save_states(self,states):
        #iteriraj po svim atributima
        for st in self.state_list:
            #trenutna vrijednost atributa
            current_state = getattr(states, st)
            # spremi u dictionary - dodaj zadnju vrijednost liste
            self.states[st].append(current_state)

class Pipe:
    def __init__(self, params, states):
        self.params = params
        self.states = states
        medium = self.params.medium
        if medium == 'lng_liq' or medium == 'lng_vap':
            self.states.M =M
        elif medium == 'glycol':
            self.states.M = glyc.M
        else:
            self.states.M = -1
            print("krivi medij!!")
        self.save = PipeSave(self.states)
        self.save.save_states(self.states)

    def calc_properties(self, T, p):
        medium = self.params.medium
        if medium == 'lng_liq':
            met_liq.set_T(T)
            self.states.ro = met_liq.get_density()
            self.states.mu = met_liq.get_dynamic_viscosity()
        elif medium == 'lng_vap':
            met_vap.set_pT(p, T)
            self.states.ro = met_vap.get_density()
            self.states.mu = met_vap.get_dynamic_viscosity()
        elif medium == 'glycol':
            self.states.ro = glyc.density(T)
            self.states.mu = glyc.viscosity(T)
        else:
            self.states.ro = -1
            self.states.mu = -1
            print("PROBLEM calc_properties")

    def calc_speed(self):
        '''
        Proračun srednje brzine strujanja

        Returns
        -------
        None.

        '''
        M = self.states.M #kG/kmol
        ro=self.states.ro #kg/m^3
        mol_flow = self.states.flow #kmol/s
        mass_flow = mol_flow * M #kg/s
        vol_flow =mass_flow / ro
        speed =vol_flow / self.params.A_cross
        self.states.speed = speed



    def pressure_out(self):
        T = 0.5 * (self.states.T_in + self.states.T_out)
        p = self.states.p_in
        self.calc_properties(T,p)
        d = self.params.d
        length = self.params.length
        roughness = self.params.roughness
        ro = self.states.ro
        mu = self.states.mu
        self.calc_speed()
        speed =self.states.speed

        friction_coeff_rez = opt.minimize_scalar(colebrook, bounds=(0.0000001, 1.0), args=(roughness, d, speed, ro, mu), method='bounded')
        friction_coeff = friction_coeff_rez.x
        self.states.friction_coeff = friction_coeff
        #Darcy weisbach
        dp = friction_coeff * (length * ro / 2) * speed** 2 / d
        self.states.dp = dp
        self.states.p_out = self.states.p_in - dp

    def pressure_in(self):
        T = 0.5 * (self.states.T_in + self.states.T_out)
        p = self.states.p_out
        self.calc_properties(T, p)
        d = self.params.d
        length = self.params.length
        roughness = self.params.roughness
        ro=self.states.ro #kg/m^3
        mu = self.states.mu
        self.calc_speed()
        speed =self.states.speed
        friction_coeff_rez = opt.minimize_scalar(colebrook, bounds=(0.0000001, 1.0), args=(roughness, d, speed, ro, mu), method='bounded')
        friction_coeff = friction_coeff_rez.x
        self.states.friction_coeff = friction_coeff
        #Darcy weisbach
        dp = friction_coeff * (length * ro / 2) * speed** 2 / d
        self.states.dp = dp
        self.states.p_in = self.states.p_out + dp

    def temperature_out(self, p):
        k = self.params.k
        medium = self.params.medium
        A_wall = self.params.A_wall
        flow = self.states.flow
        T_in = self.states.T_in
        h_in = calc_h(T_in, p, medium)
        T_env = self.params.T_env

        T_out_opt = opt.minimize_scalar(T_out_fun, bounds=(T_in, T_in+3.0),
                                    args=(T_in, T_env, p, h_in, k, A_wall, flow, medium), method='bounded')
        T_out = T_out_opt.x
        h_out = calc_h(T_out, p, medium)
        Fi = flow * (h_out - h_in)
        self.states.h_in = h_in
        self.states.h_out = h_out
        self.states.T_out = T_out
        self.states.Fi = Fi
        self.states.dT = T_out-T_in

    def temperature_in(self, p):
        k = self.params.k
        medium = self.params.medium
        A_wall = self.params.A_wall
        flow = self.states.flow
        T_out = self.states.T_out
        T_env = self.params.T_env
        h_out = calc_h(T_out, p, medium)
        T_in_opt = opt.minimize_scalar(T_in_fun, bounds=(T_out-5.0, T_out+5.0),
                                    args=(T_out, T_env, p, h_out, k, A_wall, flow, medium), method='bounded')
        T_in = T_in_opt.x
        h_in = calc_h(T_in, p, medium)
        Fi = flow * (h_out - h_in)
        self.states.h_in = h_in
        self.states.h_out = h_out
        self.states.T_in = T_in
        self.states.Fi = Fi
        self.states.dT = T_out-T_in


    def calculate(self, T_in, T_out, p_in, p_out, flow):
        p = max(p_in, p_out)
        self.states.flow = flow
        if T_in == -1 and T_out > 0.0:
            self.states.T_out = T_out
            self.temperature_in(p)
        elif T_out == -1 and T_in > 0.0:
            self.states.T_in = T_in
            self.temperature_out(p)
        else:
            #print("krive temperature - cijev!")
            #print(T_in)
            #print(T_out)
            pass
        if p_in == -1 and p_out > 0.0:
            self.states.p_out = p_out
            self.pressure_in()
        elif p_out == -1 and p_in > 0.0:
            self.states.p_in = p_in
            self.pressure_out()
        else:
            print("krivi tlakovi - cijev!")
        self.save.save_states(self.states)

    def no_flow(self, T_in, p_in):
        self.states.T_in = T_in
        self.states.T_out = T_in
        self.states.p_in = p_in
        self.states.p_out = p_in
        self.states.flow = 0.0 #kmol/s
        self.states.h_in = 0.0
        self.states.h_out = 0.0
        self.states.Fi = 0.0
        self.states.dT = 0.0
        self.states.M = 0.0
        self.states.ro = 0.0
        self.states.mu = 0.0
        self.states.dp = 0
        self.states.speed = 0
        self.states.friction_coeff = 0.0
        self.save.save_states(self.states)

def colebrook(f, roughness, d, speed, ro, mu):
    '''
    Proračun koeficijenta trenja prema colebrookovoj jednadžbi
    https://www.engineeringtoolbox.com/colebrook-equation-d_1031.html
    '''
    Re = speed * d * ro / mu
    diff = f**(-0.5) + 2 * np.log10( roughness / (d*3.72) + 2.51 / (Re * f**0.5))
    diff_sq = diff**2
    return diff_sq


def calc_h(T, p, medium):
    if medium == 'lng_liq':
        met_liq.set_T(T)
        h = met_liq.get_mol_enthalpy()*1000.0  #J/kmol
    elif medium == 'lng_vap':
        met_vap.set_pT(p, T)
        h = met_vap.get_mol_enthalpy() * 1000 #J/kmol
    elif medium == 'glycol':
        h = glyc.enthalpy(T) #J/kg
    else:
        h = -1
        print("PROBLEM calc_properties")
    return h



def T_out_fun(T_out, T_in, T_env, p, h_in, k, A, flow, medium):
    T_mean = (T_in + T_out) / 2.0
    Fi_T = k * A * (T_env - T_mean)
    h_out = calc_h(T_out, p, medium)
    Fi_h = flow * (h_out - h_in)
    res = (Fi_T - Fi_h) **2
    return res

def T_in_fun(T_in, T_out, T_env, p, h_out, k, A, flow, medium):
    T_mean = (T_in + T_out) / 2.0
    Fi_T = k * A * (T_env - T_mean)
    h_in = calc_h(T_in, p, medium)
    Fi_h = flow * (h_out - h_in)
    res = (Fi_T - Fi_h) **2
    return res
