# -*- coding: utf-8 -*-
"""
"""

from src.properties.vap_props_lowlev import MetVap
met_vap = MetVap()

class ValveParams:
    def __init__(self, Kv):
        self.Kv = Kv

class ValveStates:
    def __init__(self):
        self.p_in = 0.0 #Pa
        self.p_out = 1.0 #Pa
        self.gas_mol_flow = 0.0
        self.ro = 0.0 #kg/m3
        self.T_in = 0.0
        self.Hmol_in = 0.0
        self.T_out = 0.0




class ValveSave:
    '''
    Saving pipe states
    '''
    def __init__(self, states):
        '''
        https://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
        '''
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

class Valve:
    def __init__(self, params, states):
        self.params = params
        self.states = states
        self.save = ValveSave(self.states)
        self.save.save_states(self.states)

    def calculate(self, p_in, T_in, p_out, ro):
        self.p_in = p_in #Pa
        self.T_in = T_in #K
        self.p_out = p_out #Pa
        self.ro = ro #kg/m3
        dp = p_in -p_out #Pa
        dp_bar = dp/100000.0
        if dp <= 10.0:
            gas_flow = 0.0
            gas_mol_flow = 0.0
            self.states.T_out = T_in
        else:
            Kv = self.params.Kv
            gas_flow = Kv/0.032 * (dp_bar*ro)**0.5 #kg/h
            gas_mol_flow = gas_flow/3600 / 16.04 #kmol/s
            met_vap.set_pT(p_in, T_in)
            Hmol_in=met_vap.get_mol_enthalpy()
            met_vap.set_pHmol(p_out, Hmol_in)
            T_out = met_vap.get_T()
            self.states.T_out = T_out
        self.states.gas_mol_flow = gas_mol_flow

        self.save.save_states(self.states)
