# -*- coding: utf-8 -*-

from src.properties.liq_props_lowlev import MetLiq
M=16.04
met_liq = MetLiq()


class ThrottleValveStates:
    def __init__ (self):
        self.T_in = 0.0
        self.p_in = 0.0
        self.p_out = 0.0
        self.H_mol_in = 0.0
        self.T_out = 0.0


class ThrValSave:
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

    def save_states(self, states):
        #iteriraj po svim atributima
        for st in self.state_list:
            #trenutna vrijednost atributa
            current_state = getattr(states, st)
            # spremi u dictionary - dodaj zadnju vrijednost liste
            self.states[st].append(current_state)
            #print('aaa')

class ThrottleValve:
    def __init__(self,states):
        self.states = states
        self.save = ThrValSave(self.states)
        self.save.save_states(self.states)

    def calculate(self, T_in, p_in, p_out):
        # PROBLEM: KAD UPADNE U ZASIĆENO PODRUČJE
        self.states.T_in = T_in
        self.states.p_in = p_in
        self.states.p_out = p_out
        met_liq.set_pT_robust(p_in, T_in)
        H_in = met_liq.get_mol_enthalpy()  #kJ/kmol
        self.states.H_mol_in = H_in
        met_liq.set_pHmol_wet(p_out, H_in)
        T_out = met_liq.get_T()
        self.states.T_out = T_out
        self.save.save_states(self.states)
        
    def no_flow(self, T_in, p_in):
        self.states.T_in = T_in
        self.states.T_out = T_in
        self.states.p_in = p_in
        self.states.p_out = p_in
        self.states.H_mol_in = 0.0
        self.save.save_states(self.states)


    