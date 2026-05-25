# -*- coding: utf-8 -*-
"""
"""

from src.properties.vap_props_lowlev import MetVap
met_vap = MetVap()

class CompressorStates:
    def __init__(self, T_in, T_out, p_in, p_out):
        self.T_in = T_in
        self.T_out = T_out
        self.p_in = p_in
        self.p_out = p_out
        self.H_mol_in =0.0
        self.S_mol_in = 0.0
        self.isentropic_H = 0.0
        self.H_mol_out = 0.0


class CompressorSave:
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

class Compressor:
    def __init__(self, eta, states):
        self.eta = eta
        self.states = states
        self.save = CompressorSave(self.states)
        self.save.save_states(self.states)

    def calculate(self,T_in, p_in, p_out):
        self.states.T_in = T_in
        self.states.p_in = p_in
        self.states.p_out = p_out
        met_vap.set_pT_robust(p_in, T_in)
        H_mol_in = met_vap.get_mol_enthalpy()
        self.states.H_mol_in = H_mol_in
        S_mol_in = met_vap.get_mol_entropy()
        self.states.S_mol_in = S_mol_in
        #izentropa

        met_vap.set_pSmol(p_out, S_mol_in)

        h_isentropic = met_vap.get_mol_enthalpy()

        self.states.isentropic_H = h_isentropic
        eta = self.eta
        H_mol_out =(h_isentropic-H_mol_in) / eta + H_mol_in
        self.states.H_mol_out = H_mol_out

        met_vap.set_pHmol(p_out, H_mol_out)
        T_out = met_vap.get_T()

        self.states.H_mol_out = H_mol_out
        self.states.T_out = T_out
        self.save.save_states(self.states)

    def ne_radi(self):
        self.T_in = 0
        self.T_out = 0
        self.p_in = 0
        self.p_out = 0
        self.H_mol_in =0.0
        self.S_mol_in = 0.0
        self.isentropic_H = 0.0
        self.H_mol_out = 0.0
        self.save.save_states(self.states)
