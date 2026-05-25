# -*- coding: utf-8 -*-
"""
"""


from src.properties.vap_props_lowlev import MetVap
met_vap = MetVap()

class MixerVstates:
    def __init__(self):
        self.T_in1 = 0
        self.T_in2 = 0
        self.T_in3 = 0

        self.qn1 = 0
        self.qn2 = 0
        self.qn3 = 0

        self.p = 0
        self.T_out = 0.0
        self.Hmol1 = 0.0
        self.Hmol2 = 0.0
        self.Hmol3 = 0.0
        self.qn_tot = 0.0
        self.Htot = 0.0
        self.Hmol_out = 0.0


class MixerVsave:
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


class MixerV:
    def __init__(self, states):
        self.states = states
        self.save = MixerVsave(self.states)
        self.save.save_states(self.states)

    def calculate(self, T_in1, T_in2, T_in3, qn1, qn2, qn3, p):
        self.T_in1 = T_in1
        self.T_in2 = T_in2
        self.T_in3 = T_in3
        self.qn1 = qn1
        self.qn2 = qn2
        self.qn3 = qn3
        self.p = p
        qn_tot = qn1 + qn2 + qn3
        self.states.qn_tot = qn_tot
        met_vap.set_pT(p, T_in1)
        Hmol1 = met_vap.get_mol_enthalpy()
        self.states.Hmol1 = Hmol1
        if qn2>0.0:
            met_vap.set_pT(p, T_in2)
            Hmol2 = met_vap.get_mol_enthalpy()
        else:
            Hmol2 = 0.0
        if qn3>0.0:
            met_vap.set_pT(p, T_in3)
            Hmol3 = met_vap.get_mol_enthalpy()
        else:
            Hmol3 = 0.0
        self.states.Hmol2 = Hmol2
        self.states.Hmol3 = Hmol3

        Htot = qn1*Hmol1 + qn2*Hmol2 + qn3*Hmol3
        self.states.Htot=Htot
        if qn_tot > 0:
            Hmol_out = Htot/qn_tot
            met_vap.set_pHmol(p, Hmol_out)
            T_out = met_vap.get_T()
        else:
            Hmol_out = 0.0
            T_out = 0.0
        self.states.Hmol_out = Hmol_out
        self.states.T_out = T_out
        #print(T_out)
        self.save.save_states(self.states)

    def ne_radi(self):
        self.T_in1 = 0
        self.T_in2 = 0
        self.T_in3 = 0
        self.qn1 = 0
        self.qn2 = 0
        self.qn3 = 0
        self.p = 0
        self.T_out = 0.0
        self.Hmol1 = 0.0
        self.Hmol2 = 0.0
        self.Hmol3 = 0.0
        self.qn_tot = 0.0
        self.Htot = 0.0
        self.Hmol_out = 0.0

        self.save.save_states(self.states)