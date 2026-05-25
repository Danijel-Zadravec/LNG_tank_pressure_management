# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import scipy.optimize as opt
import src.properties.vap_props as vap
import src.properties.liq_props as liq
import src.properties.glycol_props as glyc

class HEparams:
    def __init__(self, area, k):
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
        self.A = area
        self.k = k

class HE_stream:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, T_out, mol_flow, p, stream):
        self.stream = stream
        self.T_in = T_in #K
        self.T_out = T_out #K
        self.mol_flow = mol_flow #kmol/s
        self.p = p
        self.Cmp_in = Cmp_stream(T_in, p, stream)
        self.Cmp_out = Cmp_stream(T_out, p, stream)



class StatesHE:
    def __init__(self, lng, glycol):
        self.lng = lng
        self.glycol = glycol
        self.heat_flow = self.calc_heat_flow()

    def calc_heat_flow(self):
         T_out_l = self.lng.T_out
         T_in_l = self.lng.T_in
         Cmp_in = self.lng.Cmp_in
         Cmp_out = self.lng.Cmp_out
         flow_l = self.lng.mol_flow
         heat_flow = (T_out_l*Cmp_out - T_in_l*Cmp_in) * flow_l
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

        self.lng = {i : [] for i in self.lng_st_list}
        self.glycol = {i : [] for i in self.glycol_st_list}

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
        self.states.lng.Cmp_in = Cmp_stream(lng_T_in, p, 'lng')
        self.states.lng.T_out = lng_T_out
        self.states.lng.Cmp_out = Cmp_stream(lng_T_out, p, 'lng')
        self.states.lng.mol_flow = lng_flow
        self.states.heat_flow = self.states.calc_heat_flow()

        self.states.glycol.T_in = glycol_T_in
        self.states.glycol.Cmp_in = Cmp_stream(glycol_T_in, -1, 'glycol')
        self.calc_glycol2()
        self.save.save_states(self.states) #spremanje pocetnih stanja


    def calc_glycol2(self):
        T_out_l = self.states.lng.T_out
        T_in_l = self.states.lng.T_in
        T_in_gl =self.states.glycol.T_in
        #flow_l = self.lng.mol_flow
        heat_flow = self.states.heat_flow
        k = self.params.k
        A=self.params.A
        lmtd = heat_flow / A / k
        solve_res = opt.minimize_scalar(T_gly_out, bounds=(T_in_l, T_in_gl), args=(T_in_gl, T_in_l, T_out_l, lmtd), method='bounded')
        T_out_gl = solve_res.x

        self.states.glycol.T_out = T_out_gl
        Cmp_glyc_out = Cmp_stream(T_out_gl, -1, 'glycol')
        self.states.glycol.Cmp_out = Cmp_glyc_out
        Cmp_glyc_in =self.states.glycol.Cmp_in
        flow_glyc = heat_flow / (Cmp_glyc_in*T_in_gl - Cmp_glyc_out*T_out_gl)
        self.states.glycol.mol_flow = flow_glyc





    def calc_glycol(self):
        T_out_l = self.states.lng.T_out
        T_in_l = self.states.lng.T_in
        T_in_gl =self.states.glycol.T_in
        #flow_l = self.lng.mol_flow
        heat_flow = self.states.heat_flow
        C_lng = heat_flow / (T_out_l - T_in_l) #toplinski kapacitet glikola, kJ/K
        k = self.params.k
        A=self.params.A
        pi1 = (T_in_l - T_out_l) / (T_in_l - T_in_gl) # pi 1 parametar
        print(" ")
        print("pi1: "+str(pi1))
        pi2 = k * A / C_lng # pi 2 parametar
        print(" ")
        print("pi2: "+str(pi2))
        opt_rez = opt.minimize_scalar(pi_3_func, bounds=(0, 1), args=(pi1, pi2), method='bounded')
        pi3 = opt_rez.x # pi 3 parametar
        if pi3 < 1.0:
            pass
        else:
            print("pi3>1!!!")
        print(" ")
        print("pi3: "+str(pi3))
        C_glyc = C_lng / pi3 # toplinski kapacitet glikola
        T_out_gl = T_in_gl - heat_flow / C_glyc;
        self.states.glycol.T_out = T_out_gl
        Cmp_glyc_out = Cmp_stream(T_out_gl, -1, 'glycol')
        self.states.glycol.Cmp_out = Cmp_glyc_out
        Cmp_glyc_in =self.states.glycol.Cmp_in
        flow_glyc = heat_flow / (Cmp_glyc_in*T_in_gl - Cmp_glyc_out*T_out_gl)
        self.states.glycol.mol_flow = flow_glyc


def pi_3_func(pi3, pi1, pi2):
    lhs = pi1
    exponent = np.exp(-(1-pi3)*pi2)
    rhs = (1-exponent) / (1 - pi3*exponent)
    error = (rhs-lhs)**2
    return error


def T_gly_out(T_g_out, T_g_in, T_lng_in, T_lng_out, lmtd):
    dt1 = T_g_out - T_lng_in
    dt2 = T_g_in - T_lng_out
    lmtd_calc = (dt1 - dt2) / (np.log(dt1/dt2))
    result = (lmtd_calc-lmtd) **2
    return result

def Cmp_stream(T, p, stream):
    '''
    Calculate heat capacity of a stream

    Parameters
    ----------
    T : TYPE
        DESCRIPTION.
    stream : TYPE
        DESCRIPTION.

    Returns
    -------
    Cmp : J/kmol!
        DESCRIPTION.

    '''
    if stream == 'lng':
        Cmp = vap.Cmp(T, p)*1000.0
    elif stream == 'glycol':
        Cmp = glyc.Cmp(T, -1)*1000.0
    else:
        Cmp = -1.0
    return Cmp


if __name__ == "__main__":
    #primjer koristenja:
    area = 10 #m2
    k = 1000 #W/m2K
    he_params = HEparams(area, k)
    p = 400000 #Pa
    T_lng_in = 293.0 #K
    T_lng_out = 293.0 #K
    lng_flow = 0.0 #kmol/s
    stream = 'lng'
    he_lng =HE_stream(T_lng_in, T_lng_out, lng_flow,  p, stream)

    T_glycol_in = 318.0
    T_glycol_out = 318.0
    glycol_flow = 0.0
    stream = 'glycol'

    he_glycol =HE_stream(T_glycol_in, T_glycol_out, glycol_flow, p, stream)
    he_states = StatesHE(he_lng, he_glycol)

    he = HeatExchanger(he_params, he_states)
    lng_T_in = liq.saturation_temperature(p)+2
    he.update_states(0.11, lng_T_in, T_lng_out, T_glycol_in, p)

    print(he.states.glycol.T_out)
    print(he.states.glycol.mol_flow)


