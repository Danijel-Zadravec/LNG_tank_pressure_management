# -*- coding: utf-8 -*-
"""
"""
import numpy as np
import scipy.optimize as opt

from src.properties.liq_props_lowlev import MetLiq
from src.properties.vap_props_lowlev import MetVap
from src.properties.glycol_props_lowlev import  GlycolClass
from src.rezimi import Rezimi

met_liq = MetLiq()
met_vap = MetVap()
glyc = GlycolClass(0.5) #50 % maseni udio etilen glikola
rezim = Rezimi()


class EvapParams:
    def __init__(self, m, h, t, n, s, l, W, L):
        '''
        m - broj prolaza
        '''

        self.h = h #visina, m
        self.t = t #širina, m
        self.n = n #broj labela po metru 1/
        self.p = 1/n
        self.s = s # debljina pregrade,
        self.conductivity = l #toplinska provodnost aluminij, W/mK
        self.A_prim = 2*(1-n*t) #m2/m2 - po metru povrsine sloja str 47 ALPEMA - ukupna
        self.A_sec = 2*n*(h-t) # m2/m2 - sekundarna površina - ukupna
        self.W = W # širina izmjenjivača,
        self.L = L #duljina izmjenjivača,
        self.m = m #broj prolaza (svake struje), -
        self.d_ekv =2*h*self.p / (h + self.p) #ekvivalentni promjer
        self.An = self.m * W * h * (1-n*t) #nastrujna površina
        self.b = h/self.p #omjer stranica pravokutnika


class LngEvap():
    def __init__(self, T_in, p, qn):
        self.T_in = T_in
        self.p = p
        met_liq.set_pT(p, T_in)
        self.H_molar_in = met_liq.get_mol_enthalpy()*1000.0 #J/kmol
        met_liq.set_p(p)
        T_sat = met_liq.get_T()
        self.T_sat = T_sat
        met_vap.set_p_sat(p)
        self.H_molar_out = met_vap.get_mol_enthalpy()*1000.0 #J/kmol
        self.mol_flow = qn
        self.M = met_vap.get_M()
        self.mas_flow = qn * self.M
        self.Tm = (T_sat + T_in) / 2.0
        self.ro = 0.0
        self.vis_dyn = 0.0
        self.vis_cin = 0.0
        self.l = 0.0
        self.cp = 0.0
        self.Re = 0.0
        self.Nu = 0.0
        self.alfaLiq = 0.0

        self.hLV = 0.0 #toplina isparivanja
        self.G = 0.0 #gustoca masenog protoka, kg/m2s
        self.q = 0.0 #
        self.Bo = 0.0 #Boiling number
        self.alfa = 0.0




class LNGSuperh():
    def __init__(self, T_out, p, qn):
        self.T_out = T_out
        self.p = p
        met_vap.set_p_sat(p)
        T_sat = met_vap.get_T()
        self.T_sat = T_sat
        self.H_molar_in = met_vap.get_mol_enthalpy()
        met_vap.set_pT(p, T_out)
        self.H_molar_out = met_vap.get_mol_enthalpy()
        self.mol_flow = qn

        self.Tm =(T_out+T_sat)/2.0
        met_vap.set_pT(p, self.Tm)
        self.M = met_vap.get_M()
        self.mas_flow = qn * self.M
        self.ws = 0.0 #brzina zvuka, m/s
        self.ro = met_vap.get_density()
        self.nu = met_vap.get_kinematic_viscosity()
        self.cp = met_vap.get_cp()
        self.l=met_vap.get_lambda()
        self.Pr = self.nu / (self.l / (self.ro*self.cp))
        self.acf = 0.9 #akomodacijski faktor
        self.acef = 0.9 #akomodacijski faktor energije
        self.Re = 0.0
        self.Ma = 0.0
        self.Kn = 0.0
        self.kkl = 0.0 #koef klizanja
        self.kts = 0.0 #koef. temperaturnog skoka
        self.Nu = 0.0
        self.alfa = 0.0

class GlycolEvap():
    def __init__(self, T_in, T_out, qm):
        self.T_in = T_in
        self.T_out = T_out
        self.mas_flow = qm
        self.Tm = 0.0

        self.ro = 0.0 #kg/m3
        self.vis_dyn = 0.0
        self.vis_cin = 0.0
        self.l = 0.0
        self.Re = 0.0
        self.Nu = 0.0
        self.alfa = 0.0

class GlycolSuperh():
    def __init__(self,T_in, T_out, qm):
        self.T_in = T_in
        self.T_out = T_out
        self.mas_flow = qm
        self.Tm = 0.0

        self.ro = 0.0 #kg/m3
        self.vis_dyn = 0.0
        self.vis_cin = 0.0
        self.l = 0.0
        self.Re = 0.0
        self.Nu = 0.0
        self.alfa = 0.0


class EvapCommonStates():
    def __init__(self):
        self.A_fraction =1.0
        self.A = 0.0
        self.k = 1.0
        self.heat_flow = 1.0
        self.LMTD = 1.0
        self.heat_flow_LMTD = 0.0

class SuperhCommonStates():
    def __init__(self):
        self.A_fraction =0.0
        self.A = 0.0
        self.k = 1.0
        self.heat_flow = 1.0
        self.LMTD = 1.0


class EvapStates():
    def __init__(self, evap_com, evap_lng, evap_glyc, superh_com, superh_lng, superh_glyc):
        self.evap_com = evap_com
        self.evap_lng = evap_lng
        self.evap_glyc = evap_glyc
        self.superh_com = superh_com
        self.superh_lng = superh_lng
        self.superh_glyc = superh_glyc



class SaveEvap:
    '''
    Saving heat exchanger states
    '''
    def __init__(self, states):
        #atributi u listu
        self.evap_lng_st_list = [a for a in dir(states.evap_lng) if not a.startswith('__') and not callable(getattr(states.evap_lng, a))]
        #atributi u listu
        self.evap_glycol_st_list = [a for a in dir(states.evap_glyc) if not a.startswith('__') and not callable(getattr(states.evap_glyc, a))]
        #atributi zajednicki
        self.evap_com_st_list = [a for a in dir(states.evap_com) if not a.startswith('__') and not callable(getattr(states.evap_com, a))]

        self.superh_lng_st_list = [a for a in dir(states.superh_lng) if not a.startswith('__') and not callable(getattr(states.superh_lng, a))]
        #atributi u listu
        self.superh_glycol_st_list = [a for a in dir(states.superh_glyc) if not a.startswith('__') and not callable(getattr(states.superh_glyc, a))]
        #atributi zajednoicki
        self.superh_com_st_list = [a for a in dir(states.superh_com) if not a.startswith('__') and not callable(getattr(states.superh_com, a))]

        self.evap_lng = {i : [] for i in self.evap_lng_st_list}
        self.evap_glyc = {i : [] for i in self.evap_glycol_st_list}
        self.evap_com = {i : [] for i in self.evap_com_st_list}
        self.superh_lng = {i : [] for i in self.superh_lng_st_list}
        self.superh_glyc = {i : [] for i in self.superh_glycol_st_list}
        self.superh_com = {i : [] for i in self.superh_com_st_list}


    def save_states(self,states):
        #iteriraj po svim atributima
        for lng_st in self.evap_lng_st_list:
            #trenutna vrijednost atributa
            current_state = getattr(states.evap_lng, lng_st)
            # spremi u dictionary - dodaj zadnju vrijednost liste
            self.evap_lng[lng_st].append(current_state)

        for glycol_st in self.evap_glycol_st_list:
            current_state = getattr(states.evap_glyc, glycol_st)
            self.evap_glyc[glycol_st].append(current_state)

        for common_st in self.evap_com_st_list:
            current_state = getattr(states.evap_com, common_st)
            self.evap_com[common_st].append(current_state)

        for lng_st in self.superh_lng_st_list:
            current_state = getattr(states.superh_lng, lng_st)
            # spremi u dictionary - dodaj zadnju vrijednost liste
            self.superh_lng[lng_st].append(current_state)

        for glycol_st in self.superh_glycol_st_list:
            current_state = getattr(states.superh_glyc, glycol_st)
            self.superh_glyc[glycol_st].append(current_state)

        for common_st in self.superh_com_st_list:
            current_state = getattr(states.superh_com, common_st)
            self.superh_com[common_st].append(current_state)


class Evaporator():
    def __init__(self,params, states):
        self.params = params
        self.states = states
        self.save = SaveEvap(self.states) #inicijalizacija spremanja rezultata
        self.save.save_states(self.states) #spremanje pocetnih stanja

    def calc_super_glyc_states(self):
        T_in =self.states.superh_glyc.T_in
        qm = self.states.superh_glyc.mas_flow #kg/s
        Fi = self.states.superh_com.heat_flow #W
        glyc.set_T(T_in)
        cp0 = glyc.get_cp()
        T_out0 = T_in - Fi / (qm*cp0)
        T_m_prv = (T_in+T_out0)/2
        T_m = T_m_prv-1.0
        while (abs(T_m-T_m_prv) > 0.00001):
            glyc.set_T(T_m)
            cp = glyc.get_cp()
            T_out = T_in - Fi / (qm*cp)
            #print(T_out)
            T_m_prv = T_m
            T_m = (T_in+T_out)/2
        self.states.superh_glyc.T_out = T_out
        self.states.superh_glyc.Tm = T_m
        glyc.set_T(T_m)
        self.states.superh_glyc.c = glyc.get_cp() #J/kgK
        self.states.superh_glyc.l = glyc.get_lambda() #W/mK
        self.states.superh_glyc.vis_dyn = glyc.get_dynamic_viscosity() #Pas
        self.states.superh_glyc.ro = glyc.get_density() #kg/m3
        self.states.superh_glyc.vis_cin = glyc.get_kinematic_viscosity() #m2/s

    def calc_evap_glyc_states(self):
        T_in =self.states.evap_glyc.T_in
        qm = self.states.evap_glyc.mas_flow #kg/s
        Fi = self.states.evap_com.heat_flow #W
        glyc.set_T(T_in)
        cp0 = glyc.get_cp()
        T_out0 = T_in - Fi / (qm*cp0)
        T_m_prv = (T_in+T_out0)/2
        T_m = T_m_prv-1.0
        while (abs(T_m-T_m_prv) > 0.00001):
            glyc.set_T(T_m)
            cp = glyc.get_cp()
            T_out = T_in - Fi / (qm*cp)
            #print(T_out)
            #print(T_out)
            T_m_prv = T_m
            T_m = (T_in+T_out)/2
        self.states.evap_glyc.T_out = T_out
        self.states.evap_glyc.Tm = T_m
        glyc.set_T(T_m)
        self.states.evap_glyc.c = glyc.get_cp() #J/kgK
        self.states.evap_glyc.l = glyc.get_lambda() #W/mK
        self.states.evap_glyc.vis_dyn = glyc.get_dynamic_viscosity() #Pas
        self.states.evap_glyc.ro = glyc.get_density() #kg/m3
        self.states.evap_glyc.vis_cin = glyc.get_kinematic_viscosity() #m2/s

    def calc_evap_heat_flow(self):
        h_mol_out_lng = self.states.evap_lng.H_molar_out
        h_mol_in_lng = self.states.evap_lng.H_molar_in
        qn_lng = self.states.evap_lng.mol_flow
        Fi = qn_lng*(h_mol_out_lng-h_mol_in_lng) #W
        self.states.evap_com.heat_flow = Fi
    def calc_superh_heat_flow(self):
        h_mol_out_lng = self.states.superh_lng.H_molar_out
        h_mol_in_lng = self.states.superh_lng.H_molar_in
        qn_lng = self.states.superh_lng.mol_flow
        Fi = qn_lng*(h_mol_out_lng-h_mol_in_lng) #W
        self.states.superh_com.heat_flow = Fi

    def calc_superh_alfa_glyc(self):
        # T_in =self.states.superh_glyc.T_in #K
        # T_out =self.states.superh_glyc.T_out #K
        qm = self.states.superh_glyc.mas_flow #kg/s
        b = self.params.b #omjer stranica kanala,
        An = self.params.An #površina nastrujavanja (poprečni presjek), m2
        d_ekv = self.params.d_ekv #ekvivalentni promjer, m
        ro = self.states.superh_glyc.ro #gustoća glikola, kg/m3
        nu = self.states.superh_glyc.vis_cin #kinematička viskoznost glikola, m2/s
        lmbd = self.states.superh_glyc.l
        w = qm / ro / An #brzina strujanja, m/s
        Re = d_ekv * w / nu
        self.states.superh_glyc.Re = Re
        assert (Re<2300), "Nije laminarno - glikol!"
        # Heat transfer and fluid flow in minichannels and microchannels(2014)
        #str 132, (3.48) - lam. str,  grijana 2 suprotna zida
        Nu = 7.121 - 3.627*b + 1.145 * b**2 - 0.1689 * b**3 + 0.008464 * b**4
        self.states.superh_glyc.Nu = Nu
        alfa = Nu * lmbd / d_ekv
        self.states.superh_glyc.alfa = alfa

    def calc_evap_alfa_glyc(self):
        qm = self.states.evap_glyc.mas_flow #kg/s
        b = self.params.b #omjer stranica kanala,
        An = self.params.An #površina nastrujavanja (poprečni presjek), m2
        d_ekv = self.params.d_ekv #ekvivalentni promjer, m
        ro = self.states.evap_glyc.ro #gustoća glikola, kg/m3
        nu = self.states.evap_glyc.vis_cin #kinematička viskoznost glikola, m2/s
        lmbd = self.states.evap_glyc.l
        w = qm / ro / An #brzina strujanja, m/s
        Re = d_ekv * w / nu
        self.states.evap_glyc.Re = Re
        assert (Re<2300), "Nije laminarno - glikol!"
        # Heat transfer and fluid flow in minichannels and microchannels(2014)
        #str 132, (3.48) - lam. str,  grijana 2 suprotna zida
        Nu = 7.121 - 3.627*b + 1.145 * b**2 - 0.1689 * b**3 + 0.008464 * b**4
        self.states.evap_glyc.Nu = Nu
        alfa = Nu * lmbd / d_ekv
        self.states.evap_glyc.alfa = alfa


    def calc_superh_alfa_lng(self):
        T_in = self.states.superh_lng.T_sat #K
        T_out = self.states.superh_lng.T_out #K
        p = self.states.superh_lng.p
        Tm = (T_in + T_out) / 2.0
        self.states.superh_lng.Tm = Tm
        met_vap.set_pT(p, Tm)
        ro =  met_vap.get_density()
        self.states.superh_lng.ro = ro
        nu =  met_vap.get_kinematic_viscosity()
        self.states.superh_lng.nu = nu
        cp = met_vap.get_cp()
        self.states.superh_lng.cp =cp
        lmbd = met_vap.get_lambda()
        self.states.superh_lng.l = lmbd
        Pr = nu / (lmbd/(ro*cp))
        self.states.superh_lng.Pr = Pr

        qm = self.states.superh_lng.mas_flow
        An = self.params.An #površina nastrujavanja (poprečni presjek), m2
        d_ekv = self.params.d_ekv #ekvivalentni promjer, m
        molar_mas = self.states.superh_lng.M #kg/kmol
        acf = self.states.superh_lng.acf #akomodacijski faktor
        acef =  self.states.superh_lng.acef #akomodacijski faktor energije

        w = qm / ro / An #brzina strujanja, m/s
        Re = w * d_ekv / nu
        self.states.superh_lng.Re = Re

        R = 8314 / molar_mas #J/kgK
        cv = cp-R #J/kgK
        kappa = cp/cv
        ws = (kappa*R*Tm)**0.5 #brzina zvuka, m/s
        self.states.superh_lng.ws = ws
        Ma = w/ws
        self.states.superh_lng.Ma = Ma
        Kn = (2*kappa / np.pi) **0.5 * Ma / Re
        self.states.superh_lng.Kn = Kn
        kkl = 4 *Kn * (2-acf) / acf #koeficijent klizanja
        #print(kkl)
        self.states.superh_lng.kkl = kkl
        kts =8 * (2-acef)/acef * kappa/(kappa+1) * Kn/Pr #koeficijent temp. skoka
        #print(kts)
        self.states.superh_lng.kts = kts
        Nu = 1.0 / (kts/4.0 + (17.0 + 84.0*kkl + 105.0 * kkl**2)/(140.0 * (1.0 + 3.0*kkl)**2 ) )
        #print(Nu)
        self.states.superh_lng.Nu = Nu
        alfa = Nu * lmbd / d_ekv
        #print(alfa)
        self.states.superh_lng.alfa = alfa

    def calc_evap_alfaLiq_lng(self):
         qm = self.states.evap_lng.mas_flow #kg/s
         b = self.params.b #omjer stranica kanala,
         An = self.params.An #površina nastrujavanja (poprečni presjek), m2
         d_ekv = self.params.d_ekv #ekvivalentni promjer, m
         ro = self.states.evap_lng.ro #gustoća glikola, kg/m3
         nu = self.states.evap_lng.vis_cin #kinematička viskoznost glikola, m2/s
         lmbd = self.states.evap_lng.l
         w = qm / ro / An #brzina strujanja, m/s
         Re = d_ekv * w / nu
         assert (Re<1600.0), "Re < 100 - samo ovaj model je implementiran za sad"
         self.states.evap_lng.Re = Re
         Nu = 7.121 - 3.627*b + 1.145 * b**2 - 0.1689 * b**3 + 0.008464 * b**4
         self.states.evap_lng.Nu = Nu
         alfa = Nu * lmbd / d_ekv
         self.states.evap_lng.alfaLiq = alfa


    def calc_evap_alfa_lng(self):
        An = self.params.An #površina nastrujavanja (poprečni presjek), m2
        qm = self.states.evap_lng.mas_flow #kg/s
        G = qm/An #gustoca masenog toka, kg/m2s
        self.states.evap_lng.G = G
        Fi = self.states.evap_com.heat_flow
        A_prim = self.params.A_prim
        A_sec = self.params.A_sec
        W = self.params.W #širina izmjenjiivaca, m
        L_tot = self.params.L #duljina izmjenjivača, m
        t=self.params.t
        h = self.params.h
        evap_frac = self.states.evap_com.A_fraction
        lmbd = self.states.evap_lng.l
        alfa = 10000.0
        alfa_prev = 0.0
        i=0
        while abs(alfa-alfa_prev) > 1.0:
            beta = h*(2*alfa/(lmbd*t))**0.5
            eta =  np.tanh(beta/2)*2/beta
            A_wall = evap_frac * (A_prim*W*L_tot + A_sec*eta*W*L_tot)
            q = Fi / A_wall
            hLV = self.states.evap_lng.hLV
            Bo = q/(G*hLV)
            alfaLiq = self.states.evap_lng.alfaLiq
            alfa_prev = alfa
            alfa = 1058.0 * (Bo**0.7) * alfaLiq
            i=i+1
            #print("i: " + str(i) + "  alfa: " + str(alfa))

        self.states.evap_lng.q = q
        self.states.evap_lng.Bo = Bo #boiling number
        self.states.evap_lng.alfa = alfa




    def calc_area_superh(self):
        T_lng_in = self.states.superh_lng.T_sat
        T_lng_out = self.states.superh_lng.T_out
        T_glyc_in =self.states.superh_glyc.T_in
        T_glyc_out = self.states.superh_glyc.T_out
        Fi_superh = self.states.superh_com.heat_flow
        #print(T_glyc_out)
        alpha_glyc = self.states.superh_glyc.alfa
        alpha_lng = self.states.superh_lng.alfa
        A_prim = self.params.A_prim
        A_sec = self.params.A_sec
        h = self.params.h
        t=self.params.t
        s = self.params.s
        m = self.params.m
        lmbd = self.params.conductivity
        W = self.params.W #širina izmjenjiivaca, m
        L_tot = self.params.L #duljina izmjenjivača, m
        beta_glyc = h*(2*alpha_glyc/(lmbd*t))**0.5
        eta_glyc = np.tanh(beta_glyc/2)*2/beta_glyc
        A_glyc_ = A_prim + eta_glyc*A_sec
        A_glycL = A_glyc_ * W * m
        beta_lng = h*(2*alpha_lng/(lmbd*t))**0.5
        eta_lng =  np.tanh(beta_lng/2)*2/beta_lng
        A_lng_ = A_prim + eta_lng*A_sec
        A_lngL = A_lng_ * W  * m
        thick = s+0.5*t
        R_al = thick/(lmbd*W*A_prim)
        R_lng =1/(alpha_lng*A_lngL)
        R_glyc =1/(alpha_glyc*A_glycL)
        R_sumL =R_al + R_lng + R_glyc

        dT_max = T_glyc_out - T_lng_in
        dT_min = T_glyc_in - T_lng_out

        assert (dT_max > dT_min), "dT_max < dT_min"
        LMTD = (dT_max-dT_min)/np.log(dT_max/dT_min)
        L_superh = Fi_superh * R_sumL / LMTD
        L_superh = min(L_superh, L_tot-0.00001)
        assert(L_tot > L_superh)
        A_superh_fr = L_superh/L_tot
        self.states.superh_com.A_fraction = A_superh_fr
        self.states.evap_com.A_fraction = 1.0 - A_superh_fr
        #print(LMTD)

    def calc_evap_LMTD_heat_flow(self):
        T_lng_in = self.states.evap_lng.T_in
        T_lng_out = self.states.evap_lng.T_sat
        T_glyc_in =self.states.evap_glyc.T_in
        T_glyc_out = self.states.evap_glyc.T_out
        #print(T_glyc_out)
        alpha_glyc = self.states.evap_glyc.alfa
        alpha_lng = self.states.evap_lng.alfa
        A_prim = self.params.A_prim
        A_sec = self.params.A_sec
        h = self.params.h
        t=self.params.t
        s = self.params.s
        m = self.params.m
        evap_frac = self.states.evap_com.A_fraction
        lmbd = self.params.conductivity
        W = self.params.W #širina izmjenjiivaca, m
        L_total = self.params.L #duljina izmjenjivača, m
        beta_glyc = h*(2*alpha_glyc/(lmbd*t))**0.5
        eta_glyc = np.tanh(beta_glyc/2)*2/beta_glyc
        A_glyc_ = A_prim + eta_glyc*A_sec
        A_glyc = A_glyc_ * W * L_total * m * evap_frac

        beta_lng = h*(2*alpha_lng/(lmbd*t))**0.5
        eta_lng =  np.tanh(beta_lng/2)*2/beta_lng
        A_lng_ = A_prim + eta_lng*A_sec
        A_lng = A_lng_ * W * L_total * m * evap_frac

        thick = s+0.5*t
        R_al = thick/(lmbd*W*L_total* evap_frac* A_prim)
        R_lng =1/(alpha_lng*A_lng)
        R_glyc =1/(alpha_glyc*A_glyc)

        R_sum =R_al + R_lng + R_glyc
        T_sat = self.states.evap_lng.T_sat
        dT_min = T_glyc_out - T_sat
        dT_max = T_glyc_in - T_sat
        #print([dT_max, dT_min])

        assert (dT_max > dT_min), "dT_max < dT_min"
        LMTD = (dT_max-dT_min)/np.log(dT_max/dT_min)
        self.states.evap_com.LMTD = LMTD
        Fi_LMTD = LMTD/R_sum
        self.states.evap_com.heat_flow_LMTD = Fi_LMTD

    def T_lng_out_fun(self, T_lng_out):
        self.states.superh_lng.T_out = T_lng_out
        p_lng = self.states.evap_lng.p
        met_vap.set_pT(p_lng, T_lng_out)
        self.states.superh_lng.H_molar_out = met_vap.get_mol_enthalpy()*1000.0 #J/kmol
        self.calc_superh_heat_flow()


        #glycol superheating part
        self.calc_super_glyc_states()
        self.calc_superh_alfa_glyc()
        self.calc_superh_alfa_lng()
        self.calc_area_superh()

        #glycol evaporation part
        self.states.evap_glyc.T_in = self.states.superh_glyc.T_out
        self.calc_evap_glyc_states()
        self.calc_evap_alfa_glyc()
        self.calc_evap_alfaLiq_lng()
        self.calc_evap_alfa_lng()
        self.calc_evap_LMTD_heat_flow()
        Fi = self.states.evap_com.heat_flow
        Fi_LMTD = self.states.evap_com.heat_flow_LMTD
        return (Fi-Fi_LMTD)**2



    def update_states(self, T_lng_in, T_lng_out0, p_lng, qn_lng, T_glyc_in, qm_glyc):
        if (qn_lng > 0.0):
            #lng evaporation
            self.states.evap_lng.T_in = T_lng_in
            self.states.evap_lng.p = p_lng
            self.states.evap_lng.mol_flow = qn_lng
            mas_flow_lng = qn_lng * self.states.superh_lng.M
            self.states.evap_lng.mas_flow = mas_flow_lng
            met_liq.set_p(p_lng)
            T_sat = met_liq.get_T()
            T_lng_in= min(T_sat-0.1, T_lng_in)  #TREBA NAPRAVITI S ULAZOM MOKRE PARE!!!
            met_liq.set_pT(p_lng, T_lng_in)
            self.states.evap_lng.H_molar_in = met_liq.get_mol_enthalpy()*1000.0 #J/kmol
            met_liq.set_p(p_lng)
            hL = met_liq.get_mas_enthalpy()
            T_sat = met_liq.get_T()
            self.states.evap_lng.T_sat = T_sat
            met_vap.set_p_sat(p_lng)
            hV = met_vap.get_mass_enthalpy()
            self.states.evap_lng.H_molar_out = met_vap.get_mol_enthalpy()*1000.0 #J/kmol
            self.calc_evap_heat_flow()
            Tm_lng_evap = (T_lng_in + T_sat) / 2.0
            self.states.evap_lng.Tm = Tm_lng_evap
            met_liq.set_T(Tm_lng_evap)
            self.states.evap_lng.ro = met_liq.get_density()
            self.states.evap_lng.vis_dyn = met_liq.get_dynamic_viscosity()
            self.states.evap_lng.vis_cin = met_liq.get_kinematic_viscosity()
            self.states.evap_lng.l = met_liq.get_lambda()
            self.states.evap_lng.cp = met_liq.get_cp()
            hLV = hV-hL
            self.states.evap_lng.hLV = hLV

            # lng superheating
            self.states.superh_lng.T_sat = T_sat
            self.states.superh_lng.H_molar_in = met_vap.get_mol_enthalpy()*1000.0 #J/kmol
            self.states.superh_lng.p = p_lng
            self.states.superh_lng.mol_flow = qn_lng
            self.states.superh_lng.mas_flow = mas_flow_lng

            self.states.superh_glyc.T_in = T_glyc_in
            self.states.superh_glyc.mas_flow = qm_glyc

            self.states.evap_glyc.mas_flow = qm_glyc


            opt_rez = opt.minimize_scalar(self.T_lng_out_fun, bounds=(T_sat+0.00001, T_glyc_in-0.0000001), method='bounded')
            #T_superheating = opt_rez.x
            #print(T_superheating)
            self.save.save_states(self.states) #spremanje pocetnih stanja
        elif (qn_lng == 0.0):
            self.no_flow(T_lng_in, p_lng, T_glyc_in)
        else:
            assert False, "negative flow evap"

    def no_flow(self,T_lng_in, p_lng, T_glyc_in):
        self.save.save_states(self.states) #spremanje pocetnih stanja
