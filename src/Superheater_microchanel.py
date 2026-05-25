# -*- coding: utf-8 -*-


import numpy as np
import scipy.optimize as opt
from src.properties.glycol_props_lowlev import  GlycolClass
from src.properties.vap_props_lowlev import MetVap, MetVap_slow
from src.rezimi import Rezimi

met_vap = MetVap()
glyc = GlycolClass(0.5) #50 % maseni udio etilen glikola
rezim = Rezimi()



class HEparams:
    def __init__(self, m, h, t, n, s, l, W, L):
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

        self.h = h #visina, m
        self.t = t #širina, m
        self.n = n #broj labela po metru 1/
        self.p = 1/n
        self.s = s # debljina pregrade,
        self.conductivity = l #toplinska provodnost aluminij, W/
        self.A_prim = 2*(1-n*t) #m2/m2 - po metru povrsine sloja str 47 ALPEMA
        self.A_sec = 2*n*(h-t) # m2/m2 - sekundarna površina
        self.Nu_gas = 140./17.
        self.W = W # širina izmjenjivača,
        self.L = L #duljina izmjenjivača,
        self.m = m #broj prolaza (svake struje), -
        self.d_ekv =2*h*self.p / (h + self.p) #ekvivalentni promjer
        self.An = self.m * W * h * (1-n*t) #nastrujna površina
        self.b = h/self.p #omjer stranica pravokutnika
        b=self.b
        self.Nu_gl = 8.464*10**(-3)*b**4-0.1689*b**3+1.145*b**2-3.627*b+7.121


class LngStream:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, T_out0, mol_flow, p):

        self.T_in = T_in #K
        self.T_out = T_out0 #K
        self.p = p
        self.mol_flow = mol_flow #kmol/s
        #print(mol_flow)
        met_vap.set_pT(p, T_in)
        self.H_molar_in = met_vap.get_mol_enthalpy()*1000.0 #J/kmol
        met_vap.set_pT(p, T_out0)
        self.H_molar_out =  met_vap.get_mol_enthalpy()*1000.0 #J/kmol
        self.Tm = (T_in+T_out0)/2.0
        met_vap.set_pT(p, self.Tm)
        self.M = met_vap.get_M()
        self.mas_flow = mol_flow * self.M
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



class GlycolStream:
    '''
    Struja u izmjenjivaču
    '''
    def __init__(self, T_in, mas_flow):

        self.T_in = T_in #K
        self.T_out = 0.0
        self.T_m = 0.0
        self.mas_flow = mas_flow #kg/s
        glyc.set_T(T_in)
        self.c = glyc.get_cp() #J/kgK
        self.l = glyc.get_lambda() #W/mK
        self.vis_dyn = glyc.get_dynamic_viscosity() #Pas
        self.ro = glyc.get_density() #kg/m3
        self.vis_cin = glyc.get_kinematic_viscosity() #m2/s
        self.Re = 0.0
        self.Nu = 0.0
        self.alfa = 0.0



class StatesCommonHE:
    def __init__(self):
        self.pi1=0.0
        self.pi2=0.0
        self.pi3=0.0
        self.A = 0.0
        self.heat_flow = 0.0
        self.heat_flow_LMTD = -100.0

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
        heat_flow = (H_molar_out - H_molar_in) * flow_l #W
        #print(heat_flow)
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

    def no_flow(self, lng_T_in, lng_p, glycol_T_in):
        # self.states.lng = LngStream(lng_T_in+50, lng_T_in+50, 0.0, lng_p)
        # self.states.glycol = GlycolStream(glycol_T_in, 0.0)
        # self.states.common = StatesCommonHE()
        self.save.save_states(self.states) #spremanje pocetnih stanja



    def update_states(self, lng_flow, lng_T_in, lng_T_out0, lng_p, glycol_flow, glycol_T_in):
        if (lng_flow > 0.0):
            self.states.lng.T_in = lng_T_in
            self.states.lng.p = lng_p
            met_vap.set_pT(lng_p, lng_T_in)
            self.states.lng.H_molar_in =met_vap.get_mol_enthalpy()*1000.0
            self.states.lng.T_out = lng_T_out0 #pretpostavka
            met_vap.set_pT(lng_p, lng_T_out0)
            self.states.lng.H_molar_out = met_vap.get_mol_enthalpy()*1000.0 #pretpostavka
            self.states.lng.mol_flow = lng_flow
            self.states.lng.mas_flow = lng_flow * self.states.lng.M
            self.states.common.heat_flow = self.states.calc_heat_flow()
            self.states.glycol.mas_flow = glycol_flow
            self.states.glycol.T_in = glycol_T_in
            #self.calculate_HE()
            self.calc_HE_opt()
            self.save.save_states(self.states) #spremanje pocetnih stanja

            # #self.states.glycol.h_mas_in =  glyc.enthalpy(glycol_T_in)
            ##self.calc_glycol()
            ##self.save.save_states(self.states) #spremanje pocetnih
        elif (lng_flow == 0.0):
            self.no_flow(lng_T_in, lng_p, glycol_T_in)
        else:
            assert False, "superh neg flow"

    def calc_HE_opt(self):
        T_glycol_in = self.states.glycol.T_in
        T_lng_in = self.states.lng.T_in
        opt_rez = opt.minimize_scalar(self.HE_opt, bounds=(T_lng_in+0.00001, T_glycol_in-0.0000001), method='bounded')
        T_lng_out = opt_rez.x
        p=self.states.lng.p
        met_vap.set_pT(p, T_lng_out)
        heat_flow = self.states.calc_heat_flow()
        self.states.common.heat_flow = heat_flow
        self.calc_T_out_glyc()
        #print('')
        #print('lng T_out: ' + str(opt_rez.x))
        #print('glyc T_out: ' + str(self.states.glycol.T_out))



    def HE_opt(self, T_out):
        self.states.lng.T_out = T_out
        p=self.states.lng.p
        met_vap.set_pT(p, T_out)
        self.states.lng.H_molar_out =  met_vap.get_mol_enthalpy()*1000.0 #J/kmol
        heat_flow = self.states.calc_heat_flow()
        self.states.common.heat_flow = heat_flow

        self.calc_T_out_glyc()
        self.calc_alfa_glyc()
        self.calc_alfa_lng()
        self.LMTD_heat_flow()
        heat_flow_LMTD = self.states.common.heat_flow_LMTD
        return(heat_flow-heat_flow_LMTD)**2


    def calculate_HE(self):
        ii=0
        T_glycol_in = self.states.glycol.T_in
        while abs(self.states.common.heat_flow - self.states.common.heat_flow_LMTD) > 1.0:
            self.calc_T_out_glyc()
            self.calc_alfa_glyc()
            self.calc_alfa_lng()
            self.LMTD_heat_flow()
            heat_flow = (self.states.common.heat_flow * 0.9 + self.states.common.heat_flow_LMTD * 0.1)
            self.states.common.heat_flow = heat_flow
            lng_mol_flow = self.states.lng.mol_flow
            H_molar_out = self.states.lng.H_molar_in + (heat_flow/lng_mol_flow)
            self.states.lng.H_molar_out = H_molar_out
            print(self.states.common.heat_flow)
            met_vap.set_Hmolar_p(H_molar_out, self.states.lng.p)
            T_lng_out = met_vap.get_T()
            if T_lng_out < T_glycol_in:
                print(T_lng_out)
            else:
                print('prevvisoka izlazna temp lng ' + str(T_lng_out))
                T_lng_out = T_glycol_in - 0.0001
            self.states.lng.T_out = T_lng_out
            ii=ii+1
            if ii > 2500:
                break
        self.calc_T_out_glyc()
        self.calc_alfa_glyc()
        self.calc_alfa_lng()
        self.LMTD_heat_flow()



    def LMTD_heat_flow(self):
        T_lng_in = self.states.lng.T_in
        T_lng_out = self.states.lng.T_out
        T_glyc_in =self.states.glycol.T_in
        T_glyc_out = self.states.glycol.T_out
        #print(T_glyc_out)
        alpha_glyc = self.states.glycol.alfa
        alpha_lng = self.states.lng.alfa
        A_prim = self.params.A_prim
        A_sec = self.params.A_sec
        h = self.params.h
        t=self.params.t
        s = self.params.s
        m = self.params.m
        lmbd = self.params.conductivity
        W = self.params.W #širina izmjenjiivaca, m
        L = self.params.L #duljina izmjenjivača, m
        beta_glyc = h*(2*alpha_glyc/(lmbd*t))**0.5
        eta_glyc = np.tanh(beta_glyc/2)*2/beta_glyc
        A_glyc_ = A_prim + eta_glyc*A_sec
        A_glyc = A_glyc_ * W * L * m
        #print(A_glyc)
        # print("")
        # print(beta_glyc)
        # print(eta_glyc)
        # print(A_glyc)
        # print("")
        beta_lng = h*(2*alpha_lng/(lmbd*t))**0.5
        eta_lng =  np.tanh(beta_lng/2)*2/beta_lng
        A_lng_ = A_prim + eta_lng*A_sec
        A_lng = A_lng_ * W * L * m
        # print(beta_lng)
        # print(eta_lng)
        # print(A_lng)
        thick = s+0.5*t
        R_al = thick/(lmbd*W*L*A_prim)
        R_lng =1/(alpha_lng*A_lng)
        R_glyc =1/(alpha_glyc*A_glyc)
        # print("")
        # print(R_al)
        # print(R_lng)
        # print(R_glyc)
        R_sum =R_al + R_lng + R_glyc

        dT_max = T_glyc_out - T_lng_in
        dT_min = T_glyc_in - T_lng_out
        #print([dT_max, dT_min])

        assert (dT_max > dT_min), "dT_max < dT_min"
        LMTD = (dT_max-dT_min)/np.log(dT_max/dT_min)
        #print(LMTD)
        Fi_LMTD = LMTD/R_sum
        #print(Fi_LMTD)
        self.states.common.heat_flow_LMTD = Fi_LMTD





    def calc_alfa_glyc(self):
        T_in =self.states.glycol.T_in #K
        T_out =self.states.glycol.T_out #K
        qm = self.states.glycol.mas_flow #kg/s
        b = self.params.b #omjer stranica kanala, -
        An = self.params.An #površina nastrujavanja (poprečni presjek), m2
        d_ekv = self.params.d_ekv #ekvivalentni promjer, m
        ro = self.states.glycol.ro #gustoća glikola, kg/m3
        nu = self.states.glycol.vis_cin #kinematička viskoznost glikola, m2/s
        lmbd = self.states.glycol.l
        w = qm / ro / An #brzina strujanja, m/s
        Re = d_ekv * w / nu
        self.states.glycol.Re = Re
        #print(Re)
        assert (Re<2300), "Nije laminarno - glikol!"
        # Heat transfer and fluid flow in minichannels and microchannels(2014)
        #str 132, (3.48) - lam. str,  grijana 2 suprotna zida
        Nu = 7.121 - 3.627*b + 1.145 * b**2 - 0.1689 * b**3 + 0.008464 * b**4
        #print(Nu)
        self.states.glycol.Nu = Nu
        alfa = Nu * lmbd / d_ekv
        #print(alfa)
        self.states.glycol.alfa = alfa

    def calc_alfa_lng(self):
        T_in = self.states.lng.T_in #K
        T_out = self.states.lng.T_out #K
        p = self.states.lng.p
        Tm = (T_in + T_out) / 2.0
        self.states.lng.Tm = Tm
        met_vap.set_pT(p, Tm)
        ro =  met_vap.get_density()
        self.states.lng.ro = ro
        nu =  met_vap.get_kinematic_viscosity()
        self.states.lng.nu = nu
        cp = met_vap.get_cp()
        self.states.lng.cp =cp
        lmbd = met_vap.get_lambda()
        self.states.lng.l = lmbd
        Pr = nu / (lmbd/(ro*cp))
        self.states.lng.Pr = Pr

        qm = self.states.lng.mas_flow
        An = self.params.An #površina nastrujavanja (poprečni presjek), m2
        d_ekv = self.params.d_ekv #ekvivalentni promjer, m
        molar_mas = self.states.lng.M #kg/kmol
        acf = self.states.lng.acf #akomodacijski faktor
        acef =  self.states.lng.acef #akomodacijski faktor energije

        w = qm / ro / An #brzina strujanja, m/s
        Re = w * d_ekv / nu
        self.states.lng.Re = Re
        #print(Re)

        R = 8314 / molar_mas #J/kgK
        cv = cp-R #J/kgK
        kappa = cp/cv
        #print(kappa)
        ws = (kappa*R*Tm)**0.5 #brzina zvuka, m/s
        #print(ws)
        self.states.lng.ws = ws
        Ma = w/ws
        #print(Ma)
        self.states.lng.Ma = Ma
        Kn = (2*kappa / np.pi) **0.5 * Ma / Re
        #print(Kn)
        self.states.lng.Kn = Kn
        kkl = 4 *Kn * (2-acf) / acf #koeficijent klizanja
        #print(kkl)
        self.states.lng.kkl = kkl
        kts =8 * (2-acef)/acef * kappa/(kappa+1) * Kn/Pr #koeficijent temp. skoka
        #print(kts)
        self.states.lng.kts = kts
        Nu = 1.0 / (kts/4.0 + (17.0 + 84.0*kkl + 105.0 * kkl**2)/(140.0 * (1.0 + 3.0*kkl)**2 ) )
        #print(Nu)
        self.states.lng.Nu = Nu
        alfa = Nu * lmbd / d_ekv
        #print(alfa)
        self.states.lng.alfa = alfa

    def calc_T_out_glyc(self):
        T_in =self.states.glycol.T_in
        qm = self.states.glycol.mas_flow #kg/s
        Fi = self.states.common.heat_flow #W
        glyc.set_T(T_in)
        cp0 = glyc.get_cp()
        #print(qm)
        #print(cp0)
        T_out0 = T_in - Fi / (qm*cp0)
        #print(T_out0)
        T_m_prv = (T_in+T_out0)/2
        T_m = T_m_prv-1.0
        while (abs(T_m-T_m_prv) > 0.00001):
            glyc.set_T(T_m)
            cp = glyc.get_cp()
            T_out = T_in - Fi / (qm*cp)
            #print(T_out)
            T_m_prv = T_m
            T_m = (T_in+T_out)/2
        self.states.glycol.T_out = T_out
        self.states.glycol.Tm = T_m
        glyc.set_T(T_m)
        self.states.glycol.c = glyc.get_cp() #J/kgK
        self.states.glycol.l = glyc.get_lambda() #W/mK
        self.states.glycol.vis_dyn = glyc.get_dynamic_viscosity() #Pas
        self.states.glycol.ro = glyc.get_density() #kg/m3
        self.states.glycol.vis_cin = glyc.get_kinematic_viscosity() #m2/s



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
             ("pi3>1!!!")
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
