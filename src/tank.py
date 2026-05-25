import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
#import src.properties.liq_props as liq
import src.properties.vap_props as vap
#import matplotlib.pyplot as plt
from src.heat_transfer_interface import interf_ht_fun, liq_alpha_interf, vap_alpha_interf
from src.properties.liq_props_lowlev import MetLiq
from src.properties.vap_props_lowlev import MetVap, MetVap_slow
from src.properties.equilibrium_props_lowlev import MetEq

M=16.04
met_liq = MetLiq()
met_vap = MetVap()
met_vap_slow = MetVap_slow()
met_eq = MetEq()
class Tank_params():
    '''
    Parametri spremnika
    '''
    def __init__(self, int_length, int_diameter, thickness, k_liq, k_vap, T_amb, slosh_times=[]):
        self.int_length = int_length #unutarnja duljina, m
        self.int_diameter = int_diameter #unutarnji promjer , m
        self.thickness = thickness #debljina stijenke, m
        self.k_liq = k_liq #koeficijent prolaza topline na strani kapljevitog LNG, W/m2K
        self.k_vap = k_vap #koeficijent prolaza topline na strani parovitog LNG, W/m2K
        self.T_amb = T_amb #ambient temperature, K
        self.int_volume = int_diameter ** 2. * 3.1415 / 4. * int_length
        self.surf_area_tot = 2.0 * (int_diameter+thickness) ** 2 * 3.1415 / 4.0  + (int_diameter+thickness) * 3.1415 * int_length
        self.slosh_times = slosh_times
class Tank_initial():
    '''
    Pocetna stanja spremnika
    '''
    def __init__(self, pressure, temperature_l, temperature_v, liq_vol_ratio):
        self.pressure = pressure #tlak LNG, bar
        self.temperature_l = temperature_l #temperatura LNG, K
        self.temperature_v = temperature_v #temperatura LNG, K
        self.liq_vol_ratio = liq_vol_ratio # volumni udio kapljevine u spremniku, -

class Tank_states():
    '''
    Stanja spremnika
    '''
    def __init__(self, tank_initial, tank_params):
        self.liquid = Tank_liq_states(tank_initial, tank_params) # stanja kapljevine
        self.vapor = Tank_vap_states(tank_initial, tank_params) # stanja pare
        self.common = Tank_common_states(tank_initial)

class Tank_liq_states():
    '''
    Stanja kapljevine u spremniku
    '''
    def __init__(self, tank_initial, tank_params):
        self.temperature = tank_initial.temperature_l #K
        met_liq.set_T(self.temperature)
        self.vol = tank_params.int_volume * tank_initial.liq_vol_ratio #m3
        self.vol_prv = self.vol
        self.density = met_liq.get_density() #kg/m^3
        self.quantity = self.vol * self.density / M #kmol
        self.molar_enthalpy =met_liq.get_mol_enthalpy() #kJ/kmol
        self.enthalpy = self.molar_enthalpy * self.quantity #kJ
        self.vol_ratio = tank_initial.liq_vol_ratio
        self.surf_area = 0.0
        self.surf_heat_flow = 0.0 #
        self.flow = 0.0 #kmol/s
        self.flow_quantity = 0.0 #kmol
        self.alpha_int = 0.0 #W/m2K
        self.flow_ent = 0.0
        self.Nu = 0.0
        self.Pr = 0.0
        self.Ra = 0.0

        self.slosh_f = 0.0
        self.Re_s = 0.0 #sloshing Reynolds
        self.Nu_s = 0.0 #sloshing Nusselt
        self.h_sat = 0.0

class Tank_vap_states():
    '''
    Stanja pare u spremniku
    '''
    def __init__(self, tank_initial, tank_params):
        self.temperature = tank_initial.temperature_v #K
        self.vol = tank_params.int_volume * (1.0 - tank_initial.liq_vol_ratio) #m3
        T = self.temperature
        p = tank_initial.pressure
        met_vap.set_pT(p, T)
        density = met_vap.get_density() #kg/m^3 #
        self.density = density #kg/m^3 #
        self.quantity = self.vol * self.density / M #kmol
        self.molar_enthalpy =met_vap.get_mol_enthalpy() #kJ/kmol
        self.molar_inten = met_vap.get_mol_inten()
        self.inten = self.molar_inten * self.quantity #kJ
        self.surf_area = 0.0 #površina prolaz topline - para - zrak, m2
        self.molar_density = density / M * 1000.0
        self.surf_heat_flow = 0.0
        self.flow = 0.0
        self.flow_quantity = 0.0 #kmol
        self.alpha_int = 0.0 #W/m2K
        self.flow_ent = 0.0
        self.beta = 0.0
        self.Nu = 0.0
        self.Pr = 0.0
        self.Ra = 0.0
        self.h_sat = 0.0

class Tank_common_states():
    def __init__(self, tank_initial):
        self.pressure = tank_initial.pressure #bar
        met_liq.set_p(self.pressure)
        self.T_sat = met_liq.get_T()
        self.angle = 0.0 #kut kružnog odsječka, rad
        self.segment = 0.0 # duljina kružnog odsječka, m
        self.arc = 0.0 #kružni luk, m
        self.liq_height = 0.0 #visina kapljevine u spremniku, m
        self.interface_area = 0.0 #kut dodira kapljevine i pare, m2
        self.interface_perimiter = 0.0 #opseg oplakane povrsine
        self.heat_interface = 0.0
        self.evaporation = 0.0
        self.V_interface=0.0
        self.alpha_liq_int = 0.0
        self.alpha_vap_int = 0.0
        self.T_int = 0.0
        self.k_int = 0.0
        self.interface_resistance = 0.0
        self.evaporation_heat = 0.0 #kJ - latentna toplina isparivanja
        self.evap_flow = 0.0
        self.slosh_state = False
        self.slosh_prv = False
        self.time_layer = 0.0 #s
        self.qm_slosh = 0.0
        self.Fi_slosh = 0.0
        self.hlv = 0.0
        self.delta_layer = 0.0
        self.W = 0.0

class TankSave:
    '''
    Saving tank states
    '''
    def __init__(self, states):
        '''
        https://stackoverflow.com/questions/11637293/iterate-over-object-attributes-in-python
        '''
        #atributi u listu
        self.liquid_st_list = [a for a in dir(states.liquid) if not a.startswith('__') and not callable(getattr(states.liquid, a))]
        #atributi u listu
        self.vapor_st_list = [a for a in dir(states.vapor) if not a.startswith('__') and not callable(getattr(states.vapor, a))]
        #atributi u listu
        self.common_st_list = [a for a in dir(states.common) if not a.startswith('__') and not callable(getattr(states.common, a))]

        self.liquid = {i : [] for i in self.liquid_st_list}
        self.vapor = {i : [] for i in self.vapor_st_list}
        self.common = {i : [] for i in self.common_st_list}

    def save_states(self,states):
        #iteriraj po svim atributima
        for liq_st in self.liquid_st_list:
            #trenutna vrijednost atributa
            current_state = getattr(states.liquid, liq_st)
            # spremi u dictionary - dodaj zadnju vrijednost liste
            self.liquid[liq_st].append(current_state)
        for vap_st in self.vapor_st_list:
            current_state = getattr(states.vapor, vap_st)
            self.vapor[vap_st].append(current_state)
        for com_st in self.common_st_list:
            current_state = getattr(states.common, com_st)
            self.common[com_st].append(current_state)

class Tank():
    '''
    Glavna klasa - spremnik
    '''
    def __init__(self, tank_params, tank_initial):
        self.params = tank_params #parametri
        self.initial = tank_initial #pocetna stanja
        self.states = Tank_states(tank_initial, tank_params) #stanja proracuna
        self.met_liq = MetLiq()
        self.tank_areas() #površine izmjene topline i dodira faza
        self.save = TankSave(self.states) #inicijalizacija spremanja rezultata
        self.save.save_states(self.states) #spremanje pocetnih stanja

    def tank_areas(self):
        d = self.params.int_diameter
        length = self.params.int_length
        vol_ratio_l = self.states.liquid.vol_ratio
        thickness = self.params.thickness
        area_total = self.params.surf_area_tot
        vol_ratio_v = 1 - vol_ratio_l
        area_vap = vol_ratio_v * d**2 * 3.1415 / 4.0
        radius = d/2
        rad_sq = radius**2
        const = 2 * area_vap / rad_sq
        rez = opt.minimize_scalar(angle_func, bounds=(0, 2*np.pi), args=(const), method='bounded')
        angle = rez.x
        self.states.common.angle = angle
        segment = 2.0 * np.sin(angle/2.0) * radius
        height = np.cos(angle/2.0) * radius + radius
        self.states.common.liq_height =height
        interf_area = segment * length
        self.states.common.segment = segment
        self.states.common.interface_area = interf_area
        self.states.common.interface_perimiter =2*(segment + length)
        arc = angle*(radius+thickness/2.0) #kružni luk
        self.states.common.arc = arc
        area_cylinder = arc*(length+thickness)
        area_sides = 2. * (d+thickness)**2 * 3.1415 / 4.0
        area_vapor = area_cylinder + area_sides
        self.states.vapor.surf_area = area_vapor
        self.states.liquid.surf_area = area_total - area_vapor


    def interface_heat_flow(self, dt):

        Tl = self.states.liquid.temperature
        Tv = self.states.vapor.temperature
        p=self.states.common.pressure
        area = self.states.common.interface_area
        perimeter = self.states.common.interface_perimiter
        L_ekv =area/perimeter

        if (Tv-0.0001)>(Tl+0.0001):
            rez = opt.minimize_scalar(interf_ht_fun, bounds=(Tl+0.0001, Tv-0.0001), args=(Tl, Tv, L_ekv, p, area), method='bounded')
            T_interface =rez.x
            alpha_liq, Nu_liq, Pr_liq, Ra_liq = liq_alpha_interf(Tl, T_interface, L_ekv)
            alpha_vap, Nu_vap, Pr_vap, Ra_vap = vap_alpha_interf(Tv, T_interface, p, L_ekv)
            if Pr_vap > 30.0:
                    print('Tm: ' + str((Tv + T_interface)/2) + ' p: ' + str(p) )
            resistance = (1./alpha_liq + 1./alpha_vap) / area
            heat = (Tv - Tl) / resistance * dt / 1000.0 #kJ
            k_int = 1/(1/alpha_liq+1/alpha_vap)
        else:
            heat = 0.0
            T_interface =(Tv+Tl)*0.5
            alpha_liq = 0.0
            alpha_vap = 0.0
            resistance  = 1000.0
            k_int = 0.0
            Nu_liq = 0.0
            Pr_liq = 0.0
            Ra_liq = 0.0
            Nu_vap = 0.0
            Pr_vap = 0.0
            Ra_vap = 0.0
        #print('iteracija interface = ', iteracija)
        self.states.common.heat_interface = heat
        self.states.common.T_int = T_interface
        self.states.liquid.alpha_int = alpha_liq
        self.states.vapor.alpha_int = alpha_vap
        self.states.common.interface_resistance = resistance
        self.states.common.k_int = k_int
        self.states.liquid.Nu = Nu_liq
        self.states.liquid.Pr = Pr_liq
        self.states.liquid.Ra = Ra_liq
        self.states.vapor.Nu = Nu_vap
        self.states.vapor.Pr = Pr_vap
        self.states.vapor.Ra = Ra_vap

    def interface_heat_flow_dbg(self, dt):

        Tl = self.states.liquid.temperature
        Tv = self.states.vapor.temperature
        p=self.states.common.pressure
        area = self.states.common.interface_area
        perimeter = self.states.common.interface_perimiter
        L_ekv =area/perimeter
        rez = opt.minimize_scalar(interf_ht_fun, bounds=(Tl+0.05, Tv-0.05), args=(Tl, Tv, L_ekv, p, area), method='bounded')
        T_interface =(Tl+Tv)/2
        alpha_liq = 1000.0
        alpha_vap = 30.0
        resistance = (1./alpha_liq + 1./alpha_vap) / area
        k = 3.0
        heat = (Tv - Tl) *k  * area * dt / 1000.0 #kJ
        #print('iteracija interface = ', iteracija)
        self.states.common.heat_interface = heat
        self.states.common.T_int = T_interface
        self.states.liquid.alpha_int = alpha_liq
        self.states.vapor.alpha_int = alpha_vap
        self.states.common.interface_resistance = 1/(k*area)
        self.states.common.k_int = k

    def slosh(self):
        kolicina_liq = self.states.liquid.quantity * 1000. #mol
        kolicina_vap = self.states.vapor.quantity * 1000. #mol
        entalpija_liq = self.states.liquid.enthalpy * 1000. # J
        entalpija_vap = self.states.vapor.enthalpy * 1000. # J

        kolicina_total = kolicina_liq + kolicina_vap # mol
        entalpija_total = entalpija_liq + entalpija_vap # J
        V_total = self.params.int_volume # m3

        mol_H = entalpija_total / kolicina_total #J/mol
        mol_V = V_total / kolicina_total # m3/mol
        met_eq.set_hV(mol_H, mol_V)
        T = met_eq.get_T()
        T_liq =self.states.liquid.temperature
        T_vap =self.states.vapor.temperature
        p_before = self.states.common.pressure

        #print('Tliq = ' +str(T_liq) +' K')
        #print('Tvap = ' +str(T_vap) +' K')
        #print('T = ' +str(T) +' K')
        #print('')

        p = met_eq.get_p()
        #print('p before = ' +str(p_before/100000) +' bar')
        #print('p = ' +str(p/100000) +' bar')

    def slosh_correction(self):
        L = self.params.int_length
        h = self.states.common.liq_height
        g = 9.81
        pi = np.pi
        f = ((g*pi/L)*np.tanh(h*pi/L))**0.5
        self.states.liquid.slosh_f = f
        Tliq = self.states.liquid.temperature
        met_liq.set_T(Tliq)
        nu=met_liq.get_kinematic_viscosity()
        Re_s = f*L*L/nu
        self.states.liquid.Re_s = Re_s
        Pr = met_liq.get_Pr()
        #print(Pr)
        Nu_s = (Re_s/4000.0)**0.69 *Pr**(1/3)
        self.states.liquid.Nu_s = Nu_s


    def slosh_semi_infinite(self):
        self.slosh_correction()
        Nu_s = self.states.liquid.Nu_s
        Tl = self.states.liquid.temperature
        Tv = self.states.vapor.temperature
        Ts = self.states.common.T_sat
        p=self.states.common.pressure
        Tm = (Tl + Ts) / 2.0
        met_liq.set_T(Tm)
        lmbd0 = met_liq.get_lambda()
        lmbd_ef = lmbd0 * Nu_s  #  rezonantnoi sloshing
        #t = self.states.common.time_layer
        t = 60*60*5
        a = met_liq.get_thermal_diff()
        area = self.states.common.interface_area
        delta_layer = 3.0*(a*t)**0.5
        delta_layer = min(delta_layer, 0.1)
        Fi_slosh = area*lmbd_ef*(Ts-Tl)/delta_layer / 1000.0 #kW
        met_liq.set_T(Ts)
        met_vap.set_pT_robust(p, Tv)
        hl = met_liq.get_mol_enthalpy() #kJ/kmol
        hv = met_vap.get_mol_enthalpy() #kJ/kmol
        hlv = hv-hl #kJ/kmol
        qm_slosh =Fi_slosh/hlv #kmol/s
        self.states.common.delta_layer=delta_layer
        self.states.liquid.h_sat = hl
        self.states.vapor.h_sat = hv
        self.states.common.hlv = hlv
        self.states.common.qm_slosh = qm_slosh
        self.states.common.Fi_slosh = Fi_slosh
        print('Ts-Tl slosh:' +str(Ts-Tl))

    def update_states(self, liq_flow, vap_flow, T_liq_flow, T_vap_flow, dt):
        current_time = int(self.states.common.time_layer)
        if current_time in self.params.slosh_times:
            slosh_state = True
        else:
            slosh_state = False
        self.states.common.time_layer = dt+self.states.common.time_layer
        self.states.liquid.vol_prv = self.states.liquid.vol
        self.states.liquid.flow = liq_flow #kmol/s
        self.states.vapor.flow = vap_flow #kmol/s
        liq_N = self.states.liquid.quantity #kolicina kapljevine prije koraka, kmol
        vap_N = self.states.vapor.quantity #kolicina pare prije koraka, kmol
        liq_dens = self.states.liquid.density #gustoca kapljevine prije koraka, kg/m^3
        vap_dens = self.states.vapor.density
        T_sat = self.states.common.T_sat # K
        #print("T_sat "+ str(T_sat))
        #print("")
        T_liq = self.states.liquid.temperature # temperatura kapljevine prije koraka, K
        T_vap = self.states.vapor.temperature #temperatura pare prije koraka, K
        A_interface = self.states.common.interface_area # površina dodira kapljevine i pare ( kružni odsječak)
        V_interface = 0.005 * A_interface #volumen dodira kapljevine i pare, vidi Analysis of LNG transport losses
        self.states.common.V_interface = V_interface
        #print("V_interface "+ str(V_interface))
        q_liq_ev = liq_dens / M * 0.1 *(T_liq - T_sat) / T_sat #kmol/m^3/s
        q_liq_ev = max(0.0, q_liq_ev)
        q_vap_cond = vap_dens / M * 0.1 *(T_sat - T_vap) / T_sat  #kmol/m^3/s
        q_vap_cond= max(0.0, q_vap_cond)
        q_ev = q_liq_ev - q_vap_cond
        #print(q_ev)
        liq_evap = q_ev * V_interface * dt # količina ishlapljene kapljevine, kmol
        #print("")
        #print("liq_evap "+ str(liq_evap))
        if T_liq_flow == -1: #odvodi se kapljevina stanja spremnika
            liq_flow_mol_ent = self.states.liquid.molar_enthalpy
        else: #dovodi se kapljevina drugog stanja
            met_liq.set_T(T_liq_flow)
            liq_flow_mol_ent = met_liq.get_mol_enthalpy()
        if T_vap_flow == -1: #odvodi se para stanja spremnika
            vap_flow_mol_ent = self.states.vapor.molar_enthalpy
        else: #dovodi se para drugog stanja
            p_vap = self.states.common.pressure
            if T_vap_flow < 190.0:
                met_vap.set_T_sat(T_vap_flow) # dolazi szp s PBU -
            else:
                met_vap.set_pT(p_vap, T_vap_flow)
            vap_flow_mol_ent = met_vap.get_mol_enthalpy()
        met_vap.set_T_sat(T_sat)
        sat_vap_mol_ent =  met_vap.get_mol_enthalpy() #specifična entalpija szp kod isparavanja, kJ/kmol
        liq_flow_quantity = liq_flow * dt #kmol
        vap_flow_quantity = vap_flow * dt #kmol
        self.states.liquid.flow_quantity = liq_flow_quantity
        self.states.vapor.flow_quantity = vap_flow_quantity
         #entalpija kapljevine
        liq_surf_area = self.states.liquid.surf_area # povrsina oplakana kapljevinom, m^2
        vap_surf_area =  self.states.vapor.surf_area #  povrsina oplakana parom, m^2
        heat_gain_liq =liq_surf_area * self.params.k_liq * (self.params.T_amb -  T_liq) * dt / 1000.0 #toplina izvana prema kapljevini, kJ
        #print("")
        #print("heat_gain_liq "+ str(heat_gain_liq))
        heat_gain_vap = vap_surf_area * self.params.k_vap * (self.params.T_amb -  T_vap) * dt / 1000.0 #toplina izvana prema pari, kJ
        #print("")
        #print("heat_gain_vap "+ str(heat_gain_vap))
        self.interface_heat_flow(dt)
        if slosh_state:
            print('time slosh: ' +str(current_time))
            self.slosh_semi_infinite()
            #heat_interface =0.0
            #liq_evap = 0.0
        else:
            self.states.liquid.slosh_f = 0.0
            self.states.liquid.Re_s = 0.0
            self.states.liquid.Nu_s = 1.0
            self.states.common.qm_slosh = 0.0
            self.states.common.Fi_slosh = 0.0
            self.states.common.delta_layer = 0.0
            met_liq.set_T(T_sat)
            met_vap.set_T_sat(T_sat)
            hl = met_liq.get_mol_enthalpy()
            hv = met_vap.get_mol_enthalpy()
            hlv = hv-hl
            self.states.liquid.h_sat = hl
            self.states.vapor.h_sat = hv
            self.states.common.hlv = hlv

        heat_interface = self.states.common.heat_interface


        #print("")
        #print("heat_interface "+ str(heat_interface))
        # entalpija kapljevine
        liq_flow_heat = liq_flow*liq_flow_mol_ent * dt
        vap_flow_heat = vap_flow*vap_flow_mol_ent * dt
        evaporation_heat = liq_evap*sat_vap_mol_ent
        Q_slosh = self.states.common.Fi_slosh * dt
        m_slosh = self.states.common.qm_slosh * dt
        hl_sat = self.states.liquid.h_sat
        hv_sat = self.states.vapor.h_sat
        liq_quantity = liq_N + liq_flow_quantity - liq_evap + m_slosh #količina kapljevine nakon isparavanja i oduzimanja, kmol
        vap_quantity = vap_N + vap_flow_quantity + liq_evap - m_slosh #količina pare nakon isparavanja i oduzimanja, kmol

        p = self.states.common.pressure
        liq_vol_prv = self.states.liquid.vol_prv
        W_prv = -100.0
        W=0.0
        while abs(W-W_prv) > 0.001:
            liq_ent = self.states.liquid.enthalpy + W + liq_flow_heat - evaporation_heat + heat_gain_liq + heat_interface + Q_slosh + m_slosh*hl_sat
            liq_mol_ent = liq_ent / liq_quantity #kJ/kmol
            met_liq.set_molH(liq_mol_ent)
            liq_temperature = met_liq.get_T()
            liq_density = met_liq.get_density() #kg/m^3
            liq_volume = liq_quantity * M / liq_density # m^3
            W_prv = W
            W = (liq_vol_prv-liq_volume)*p/1000.0





        #entalpija pare
        vap_inten = self.states.vapor.inten - W + vap_flow_heat + evaporation_heat + heat_gain_vap - heat_interface - m_slosh*hv_sat
        vap_mol_inten = vap_inten / vap_quantity #kJ/kmol

        vapor_volume = self.params.int_volume - liq_volume
        vap_mol_dens =  (vap_quantity * 1000.0) / vapor_volume #mol/m3

        #liq volime - vap volume - vap density + vap molar ehthalpy - pressure i temperatura
        #if slosh_state:
            #print([vap_inten, vap_quantity])
            #print([vap_mol_inten, vap_mol_dens])
        met_vap_slow.set_UD_molar(vap_mol_inten, vap_mol_dens)
        vap_phase = met_vap_slow.get_phase()
        vap_temperature, pressure = met_vap_slow.get_Tp()
        if (vap_phase == 5) or (vap_phase == 2) : #gas
            x_vap = 1.0
            liq_cond = 0.0
            vap_mol_ent = met_vap_slow.get_mol_enthalpy()
        elif vap_phase == 6: #2 phases
            #print('mokra paraaa!!!!!!!!!!!!')
            ######  RAD KOD KONDENZACIJE??
            #Wcond_prv = -100.0
            #Wcond=0.0
            #while abs(Wcond-Wcond_prv) > 0.01:
            x_vap = met_vap_slow.get_x() #racunaj udio pare
            liq_cond = (1.0-x_vap) * vap_quantity #kondenzirana kapljevina, kmol
            #print(x_vap)
            liq_quantity = liq_quantity + liq_cond #ukupna kol. kapljevine s kondenzatom, kmol
            met_liq.set_p(pressure) #stanje kondenzata
            cond_enthalpy = met_liq.get_mol_enthalpy() #mol entalpija kondentara, kJ/kmol
            V_liq_before_cond = liq_volume


            liq_ent = liq_ent + liq_cond*cond_enthalpy # + Wcond #ukupna entalpija kapljevine s kondenzatom, kJ
            liq_mol_ent = liq_ent / liq_quantity  #kJ/kmol
            met_liq.set_molH(liq_mol_ent)
            liq_temperature = met_liq.get_T()
            liq_density = met_liq.get_density() #kg/m^3
            liq_volume = liq_quantity * M / liq_density # m^3
            #Wcond_prv = Wcond
            #Wcond = (liq_volume-V_liq_before_cond)*pressure/1000.0

            vap_quantity = x_vap*vap_quantity #ostatak je para, kmol
            met_vap.set_p_sat(pressure) #stanje pare - zasicena para
            vap_mol_ent = met_vap.get_mol_enthalpy() #molarna entalpija pare, kJ/kmol
            vap_mol_inten = met_vap.get_mol_inten()
            vap_inten = vap_mol_inten*vap_quantity #ukupna entalpija pare, kJ
            vapor_volume = self.params.int_volume - liq_volume
            vap_mol_dens =  (vap_quantity * 1000.0) / vapor_volume #mol/m3
        else:
            print("WTF??")
            print(vap_phase)

        #print(pressure)
        met_liq.set_p(pressure)
        T_sat = met_liq.get_T()
        self.states.common.evap_flow = liq_evap
        self.states.liquid.vol_ratio = liq_volume / self.params.int_volume
        self.states.vapor.vol = vapor_volume
        self.states.liquid.vol = liq_volume
        self.states.liquid.temperature = liq_temperature
        self.states.liquid.density = liq_density
        self.states.liquid.quantity = liq_quantity
        self.states.vapor.quantity = vap_quantity
        self.states.vapor.molar_density = vap_mol_dens
        self.states.vapor.molar_inten = vap_mol_inten
        self.states.vapor.molar_enthalpy = vap_mol_ent

        self.states.liquid.enthalpy =liq_ent
        self.states.vapor.inten = vap_inten
        self.states.liquid.molar_enthalpy = liq_mol_ent
        self.states.vapor.temperature = vap_temperature
        self.states.common.pressure = pressure
        self.states.common.T_sat = T_sat
        self.states.common.evaporation = liq_evap
        self.states.common.heat_interface = heat_interface
        self.states.liquid.surf_heat_flow = heat_gain_liq #kJ
        self.states.vapor.surf_heat_flow = heat_gain_vap #KJ
        self.states.liquid.flow_ent = liq_flow_heat #KJ - entalpija struje
        self.states.vapor.flow_ent = vap_flow_heat #KJ - entalpija struje
        self.states.common.evaporation_heat = evaporation_heat
        self.tank_areas()
        self.states.common.W = W
        self.save.save_states(self.states) # spremi stanja

def angle_func(alpha, const):
    er = alpha - np.sin(alpha) - const
    error = er**2
    return error

#%%
if __name__ == "__main__":
#     #primjer koristenja:



    int_length = 14.5 #m
    int_diameter = 5.9 #m
    thickness = 0.05 #m
    k_liq = (1+0.4)*1/(0.2/0.0225+1/2000+1/20) #W/mK #excel
    k_vap = (1+0.4)*1/(0.3/0.015+1/200+1/20) #W/mK #excek
    pressure = 283000. #Pa
    temperature_l = 125.7 #K
    temperature_v = 125.9 #K

    liq_vol_ratio = 0.75 #
    T_amb = 293.15 #K
    tank_params = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb)
    tank_initial = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
    tank = Tank(tank_params,tank_initial)
    vap_flow = 0.0
    for i in range(450):
        # if tank.states.common.pressure > 700000.0:
        #     vap_flow = -0.01
        # elif tank.states.common.pressure > 600000.0:
        #     vap_flow = -0.005
        # elif tank.states.common.pressure < 300000.0:
        #      vap_flow = 0.0
        if i != 200 :
            tank.update_states(0.0, vap_flow, -1, -1, 10)
        else:
            print(i)
            tank.update_states(0.0, vap_flow, -1, -1, 10, True)

        if i%1000 == 0:
            print("               " + str(i) + "                    ")

#%%

    times =np.array(range(450)) *10 / 60 / 60
    pressures = np.array(tank.save.common['pressure']) / 1000.0
    pressures=pressures[1:]
    plt.plot(times, pressures)
    plt.plot(tank.save.common['pressure'])
    plt.figure()
    plt.plot(tank.save.liquid['temperature'])
    plt.figure()
    plt.plot(tank.save.vapor['temperature'])
    plt.figure()
    plt.plot(tank.save.common['T_sat'])
    plt.figure()
    plt.plot(tank.save.vapor['flow'])
    plt.figure()
    plt.plot(tank.save.liquid['vol_ratio'])
    plt.figure()
    plt.plot(tank.save.common['evaporation'])

    plt.plot(tank.save.common['qm_slosh'])

    plt.plot(tank.save.vapor['enthalpy'])

    plt.figure()
    plt.plot(tank.save.vapor['quantity'])

    fig, ax1 = plt.subplots()
    ax1.plot(tank.save.common['pressure'])
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(tank.save.common['evaporation'],  color='g')


    fig, ax1 = plt.subplots()
    ax1.plot(tank.save.vapor['temperature'])
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(tank.save.common['evaporation'], color='r')

    fig, ax1 = plt.subplots()
    ax1.plot(tank.save.liquid['temperature'])
    ax1.plot(tank.save.common['T_sat'], color='r')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(tank.save.common['evaporation'], color='g')

    plt.figure()
    plt.plot(tank.save.vapor['molar_density'])
    tank.save.vapor['molar_enthalpy'][-1]

    plt.plot(tank.save.liquid['density'])

    plt.figure()
    plt.plot(tank.save.vapor['vol'])
    hf = np.array(tank.save.liquid['surf_heat_flow'])*1000/10
    hfv = np.array(tank.save.vapor['surf_heat_flow'])*1000/10
    hf_tot = hf + hfv
    plt.plot(hf_tot)


    dens = tank.states.vapor.molar_density
    #liq volime - vap volume - vap density + vap molar ehthalpy - pressure i temperatura
    ent = tank.states.vapor.molar_enthalpy
    vol = tank.states.vapor.vol


    plt.plot(tank.save.common['k_int'])




    fig, ax1 = plt.subplots()
    ax1.plot(tank.save.common['evaporation_heat'], color='r')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(tank.save.common['evaporation'], color='g')




#%%
if __name__ == "__main__":
    int_length = 0.909 #m
    int_diameter = 0.6 #m
    thickness = 0.05 #m
    k_liq = 1/(0.3/0.015+1/2000+1/20) #W/mK #excel
    k_vap = 1/(0.3/0.015+1/200+1/20) #W/mK #excek
    pressure = 283000. #Pa
    temperature_l = 125.7 #K
    temperature_v = 125.8 #K

    liq_vol_ratio = 0.75 #
    T_amb = 293.15 #K
    tank_params = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb)
    tank_initial = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
    tank = Tank(tank_params,tank_initial)
    vap_flow = 0.0
    for i in range(45000):
        # if tank.states.common.pressure > 700000.0:
        #     vap_flow = -0.01
        # elif tank.states.common.pressure > 600000.0:
        #     vap_flow = -0.005
        # elif tank.states.common.pressure < 300000.0:
        #      vap_flow = 0.0
        tank.update_states(0.0, vap_flow, -1, -1, 10)
        if i%1000 == 0:
            print("               " + str(i) + "                    ")

    times =np.array(range(45000)) *10 / 60 / 60
    pressures = np.array(tank.save.common['pressure']) / 1000.0
    pressures=pressures[1:]
    plt.plot(times, pressures)

    plt.figure()
    plt.plot(tank.save.vapor['Pr'])
    plt.plot(tank.save.vapor['Ra'])

#%%
    plt.plot(tank.save.liquid['temperature'], color='b')
    plt.plot(tank.save.vapor['temperature'], color='r')
    plt.plot(tank.save.common['T_int'], color='g')

    fig, ax1 = plt.subplots()
    ax1.plot(tank.save.liquid['alpha_int'], color='b')
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.plot(tank.save.vapor['alpha_int'], color='r')

    plt.plot(tank.save.common['k_int'], color='g')
    plt.plot(tank.save.common['T_int'], color='g')

    plt.plot(tank.save.liquid['alpha_int'], color='b')
    plt.plot(tank.save.vapor['alpha_int'], color='r')


#%% #trajanje


    int_length = 14.5 #m
    int_diameter = 5.9 #m
    thickness = 0.05 #m
    k_liq = (1+0.4)*1/(0.19/0.023+1/2000+1/20) #W/mK #excel
    k_vap = (1+0.4)*1/(0.19/0.023+1/200+1/20) #W/mK #excek
    pressure = 100000. #Pa
    temperature_l = 112.9 #K
    temperature_v = 112.91 #K

    liq_vol_ratio = 0.83 #
    T_amb = 293.15 #K
    tank_params = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb)
    tank_initial = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
    tank = Tank(tank_params,tank_initial)
    vap_flow = 0.0
    times=0
    while (tank.states.common.pressure < 810000):
        # if tank.states.common.pressure > 700000.0:
        #     vap_flow = -0.01
        # elif tank.states.common.pressure > 600000.0:
        #     vap_flow = -0.005
        # elif tank.states.common.pressure < 300000.0:

        tank.update_states(0.0, vap_flow, -1, -1, 10)
        times = times+10
        #print(tank.states.common.pressure)
        if times%10000 == 0:
            print("               " + str(times) + "          " + str(tank.states.common.pressure))
    timed = times/60/60/24
    print(timed)
    plt.plot(tank.save.common["pressure"])
