import numpy as np
import scipy.optimize as opt
import src.properties.liq_props as liq
import src.properties.vap_props as vap
import matplotlib.pyplot as plt
from src.heat_transfer_interface import interf_ht_fun, liq_alpha_interf, vap_alpha_interf
class Tank_params():
    '''
    Parametri spremnika
    '''
    def __init__(self, int_length, int_diameter, thickness, k_liq, k_vap, k_int, T_amb):
        self.int_length = int_length #unutarnja duljina, m
        self.int_diameter = int_diameter #unutarnji promjer , m
        self.thickness = thickness #debljina stijenke, m
        self.k_liq = k_liq #koeficijent prolaza topline na strani kapljevitog LNG, W/m2K
        self.k_vap = k_vap #koeficijent prolaza topline na strani parovitog LNG, W/m2K
        self.k_int = k_int #  #koeficijent prolaza topline na dodiru pare i kapljevine, W/m2K (interface)
        self.T_amb = T_amb #ambient temperature, K
        self.int_volume = int_diameter ** 2. * 3.1415 / 4. * int_length
        self.surf_area_tot = 2.0 * (int_diameter+thickness) ** 2 * 3.1415 / 4.0  + (int_diameter+thickness) * 3.1415 * int_length

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
        self.vol = tank_params.int_volume * tank_initial.liq_vol_ratio #m3
        self.density = liq.density(self.temperature) #kg/m^3
        self.quantity = self.vol * self.density / liq.M #kmol
        self.molar_enthalpy =liq.enthalpy_sat(self.temperature) #kJ/kmol
        self.enthalpy = self.molar_enthalpy * self.quantity #kJ
        self.vol_ratio = tank_initial.liq_vol_ratio
        self.surf_area = 0.0
        self.surf_heat_flow = 0.0 #
        self.flow = 0.0 #kmol/s
        self.flow_quantity = 0.0 #kmol
        self.alpha_int = 0.0 #W/m2K
        self.flow_ent = 0.0

class Tank_vap_states():
    '''
    Stanja kapljevine u spremniku
    '''
    def __init__(self, tank_initial, tank_params):
        self.temperature = tank_initial.temperature_v #K
        self.vol = tank_params.int_volume * (1.0 - tank_initial.liq_vol_ratio) #m3
        T = self.temperature
        p = tank_initial.pressure
        density = vap.density(T, p) #kg/m^3 #
        self.density = density #kg/m^3 #
        self.quantity = self.vol * self.density / liq.M #kmol
        self.molar_enthalpy =vap.enthalpy(T, p) #kJ/kmol
        self.enthalpy = self.molar_enthalpy * self.quantity #kJ
        self.surf_area = 0.0 #površina prolaz topline - para - zrak, m2
        self.molar_density = density / vap.M * 1000.0
        self.surf_heat_flow = 0.0
        self.flow = 0.0
        self.flow_quantity = 0.0 #kmol
        self.alpha_int = 0.0 #W/m2K
        self.flow_ent = 0.0

class Tank_common_states():
    def __init__(self, tank_initial):
        self.pressure = tank_initial.pressure #bar
        self.T_sat = liq.saturation_temperature(self.pressure)
        self.angle = 0.0 #kut kružnog odsječka, rad
        self.segment = 0.0 # duljina kružnog odsječka, m
        self.arc = 0.0 #kružni luk, m
        self.interface_area = 0.0 #kut dodira kapljevine i pare, m2
        self.interface_perimiter = 0.0 #opseg oplakane povrsine
        self.heat_interface = 0.0
        self.evaporation = 0.0
        self.V_interface=0.0
        self.alpha_liq_int = 0.0
        self.alpha_vap_int = 0.0
        self.T_int = 0.0
        self.interface_resistance = 0.0
        self.evaporation_heat = 0.0 #kJ - latentna toplina isparivanja



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
        rez = opt.minimize_scalar(interf_ht_fun, bounds=(Tl, Tv), args=(Tl, Tv, L_ekv, p, area), method='bounded')
        T_interface =rez.x
        alpha_liq = liq_alpha_interf(Tl, T_interface, L_ekv)
        alpha_vap = vap_alpha_interf(Tv, T_interface, p, L_ekv)
        resistance = (1./alpha_liq + 1./alpha_vap) / area
        heat = (Tv - Tl) / resistance * dt / 1000.0 #kJ
        #print('iteracija interface = ', iteracija)
        self.states.common.heat_interface = heat
        self.states.common.T_int = T_interface
        self.states.liquid.alpha_int = alpha_liq
        self.states.vapor.alpha_int = alpha_vap
        self.states.common.interface_resistance = resistance


    def update_states(self, liq_flow, vap_flow, T_liq_flow, T_vap_flow, dt):
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
        q_liq_ev = liq_dens / liq.M * 0.1 *(T_liq - T_sat) / T_sat #kmol/m^3/s
        q_liq_ev = max(0.0, q_liq_ev)
        q_vap_cond = vap_dens / vap.M * 0.1 *(T_sat - T_vap) / T_sat  #kmol/m^3/s
        q_vap_cond= max(0.0, q_vap_cond)
        q_ev = q_liq_ev - q_vap_cond
        liq_evap = q_ev * V_interface * dt # količina ishlapljene kapljevine, kmol
        #print("")
        #print("liq_evap "+ str(liq_evap))
        if T_liq_flow == -1: #odvodi se kapljevina stanja spremnika
            liq_flow_mol_ent = self.states.liquid.molar_enthalpy
        else: #dovodi se kapljevina drugog stanja
            liq_flow_mol_ent = liq.enthalpy(T_liq_flow)
        if T_vap_flow == -1: #odvodi se para stanja spremnika
            vap_flow_mol_ent = self.states.vapor.molar_enthalpy
        else: #dovodi se para drugog stanja
            vap_flow_mol_ent = vap.sat_vap_enthalpy(T_vap_flow) #  TO DOOOO!!!
        sat_vap_mol_ent =vap.sat_vap_enthalpy(T_sat) #specifična entalpija szp kod isparavanja, kJ/kmol
        liq_flow_quantity = liq_flow * dt #kmol
        vap_flow_quantity = vap_flow * dt #kmol
        self.states.liquid.flow_quantity = liq_flow_quantity
        self.states.vapor.flow_quantity = vap_flow_quantity
        liq_quantity = liq_N + liq_flow_quantity - liq_evap #količina kapljevine nakon isparavanja i oduzimanja, kmol
        vap_quantity = vap_N + vap_flow_quantity + liq_evap #količina pare nakon isparavanja i oduzimanja, kmol
        #entalpija kapljevine
        liq_surf_area = self.states.liquid.surf_area # povrsina oplakana kapljevinom, m^2
        vap_surf_area =  self.states.vapor.surf_area #  povrsina oplakana parom, m^2
        heat_gain_liq =liq_surf_area * self.params.k_liq * (self.params.T_amb -  T_liq) * dt / 1000.0 #toplina izvana prema kapljevini, kJ
        #print("")
        #print("heat_gain_liq "+ str(heat_gain_liq))
        heat_gain_vap =vap_surf_area * self.params.k_vap * (self.params.T_amb -  T_vap) * dt / 1000.0 #toplina izvana prema pari, kJ
        #print("")
        #print("heat_gain_vap "+ str(heat_gain_vap))
        self.interface_heat_flow(dt)
        heat_interface = self.states.common.heat_interface
        #print("")
        #print("heat_interface "+ str(heat_interface))
        # entalpija kapljevine
        liq_flow_heat = liq_flow*liq_flow_mol_ent * dt
        vap_flow_heat = vap_flow*vap_flow_mol_ent * dt
        evaporation_heat = liq_evap*sat_vap_mol_ent
        liq_ent = self.states.liquid.enthalpy + liq_flow_heat - evaporation_heat + heat_gain_liq + heat_interface
        #entalpija pare
        vap_ent = self.states.vapor.enthalpy + vap_flow_heat + evaporation_heat + heat_gain_vap - heat_interface
        liq_mol_ent = liq_ent / liq_quantity #kJ/kmol
        vap_mol_ent = vap_ent / vap_quantity #kJ/kmol
        liq_temperature = liq.temp_ent(liq_mol_ent)
        liq_density = liq.density(liq_temperature) #kg/m^3
        self.states.liquid.temperature = liq_temperature
        self.states.liquid.density = liq_density
        liq_volume = liq_quantity * liq.M / liq_density # m^3
        self.states.liquid.vol = liq_volume
        vapor_volume = self.params.int_volume - liq_volume
        self.states.vapor.vol = vapor_volume
        self.states.liquid.vol_ratio = liq_volume / self.params.int_volume
        self.states.liquid.quantity = liq_quantity
        self.states.vapor.quantity = vap_quantity
        vap_mol_dens =  (vap_quantity * 1000.0) / vapor_volume

        self.states.vapor.molar_density = vap_mol_dens
        #liq volime - vap volume - vap density + vap molar ehthalpy - pressure i temperatura
        self.states.vapor.molar_enthalpy = vap_mol_ent
        vap_temperature, pressure = vap.states_ent_dens(vap_mol_ent, vap_mol_dens)
        T_sat = liq.saturation_temperature(pressure)


        self.states.liquid.enthalpy =liq_ent
        self.states.vapor.enthalpy = vap_ent
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

        self.save.save_states(self.states) # spremi stanja

def angle_func(alpha, const):
    er = alpha - np.sin(alpha) - const
    error = er**2
    return error

#%%
if __name__ == "__main__":
#     #primjer koristenja:



    int_length = 0.909 #m
    int_diameter = 0.6 #m
    thickness = 0.05 #m
    k_liq = 1/(0.3/0.015+1/2000+1/20) #W/mK #excel
    k_vap = 1/(0.3/0.015+1/200+1/20) #W/mK #excek
    k_int = 3.0
    pressure = 283000. #Pa
    temperature_l = 125.7 #K
    temperature_v = 125.8 #K

    liq_vol_ratio = 0.75 #
    T_amb = 293.15 #K
    tank_params = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, k_int, T_amb)
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

#%%

    times =np.array(range(45000)) *10 / 60 / 60
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

    plt.figure()
    plt.plot(tank.save.vapor['vol'])
    hf = np.array(tank.save.liquid['surf_heat_flow'])*1000/10
    hfv = np.array(tank.save.vapor['surf_heat_flow'])*1000/10
    hf_tot =hf + hfv
    plt.plot(hf_tot)


    dens = tank.states.vapor.molar_density
    #liq volime - vap volume - vap density + vap molar ehthalpy - pressure i temperatura
    ent = tank.states.vapor.molar_enthalpy
    vol = tank.states.vapor.vol
