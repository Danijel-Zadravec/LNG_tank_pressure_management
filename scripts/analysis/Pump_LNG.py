# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 11:05:54 2022

@author: faksH
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
	sys.path.insert(0, str(ROOT))


import src.properties.liq_props as liq
from src.System_Pbu import Engine,create_engine_inputs
from src.rezimi import Rezimi
from configurations.system_Pump import create_system_Pump

rezim = Rezimi()

#           port              man 1              sail             man 2             port               man 3             sail             man 4              port              man 5            sail              man 6              port
times = [0.01, 1*24*60*60-60*30, 1*24*60*60+30*60, 2*24*60*60-60*30, 2*24*60*60+60*30, 3*24*60*60-60*30, 3*24*60*60+60*30, 4*24*60*60-60*30, 4*24*60*60+60*30, 5*24*60*60-60*30, 5*24*60*60+60*30, 6*24*60*60-60*30,  6*24*60*60+60*30]
#times = [0.01, 200, 500, 2*60*60*2, 2*60*60*3]
demands = [rezim.port_winter, rezim.manouvering_winter, rezim.service_winter, rezim.manouvering_winter, rezim.port_winter, rezim.manouvering_winter, rezim.service_winter, rezim.manouvering_winter, rezim.port_winter,  rezim.manouvering_winter, rezim.service_winter, rezim.manouvering_winter] #2. 0.03

demands = [i*0.59 for i in demands]

pres = 600000.0 #6 bar -prema dokumentu  P1381 Inquiry No.02-707-000 LNG System-Final.pdf
pressures = [pres]*len(demands)
temp = liq.saturation_temperature(pres) + 5.0
temperatures_glyc = [273.15+30.0]*len(demands) #  -prema dokumentu  P1381 Inquiry No.02-707-000 LNG System-Final.pdf
(times_inp, demands_inp, pressures_inp, temperatures_inp) = create_engine_inputs(times, demands, pressures, temperatures_glyc)


slosh_times = []
evap_flow = 11.0
super_flow = 12.0

system_p = create_system_Pump(slosh_times, evap_flow, super_flow)

time_end = 24*60*6*3-1
engine = Engine(demands_inp[0:], pressures_inp[0:], temperatures_inp[0:], times_inp[0:])

system_p.calculate(engine)

res_all_p = system_p.results
res_all_p_head = res_all_p.head()
times_used_p = res_all_p[' ']['time'].to_numpy()
max_time_p = times_used_p[-1]

times = times_used_p/60/60/24 #dan
pressures = np.array(system_p.tank.save.common['pressure']) / 100000.0 #bar
fill_pct = np.array(system_p.tank.save.liquid['vol_ratio']) *100 #%
flow_liq = np.array(system_p.tank.save.liquid['flow'])* 16.04*3600 #kg/h
flow_vap = np.array(system_p.tank.save.vapor['flow']) * 16.04*3600 #kg/h
BOG_ex = system_p.BOG_excess*16.04

# --- Save results to CSV
BOG_ex = np.array(BOG_ex)
df = pd.DataFrame({
	'time_s': times_used_p,
	'day': times,
	'pressure_bar': pressures,
	'fill_pct': fill_pct,
	'liquid_flow_kg_per_h': flow_liq,
	'vapor_flow_kg_per_h': flow_vap,
	'bog_excess_kg': BOG_ex,
})

out_dir = ROOT / 'results'
out_dir.mkdir(parents=True, exist_ok=True)
out_file = out_dir / 'Pump_LNG_results.csv'
df.to_csv(out_file, index=False)
print(f"Saved results to {out_file}")
