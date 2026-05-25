
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
	sys.path.insert(0, str(ROOT))
     
from src.tank import Tank, Tank_initial, Tank_params


int_length = 0.909 #m
int_diameter = 0.6 #m
thickness = 0.03 #m
k_liq = 0.056491 #W/mK #excel
k_vap = 0.056491 #W/mK #excek
pressure = 283000. #Pa
temperature_l = 125.7 #K
temperature_v = 125.8 #K

liq_vol_ratio = 0.75 #
T_amb = 297.0 #K
tank_params = Tank_params(int_length, int_diameter, thickness, k_liq, k_vap, T_amb)
tank_initial = Tank_initial(pressure, temperature_l, temperature_v, liq_vol_ratio)
tank = Tank(tank_params,tank_initial)
vap_flow = 0.0

for i in range(50400):
    tank.update_states(0.0, vap_flow, -1, -1, 10)
    if i%1000 == 0:
        print("               " + str(i) + "                    ")

times =np.array(range(50400)) *10 / 60 / 60
pressures = np.array(tank.save.common['pressure']) / 1000.0
pressures=pressures[1:]




# Import the CSV file
data = pd.read_csv('scripts\\validation\\validation_data.csv')


ts = data["t"].values
ps = data["p"].values

# Plot all on same figure
plt.figure()
plt.plot(times, pressures, label='Simulation')
plt.scatter(ts, ps, label='Experiment')
plt.grid()
plt.xlabel('Time (hours)')
plt.ylabel('Pressure (kPa)')
plt.legend()
plt.show()

# Interpolate simulation pressures at experimental times
sim_pressures_at_ts = np.interp(ts, times, pressures)

# Calculate relative deviation (%)
relative_deviation = (ps - sim_pressures_at_ts) / sim_pressures_at_ts * 100

print(relative_deviation)

avg_rel_deviation = np.mean(np.abs(relative_deviation))
print(f"Average relative deviation: {avg_rel_deviation:.2f}% for k = {k_liq:.2f}")

# Create figure
plt.rcParams['font.size'] = 11
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['grid.linewidth'] = 0.8
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2

fig, ax = plt.subplots(figsize=(8, 5), dpi=100)
ax.plot(times, pressures, label='Simulation', linewidth=2, color='#1f77b4')
ax.scatter(ts, ps, label='Experiment', s=50, color='#ff7f0e', alpha=0.7, edgecolors='black', linewidth=0.5)

ax.set_xlabel('Time (hours)', fontsize=12, fontweight='bold')
ax.set_ylabel('Pressure (kPa)', fontsize=12, fontweight='bold')
ax.set_title('Tank Pressure Comparison', fontsize=13, fontweight='bold', pad=15)
ax.axhline(y=1585, color='red', linestyle='--', linewidth=1.5, label='MAWP', alpha=0.8)
ax.grid(True, alpha=0.3, linestyle='--')
ax.legend(loc='best', frameon=True, shadow=True, fontsize=11)


plt.tight_layout()

# Save at 600 dpi
plt.savefig('scripts\\validation\\tank_pressure_comparison.png', dpi=600, bbox_inches='tight')
plt.savefig('scripts\\validation\\tank_pressure_comparison.pdf', dpi=600, bbox_inches='tight')
print("Figure saved as: tank_pressure_comparison.png and tank_pressure_comparison.pdf")
plt.show()

