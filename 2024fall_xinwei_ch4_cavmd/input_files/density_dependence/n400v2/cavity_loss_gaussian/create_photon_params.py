import json
import numpy as np

nmol=400
single_charges = [-0.24, 0.06, 0.06, 0.06, 0.06]

data = {}
data["apply_photon"] = True
#data["n_modes"] = 1
data["eff_mass"] = 1.0
data["freqs_cm"] = 1311.2
data["E0"] = 0.0000
charges = single_charges*nmol

print("charges length = ", len(charges))
# I need partial charge at each atom (C, H, H, H, H)
charges = [float("%.5f" %x) for x in charges]

data["charge_array"] = charges
'''
freq = np.loadtxt('/work/taoeli/users/astolfo/project/ch4_cavity/origin_input/freq.txt')
for e0 in range(np.shape(freq)[0]):
    if abs(freq[e0,0] - data["E0"]) < 1e-8 : 
        UP_freq = freq[e0,2]
        print(f'For E0={freq[e0,0]}, exciting with UP_freq={UP_freq}')
# Also adding parameters to define the external laser
data["add_pulse_photon"] = True
data["add_pulse_direction"] = 0
data["pulse_params"] = [0.0006, 500, UP_freq, 1.14514, 10.0]
data["pulse_atoms"] = []
data["transition_photon_charge"] = 0.028
data["dt"] = 0.5
'''
with open('./photon_params.json', 'w') as outfile:
    json.dump(data, outfile)
