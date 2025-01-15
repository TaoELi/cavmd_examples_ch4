import json
import numpy as np

nmol=400
freqidx=0
single_charges = [-0.24, 0.06, 0.06, 0.06, 0.06]

freq = np.loadtxt('./freq_fix.txt')
data = {}
data["apply_photon"] = True
#data["n_modes"] = 1
data["eff_mass"] = 1.0
data["freqs_cm"] = freq[freqidx,0]
data["E0"] = freq[freqidx,1] / np.sqrt(nmol//100)
charges = single_charges*nmol

print("charges length = ", len(charges))
# I need partial charge at each atom (C, H, H, H, H)
charges = [float("%.5f" %x) for x in charges]

data["charge_array"] = charges
'''
# Also adding parameters to define the external laser
data["add_pulse_photon"] = True
data["add_pulse_direction"] = 0
data["pulse_params"] = [0.0006, 500, 1619.8, 1.14514, 10.0]
data["pulse_atoms"] = []
data["transition_photon_charge"] = 0.028
data["dt"] = 0.5
'''
with open('./photon_params.json', 'w') as outfile:
    json.dump(data, outfile)
