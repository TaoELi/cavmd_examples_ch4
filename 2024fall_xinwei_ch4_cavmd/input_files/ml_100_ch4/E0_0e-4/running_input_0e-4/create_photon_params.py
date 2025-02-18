import json
import numpy as np

nmol=100
single_charges = [-0.24, 0.06, 0.06, 0.06, 0.06]

data = {}
data["apply_photon"] = True
#data["n_modes"] = 1
data["eff_mass"] = 1.0
data["freqs_cm"] = 1254.4
data["E0"] = 0.0
charges = single_charges*nmol

print("charges length = ", len(charges))
# I need partial charge at each atom (C, H, H, H, H)
charges = [float("%.5f" %x) for x in charges]

data["charge_array"] = charges
'''
# Also adding parameters to define the external laser
data["add_cw"] = True
data["add_cw_direction"] = 0
data["cw_params"] = [0.0006, 1581.4, 1.14514, 100.0, 600.0]
data["cw_atoms"] = [-1]
data["dt"] = 0.5
'''
with open('./photon_params.json', 'w') as outfile:
    json.dump(data, outfile)
