#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import re


# Where to find the SpeciesParam subclass
file_dir = os.path.join('src',  'simulator', 'agent', 'zoo')
agent_file_name = 'Clostridium'
param_file_name = 'ClostridiumParam'
file_suffix = '.java'


# Set this to true when you're happy the values should be changed on file
overwrite_file = False


# Units you want to use as defaults, and the relationship between default
# units required and SI units (g, s, M, L)
mass_unit = 'fg'
mass_scalar = 1e-15
time_unit = 'h'
time_scalar = 3600
concn_unit = ' M'
concn_scalar = 1e0
vol_unit = ' L'
vol_scalar = 1e0

#############################################################################

convert = {
# Mass
' g' : [mass_scalar, mass_unit],
'mg' : [mass_scalar*1e3, mass_unit],
'ug' : [mass_scalar*1e6, mass_unit],
'ng' : [mass_scalar*1e9, mass_unit],
'pg' : [mass_scalar*1e12, mass_unit],
'fg' : [mass_scalar*1e15, mass_unit],

# Time
'h' : [time_scalar/3600, time_unit],
's' : [time_scalar, time_unit],
'ms' : [time_scalar*1e3, time_unit],

# Concentration
' M' : [concn_scalar, concn_unit],
'mM' : [concn_scalar*1e3, concn_unit],
'uM' : [concn_scalar*1e6, concn_unit],
'nM' : [concn_scalar*1e9, concn_unit],
'pM' : [concn_scalar*1e12, concn_unit],
'fM' : [concn_scalar*1e15, concn_unit],

# Volume
' L' : [vol_scalar, vol_unit],
'mL' : [vol_scalar*1e3, vol_unit],
'uL' : [vol_scalar*1e6, vol_unit],
'nL' : [vol_scalar*1e9, vol_unit],
'pL' : [vol_scalar*1e12, vol_unit],
'fL' : [vol_scalar*1e15, vol_unit],

}

def convert_param(value, units):
    '''
    Both current_value and current_units should be given as Strings, and will
    be returned as such. 
    '''
    units = units.replace('\n', '')
    float_val = float(value)
    for unit, [scalar, default] in convert.iteritems():
        if unit == default:
            continue
        power = ''
        try:
            index = units.index(unit) + len(unit)
            while index < len(units) and units[index] is not ' ':
                power += units[index]
                index += 1
            if power is '':
                power = '1'
            float_val /= scalar**float(power)
            units = units.replace(unit, default)
        except ValueError:
            pass
    return ' '+str(float_val), units+'\n'



''' Parameter file '''

param_file_path = os.path.join(file_dir, param_file_name+file_suffix)

with open(param_file_path, 'r') as f:
    lines = f.readlines()


for i in range(len(lines)):
    line = lines[i]
    if param_file_name+'()' in line:
        break
    if 'public double' in line.lower():
        [start, name, value, units] = re.split('ouble|=|;', line)
        new_val, new_units = convert_param(value, units)
        if new_units != units:
            print('%s has changed\n  Old: %s %s\n  New: %s %s\n'
                    %(name, value, units.replace('\n',''),
                      new_val, new_units.replace('\n','')))
            lines[i] = line.replace(units, new_units).replace(value, new_val)


if overwrite_file:
    print 'Overwriting parameter file...'
    with open(param_file_path, 'w') as f:
        f.write(''.join(lines))



''' Agent file '''

agent_file_path = os.path.join(file_dir, agent_file_name+file_suffix)

with open(agent_file_path, 'r') as f:
    lines = f.readlines()

current_name = ''
for i in range(len(lines)):
    line = lines[i]
    if '}' in line:
        break
    if 'All protein concentrations are in:' in line:
        [start, old_concn] = re.split(':', line)
        old_concn = old_concn.replace('\n', '')
        if old_concn == concn_unit:
            break
        else:
             lines[i] = line.replace(old_concn, concn_unit)
    if '_proteinNames[' in line:
        [start, value, units] = re.split('=|;', line)
        current_name = value
    if '_proteinLevels[' in line:
        [start, value, units] = re.split('=|;', line)
        units = old_concn
        new_val, new_units = convert_param(value, units)
        if new_units != units:
            print('%s has changed\n  Old: %s %s\n  New: %s %s\n'
                    %(current_name, value, old_concn.replace('\n',''),
                      new_val, new_units.replace('\n','')))
            lines[i] = line.replace(units, new_units).replace(value, new_val)


if overwrite_file:
    print 'Overwriting agent file...'
    with open(agent_file_path, 'w') as f:
        f.write(''.join(lines))
