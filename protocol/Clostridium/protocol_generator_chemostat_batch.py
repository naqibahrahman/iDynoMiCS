#!/usr/bin/python
from __future__ import division
from __future__ import with_statement

# Standard ki, in Molar
standard_ki = 0.1
# Fraction change for the mutant
mutant_change = 0.1

change_dict = {

'$$$OutputPeriodHours$$$' : 1.0,

'$$$TimeStepHours$$$' : 0.1,

'$$$FinalTimeDays$$$' : 7.0,

'$$$DilutionPerHour$$$' : 0.2,

'$$$InflowMolar$$$' : 0.1,

'$$$KiAcidMolar$$$' : standard_ki * (1 + mutant_change),

'$$$KiSolventMolar$$$' : standard_ki * (1 - mutant_change),

}

def prettify(key):
    return str(change_dict[key]).replace('.', 'p')

out_name = 'WTvsMut_%sDil_%sGluIn_%sKiAcid_%sKiSolv.xml' \
    %(prettify('$$$DilutionPerHour$$$'), prettify('$$$InflowMolar$$$'),
      prettify('$$$KiAcidMolar$$$'), prettify('$$$KiSolventMolar$$$'))


with open('protocol_template_chemostat_batch.xml', 'r') as f:
    script = f.read()

for key, item in change_dict.iteritems():
    script = script.replace(key, str(item))

with open(out_name, 'w') as f:
    f.write(script)