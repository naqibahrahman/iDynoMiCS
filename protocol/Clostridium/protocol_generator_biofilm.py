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

'$$$BulkMolar$$$' : 0.001,

'$$$KiAcidMolar$$$' : standard_ki * (1 + mutant_change),

'$$$KiSolventMolar$$$' : standard_ki * (1 - mutant_change),

'$$$RandomSeed$$$' : 42,

}

def prettify(key):
    return str(change_dict[key]).replace('.', 'p')

out_name = 'WTvsMut_%sGluBulk_%sKiAcid_%sKiSolv_%sRS.xml' \
    %(prettify('$$$BulkMolar$$$'), prettify('$$$KiAcidMolar$$$'),
          prettify('$$$KiSolventMolar$$$'), prettify('$$$RandomSeed$$$'))


with open('protocol_template_biofilm.xml', 'r') as f:
    script = f.read()

for key, item in change_dict.iteritems():
    script = script.replace(key, str(item))

with open(out_name, 'w') as f:
    f.write(script)