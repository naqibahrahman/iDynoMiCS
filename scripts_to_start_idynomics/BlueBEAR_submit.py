#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
from optparse import OptionParser
import os


### Parse options
parser = OptionParser()
parser.add_option("-H", "--Hours", dest="hours", default=0,
                                             type="int", help="hours walltime")
parser.add_option("-M", "--Minutes", dest="minutes", default=1,
                                           type="int", help="minutes walltime")
parser.add_option("-q", "--Queue", dest="queue", default="bbdefault",
                                                         help="BlueBEAR queue")
(options, args) = parser.parse_args()

for protocol in args:
    if not os.path.isfile(protocol):
        print "Protocol file doesn't exist! %s" %(protocol)
        continue
    protocol = os.path.abspath(protocol)
    protocol_dir, protocol_file = os.path.split(protocol)
    if not (isinstance(options.hours, int) or isinstance(options.minutes, int)):
        print "Time arguments must be integers!\nHours = %s\nMinutes = %s" \
            %(options.hours, options.minutes)
        print type(options.hours), type(options.minutes)
        exit()
    cores = 8 if (options.queue == 'bbtest') else 16
    script =  '''#!/bin/bash
#MOAB -N elmerJob
#MOAB -A kreftj01
#MOAB -j oe
#MOAB -q %s
#MOAB -l "walltime=%d:%d:0,nodes=1:ppn=%d"
module load apps/python2/v2.7.3
python2 RunIdyno.py %s
    ''' %(options.queue, options.hours, options.minutes, cores, protocol)
    script_path = os.path.join(protocol_dir,
                               'msub_'+protocol_file.replace('.xml', '')+'.sh')
    with open(script_path, 'w') as f:
        f.write(script)
    os.system('msub %s' %(script_path))
