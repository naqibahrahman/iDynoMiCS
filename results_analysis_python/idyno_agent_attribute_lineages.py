#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
from optparse import OptionParser
import toolbox_idynomics
import toolbox_plotting


parser = OptionParser()
parser.add_option("-a", "--Attribute", dest="attribute", default="biomass",
                  help="name of the attribute wanted")
parser.add_option("-u", "--Units", dest="units", default="fg", 
                                      help="units of the attribute wanted")
parser.add_option("-z", "--ZeroYAxis", dest="zero_y_axis", default=False,
                    action="store_true",
                    help="forces the lower limit of the y axis to zero")
(options, args) = parser.parse_args()

def isSameCell(early, late):
        return ( early.vars['family'] == late.vars['family'] ) and \
                ( early.vars['genealogy'] == late.vars['genealogy'] )
    
def isPredecessor(early, late):
    if not ( early.vars['family'] == late.vars['family'] ):
        return False
    early_genealogy = int(early.vars['genealogy'])
    late_genealogy = int(late.vars['genealogy'])
    late_generation = int(late.vars['generation'])
    return (late_genealogy == early_genealogy + 2**(late_generation))

for path in args:
    sim_dir = toolbox_idynomics.SimulationDirectory(path)
    color_dict_path = os.path.join(sim_dir.figures_dir, 'color_info.txt')
    if os.path.isfile(color_dict_path):
        species_color_dict= toolbox_idynomics.read_color_dict(color_dict_path)
    else:
        species_color_dict = \
                         toolbox_idynomics.get_default_species_colors(sim_dir)
        toolbox_idynomics.save_color_dict(species_color_dict, color_dict_path)
    fig = toolbox_plotting.PlosFigure()
    axis = fig.add_subplot('', 111)
    max_time = 0.0
    max_val = 0.0
    min_val = 0.0
    cells = []
    for iter_info in sim_dir.get_iterate_information():
        max_time = max(max_time, iter_info.time)
        for cell in iter_info.agent_output.get_all_cells():
            family = int(cell.vars['family'])
            genealogy = int(cell.vars['genealogy'])
            attribute_val = float(cell.vars[options.attribute])
            max_val = max(max_val, attribute_val)
            min_val = min(min_val, attribute_val)
            predecessor = None
            for other in cells:
                if isSameCell(cell, other):
                    predecessor = other
                    break
            if predecessor == None:
                cell.temp_values = [attribute_val]
                cell.times = [iter_info.time]
                cell.color = species_color_dict[cell.species]
                cells.append(cell)
            else:
                predecessor.temp_values.append(attribute_val)
                predecessor.times.append(iter_info.time)
    print len(cells)
    for cell in cells:
        if cell.vars['genealogy'] == '0':
            continue
        birthday = float(cell.vars['birthday'])
        for other in cells:
            if isPredecessor(other, cell):
                val, index = max((t, i) for (i, t) in \
                                    enumerate(other.times) if t < birthday)
                cell.times.insert(0, other.times[index])
                cell.temp_values.insert(0, other.temp_values[index])
    for cell in cells:
        axis.plot(cell.times, cell.temp_values, color=cell.color, alpha=0.3)
    if options.zero_y_axis:
        min_val = 0.0
    if max_time > 0.0:
        axis.set_xlim([0.0, max_time])
    if max_val > 0.0:
        axis.set_ylim([min_val, max_val])
    axis.set_xlabel('Time (h)')
    axis.set_ylabel('%s (%s)'%(options.attribute, options.units))
    fig.process_subplots()
    fig.subplots_adjust(left=0.23, bottom=0.2, right=0.98, top=0.95)
    fig.save(os.path.join(path, 'figures',
                              'agent_lineage_'+options.attribute+'.png'))