#!/usr/bin/python
from __future__ import division
from __future__ import with_statement
import os
import sys
import toolbox_idynomics
import toolbox_plotting


path = sys.argv[1]
sim_dir = toolbox_idynomics.SimulationDirectory(path)

color_dict_path = os.path.join(sim_dir.figures_dir, 'color_info.txt')
if os.path.isfile(color_dict_path):
    species_color_dict = toolbox_idynomics.read_color_dict(color_dict_path)
else:
    species_color_dict = toolbox_idynomics.get_default_species_colors(sim_dir)
    toolbox_idynomics.save_color_dict(species_color_dict, color_dict_path)

fig = toolbox_plotting.PlosFigure()
axis = fig.add_subplot('', 111)

plot_left, plot_right, text_left = 30, 35, 38
y, y_diff = 100, -60

max_time = 1
max_abundance = 1

for species_name in sim_dir.get_species_names():
    time = []
    glycogens = []
    solventogens  = []
    sporulatings = []
    spores = []
    for iter_info in sim_dir.get_iterate_information():
        time.append(iter_info.time)
        max_time = max(max_time, iter_info.time)
        species = iter_info.agent_output.get_species_by_name(species_name)
        glycogens.append(0)
        solventogens.append(0)
        sporulatings.append(0)
        spores.append(0)
        for cell in species.members:
            spStatus = cell.vars['sporeStatus']
            mbStatus = cell.vars['metabolismStatus']
            if spStatus == 'spore':
                spores[-1] += 1
                continue
            if spStatus == 'sporulating':
                sporulatings[-1] += 1
                continue
            if mbStatus == 'glycolysis':
                glycogens[-1] += 1
                continue
            if mbStatus == 'solventogenesis':
                solventogens[-1] += 1
                continue
    max_glycogens = max(glycogens)
    max_solventogens = max(solventogens)
    max_sporulatings = max(sporulatings)
    max_spores = max(spores)
    max_abundance = max(max_abundance, max_glycogens, max_solventogens,
                                        max_sporulatings, max_spores)
    if max_glycogens > 0:
        axis.plot(time, glycogens, 
                      color=species_color_dict[species_name+"_glycolysis"])
    if max_solventogens > 0:
        axis.plot(time, solventogens, 
                    color=species_color_dict[species_name+"_solventogenesis"])
    if max_sporulatings > 0:
        axis.plot(time, sporulatings, 
                    color=species_color_dict[species_name+"_sporulating"])
    if max_spores > 0:
        axis.plot(time, spores, 
                    color=species_color_dict[species_name+"_spore"])

axis.set_xlim([0, max_time])
axis.set_ylim([0, max_abundance])

axis.set_xlabel('Time (h)')
axis.set_ylabel('Number of cells')

fig.process_subplots()
fig.subplots_adjust(left=0.2, bottom=0.15, right=0.98, top=0.98)
fig.save(os.path.join(path, 'figures', 'total_species_abundance.png'))
