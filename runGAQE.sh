#!/bin/bash

elementName=Cu
atomicWeight=26.981
pseudoDir='../'
pseudoPotential='Cu.pbe-dn-rrkjus_psl.1.0.0.UPF'
prefixName='Cu-4-02'
#zeropopu='../history_ga-4.json'
ismagnetic=1

numIndvs=5
numAtoms=4
numChildren=4
mutationRate=0.2
nGenerations=2
rScale=1.5
rhoScaling=3
numNodes=1

#python3 espresso_ga.py -prefix ${prefixName} -getZeroPopu ${zeropopu} -num_indvs ${numIndvs} -num_atoms ${numAtoms} -num_children ${numChildren} -mutation_rate ${mutationRate} -n_generations ${nGenerations} -r_scale ${rScale} -rho_scaling ${rhoScaling} -ele_name ${elementName} -atom_weight ${atomicWeight} -pseudoDir ${pseudoDir} -pseudo ${pseudoPotential} -nNodes ${numNodes} -ifMagnetic ${ismagnetic} > salida_nohup.dat
python3 espresso_ga.py -prefix ${prefixName} -num_indvs ${numIndvs} -num_atoms ${numAtoms} -num_children ${numChildren} -mutation_rate ${mutationRate} -n_generations ${nGenerations} -r_scale ${rScale} -rho_scaling ${rhoScaling} -ele_name ${elementName} -atom_weight ${atomicWeight} -pseudoDir ${pseudoDir} -pseudo ${pseudoPotential} -nNodes ${numNodes} -ifMagnetic ${ismagnetic} > salida_nohup.dat