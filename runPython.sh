#!/bin/bash

elementName=Cu
atomicWeight=63.546
psuedoPotential='Cu.pbe-dn-rrkjus_psl.1.0.0.UPF'

numIndvs=5
numAtoms=4
numChildren=4
mutationRate=0.3
nGenerations=2
rScale=1.5
rhoScaling=3
numNodes=4

python3 espresso_ga.py -num_indvs ${numIndvs} -num_atoms ${numAtoms} -num_children ${numChildren} -mutation_rate ${mutationRate} -n_generations ${nGenerations} -r_scale ${rScale} -rho_scaling ${rhoScaling} -ele_name ${elementName} -atom_weight ${atomicWeight} -pseudo ${psuedoPotential} -nNodes ${numNodes} > salida_nohup.dat
