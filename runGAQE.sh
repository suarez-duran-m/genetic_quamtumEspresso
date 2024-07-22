#!/bin/bash

elementName=Cu
atomicWeight=63.546
pseudoDir='../'
pseudoPotential='Cu.pbe-dn-rrkjus_psl.1.0.0.UPF'
prefixName='Cu-10-A'

numIndvs=10
numAtoms=5
numChildren=8
mutationRate=0.1
nGenerations=20
rScale=1.5
rhoScaling=3
numNodes=20

python3 espresso_ga.py -prefix ${prefixName} -num_indvs ${numIndvs} -num_atoms ${numAtoms} -num_children ${numChildren} -mutation_rate ${mutationRate} -n_generations ${nGenerations} -r_scale ${rScale} -rho_scaling ${rhoScaling} -ele_name ${elementName} -atom_weight ${atomicWeight} -pseudoDir ${pseudoDir} -pseudo ${pseudoPotential} -nNodes ${numNodes} > salida_nohup.dat
