#!/bin/bash

elementName=Al
atomicWeight=26.981
pseudoDir='../../'
pseudoPotential='Al.pbe-n-rrkjus_psl.1.0.0.UPF'
prefixName='Al-35-02'

numIndvs=5
numAtoms=35
numChildren=4
mutationRate=0.2
nGenerations=60
rScale=1.5
rhoScaling=3
numNodes=20

python3 espresso_ga.py -prefix ${prefixName} -num_indvs ${numIndvs} -num_atoms ${numAtoms} -num_children ${numChildren} -mutation_rate ${mutationRate} -n_generations ${nGenerations} -r_scale ${rScale} -rho_scaling ${rhoScaling} -ele_name ${elementName} -atom_weight ${atomicWeight} -pseudoDir ${pseudoDir} -pseudo ${pseudoPotential} -nNodes ${numNodes} > salida_nohup.dat
