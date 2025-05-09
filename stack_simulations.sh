#!/bin/bash

for config in case2_VCN/leftover/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

#case2_VCN/dominance_dilution/01dd_configs/*

#for config in case2_VCN/leftover/*; do
#	echo "Processing: $config"
#	InsertionSimulation/InsertionGenomeBuilder.py "$config"
#	InsertionSimulation/simplot.py "$config"
#done
#
#for config in case2_VCN/dominance_dilution/1dd_configs/*; do
#	echo "Processing: $config"
#	InsertionSimulation/InsertionGenomeBuilder.py "$config"
#	InsertionSimulation/simplot.py "$config"
#done

