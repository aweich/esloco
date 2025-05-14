#!/bin/bash

for config in case2_VCN/5_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/7_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/10_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/12_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/15_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/dominance_experiment/10dd_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/dominance_experiment/1dd_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

for config in case2_VCN/dominance_experiment/01dd_configs/*; do
	echo "Processing: $config"
	esloco "$config"
	plot_esloco "$config"
done

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

