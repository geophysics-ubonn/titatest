#!/bin/bash

function run_dir(){
	dire=$1
	cd ${dire}
	clean_sim.py; cd exe; CRMod; CRTomo_dev_titan; cd ..; ../plot.sh
	td_create_overview.py
	beta=`tail -1 exe/crtomo.cfg`
	cp overview.png ../overview_${dire}_${beta}.png
	cd ..
}


run_dir 01_smooth
run_dir 02_mgs
run_dir 03_mgs
run_dir 04_mgs
run_dir 05_mgs

