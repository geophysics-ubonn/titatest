#!/bin/bash
#----------------------------------------------------------------------------------------
#	Name		: gen_pbuilder_packages.sh
#	Description	: Generate griev .deb packages for various distributions. Use the alien package to convert those .debs to .rpms
#	Input		:
#	Output		:
#	Return		:
#----------------------------------------------------------------------------------------
# create griev packages with pbuilder
#pbuilder_target_loc="/var/cache/pbuilder"
pbuilder_target_loc="/home/mweigand/BIG_DATA/pbuilder"

if [ ! -e "$pbuilder_target_loc"/squeeze-base.tgz ]; then
	DIST=squeeze pbuilder create --basetgz  "$pbuilder_target_loc"/squeeze-base.tgz
else
	DIST=squeeze pbuilder update --basetgz  "$pbuilder_target_loc"/squeeze-base.tgz
fi

# ubuntu distributions
#for distri in karmic lucid # hardy 
#do
#	# create if necessary
#	if [ ! -e "$pbuilder_target_loc"/$distri-base.tgz ]; then
#		echo "###################################"
#		echo "## creating $distri pbuilder image ##"
#		echo "###################################"
#	
#		DIST=$distri pbuilder create --basetgz  "$pbuilder_target_loc"/$distri-base.tgz  --distribution $distri --mirror "http://archive.ubuntu.com/ubuntu/" --othermirror "deb http://archive.ubuntu.com/ubuntu/ $distri universe multiverse"
#
#	# or update
#	else
#		echo "###################################"
#		echo "## updating $distri pbuilder image ##"
#		echo "###################################"
#	
#		DIST=$distri pbuilder update --basetgz  "$pbuilder_target_loc"/$distri-base.tgz
#
#	fi
#done


## test if package was created yet
if [ -d package ]; then
	echo "## deleteting and recreating package directory ##"
	rm -rf package
fi

./prep_deb_package.sh
cd package/crtomomod*
dpkg-buildpackage
cd ../../


## build for various distributions
cd package

for distri in squeeze # karmic lucid
do
	echo "cleaning $distri build dir"
	pbuilder --clean --basetgz  "$pbuilder_target_loc"/$distri-base.tgz --buildresult  "$pbuilder_target_loc"/$distri/result
	echo "building package"
	pbuilder --build --basetgz  "$pbuilder_target_loc"/$distri-base.tgz --buildresult  "$pbuilder_target_loc"/$distri/result *.dsc

	# create rpm
	cur_dir=`pwd`
	cd  "$pbuilder_target_loc"/$distri/result
	alien --to-rpm *.deb
	cd "$cur_dir"

	# scp to malm
	scp  "$pbuilder_target_loc"/$distri/result/crtomomod* mweigand@malm:/users/mweigand/www/$distri

done


cd ..

