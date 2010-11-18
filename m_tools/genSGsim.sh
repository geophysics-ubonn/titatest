#!/bin/bash
####################################
# variables for SGsim:
###################################
#########
# condition contains file name with conditoning data 
#		(x y log_10(rho) log_10(rho)-mean)
#########
# Grid specefications
# nx number of x cells
# x0 x-offset (m)
# dx x increment (m)
# ny number of y cells
# y0 y-offset (m)
# dy y increment (m)
#########
# seed random number seed
######### 
# nugget nugget effect for variogram model
# it variogram model (1 = exponential, 2 = spherical, 3 = gauss)
# variance variance of the model
# angle dip angle (for anisotropy..)
# ax x-correlation length (m)
# ay y-correlation length (m)
########
write_sgsim(){
	echo "  Parameters for SGSIM" >> $tmpfl 
	echo " ********************"  >> $tmpfl 
	echo " "  >> $tmpfl 
	echo " START OF PARAMETERS:"  >> $tmpfl 
	echo " $condition  -condition data file"  >> $tmpfl 
	echo " 1  2  3  4  0  0              -  columns for X,Y,Z,vr,wt,sec.var."  >> $tmpfl 
	echo " -1.e20       1.0e21             -  trimming limits"  >> $tmpfl 
	echo " 0                             -transform the data (0=no, 1=yes)"  >> $tmpfl 
	echo " sgsim.trn                     -  file for output trans table"  >> $tmpfl 
	echo " 0                             -  consider ref. dist (0=no, 1=yes)"  >> $tmpfl 
	echo " histsmth.out                  -  file with ref. dist distribution"  >> $tmpfl 
	echo " 1  2                          -  columns for vr and wt"  >> $tmpfl 
	echo " 0.0    15.0                   -  zmin,zmax(tail extrapolation)"  >> $tmpfl 
	echo " 1       0.0                   -  lower tail option, parameter"  >> $tmpfl 
	echo " 1      15.0                   -  upper tail option, parameter"  >> $tmpfl 
	echo " 1                             -debugging level: 0,1,2,3"  >> $tmpfl 
	echo " sgsim.dbg                     -file for debugging output"  >> $tmpfl 
	echo " sgsim.out                     -file for simulation output"  >> $tmpfl 
	echo " 1                             -number of realizations to generate"  >> $tmpfl 
	echo " $nx    $x0    $dx              -nx,xmn,xsiz"  >> $tmpfl 
	echo " $ny    $y0    $dy              -ny,ymn,ysiz"  >> $tmpfl 
	echo " 1     0.5    1.0              -nz,zmn,zsiz"  >> $tmpfl 
	echo " $seed                         -random number seed"  >> $tmpfl 
	echo " 0     7                       -min and max original data for sim"  >> $tmpfl 
	echo " 12                            -number of simulated nodes to use"  >> $tmpfl 
	echo " 1                             -assign data to nodes (0=no, 1=yes)"  >> $tmpfl 
	echo " 0     3                       -multiple grid search (0=no, 1=yes),num"  >> $tmpfl 
	echo " 10                             -maximum data per octant (0=not used)"  >> $tmpfl 
	echo " 1.0  1.0  0.0              -maximum search radii (hmax,hmin,vert)"  >> $tmpfl 
	echo "  0.0   0.0   0.0              -angles for search ellipsoid"  >> $tmpfl 
	echo " 0     0.60   1.0              -ktype: 0=SK,1=OK,2=LVM,3=EXDR,4=COLC"  >> $tmpfl 
	echo " ../data/ydata.dat             -  file with LVM, EXDR, or COLC variable"  >> $tmpfl 
	echo " 4                             -  column for secondary variable"  >> $tmpfl 
	echo " 1    $nugget                      -nst, nugget effect"  >> $tmpfl 
	echo " $it    $variance  $angle   0.0   0.0     -it,cc,ang1,ang2,ang3"  >> $tmpfl 
	echo "          $ax  $ay  0.0     -a_hmax, a_hmin, a_vert" >> $tmpfl 
}

get_cond_data(){
	echo $1
	cp cur_crmod_mag.pdf cond
	cd $cur2/$d
	fn_in=`ls *$d.mag`
	fn_out='Cond_'$d'.dat'
	rm -f $fn_out
	if [ -e "$fn_in" ];then
		echo "$fn_in exists"
	else
		echo "$fn_in is no valid file"
	fi
	for x in $xlist;do
		for y in $ylist;do
			echo "cheking value $x $y"
			val=`awk -v x=$x -v dx=$dx -v y=$y -v dy=$dy -v m=$mean '{ if( ( sqrt( ($1-x)^2 ) < dx ) && (sqrt( ($2-y)^2 ) < dy ) ) { print $1"  "$2"  "$3"  "$3-m } }' $fn_in`
		        if [ -n "$val" ];then
			       echo "found $val"
			       echo "$val" >> tmp.dat # remove mean from log values
			fi
		done
	done
	n=`wc -l tmp.dat|awk '{print $1}'`
       	echo "found $n matches"
	echo $n	> tmp.head
	echo '4' >> tmp.head
	echo 'x-coord' >> tmp.head
	echo 'y-coord' >> tmp.head
	echo 'z-coord' >> tmp.head
	echo 'resisitivity (Ohm m)' >> tmp.head
	cat tmp.head tmp.dat > $cur2/$fn_out
	rm tmp.head tmp.dat
	cd $cur2
}

seed1=`date +%N|awk '{printf("%d",$1*2+1)}'` # general too big number
seed2=`date +%N|awk '{printf("%d",$1/1.e3)}'` # truncate some digits
seed2=`echo $seed2|awk '{printf("%d",$1*2+1)}'` # odd
seed=$seed2
echo $seed1 $seed2
