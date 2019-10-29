#!/bin/bash
#ponychen

ini_strain=1.0
fin_strain=1.5
num=20

step=`echo "scale=2;($fin_strain-$ini_strain)/($num-1)"|bc -l`
for i in $(seq $ini_strain $step $fin_strain)
do  
	if [ ! -d $i ];then
	    mkdir $i
	else
		rm -r $i/*
	fi
	cp INCAR KPOINTS POTCAR POSCAR run_vasp.sh ./$i
	sed -in "2c $i" $i/POSCAR
	rm $i/POSCARn
done

echo "done!"
