#!/bin/bash
#a simple bash shell to perform ideal tansile or shear process,and get the
#stress-strain curve
#usage:prepare needed file for VASP, set ISIF=2, run this script, use stress_strain to plot curve
#warning: you should use the direct coordition 
#warning: your cell should better to be orthogonallty
#ponychen 2019/05/23
#email: 18709821294@outlook.com

#change following parameters as you like
orientation="XX" #set which type of strain
initial=0.0      #set the initial strain
step=0.01         #set the step size
num=50            #how much strain to apply
mpiexec="vasp_std"

#do not change following codes unless you know what you are doing
deform(){
#this function apply the specific strain on cell
local newstrain
newstrain=(` echo "$@" `)

eval $(awk -v arr1="${newstrain[*]}" '
           BEGIN{split(arr1, adds, " ");
	             adds[1]+=1;
			     adds[2]+=1;
			     adds[3]+=1;}
		   NR>=3 && NR<=5 {r0[NR]=$1;s0[NR]=$2;t0[NR]=$3}
		   END{r1[3]=adds[1]*r0[3]+adds[4]*r0[4]+adds[5]*r0[5];
		       s1[3]=adds[1]*s0[3]+adds[4]*s0[4]+adds[5]*s0[5];
			   t1[3]=adds[1]*t0[3]+adds[4]*t0[4]+adds[5]*t0[5];
			   r1[4]=adds[4]*r0[3]+adds[2]*r0[4]+adds[6]*r0[5];
			   s1[4]=adds[4]*so[3]+adds[2]*s0[4]+adds[6]*s0[5];
			   t1[4]=adds[4]*t0[3]+adds[2]*t0[4]+adds[6]*t0[5];
			   r1[5]=adds[5]*r0[3]+adds[6]*r0[4]+adds[3]*r0[5];
			   s1[5]=adds[5]*s0[3]+adds[6]*s0[4]+adds[3]*s0[5];
			   t1[5]=adds[5]*t0[3]+adds[6]*t0[4]+adds[3]*t0[5];
			   for(i=3;i<=5;i++){
				   printf("newpos[%s]=\" %9.6f %9.6f %9.6f\n\";",i,r1[i],s1[i],t1[i])}
			   }' POSCAR)
sed -i "3c${newpos[3]}" POSCAR
sed -i "4c${newpos[4]}" POSCAR
sed -i "5c${newpos[5]}" POSCAR
}

cp POSCAR POSCAR.orig
echo "starting calculation, hold on, drink coffee and sleeping..."

for((i=0;i<=num;i++))
do
	eval $(awk -v i=$i -v step=$step -v initial=$initial '
	BEGIN {strain=initial+i*step;
	       printf("strain=%f;",strain)}'
	)
	case $orientation in
		"XX")
			Setstrain=( $strain 0.0 0.0 0.0 0.0 0.0 )
			col=3
			;;
		"YY")
			Setstrain=( 0.0 $strain 0.0 0.0 0.0 0.0 )
			col=4
			;;
		"ZZ")
			Setstrain=( 0.0 0.0 $strain 0.0 0.0 0.0 )
			col=5
			;;
		"XY")
			Setstrain=( 0.0 0.0 0.0 $strain 0.0 0.0 )
			col=6
			;;
		"XZ")
			Setstrain=( 0.0 0.0 0.0 0.0 $strain 0.0 )
			col=7
			;;
		"YZ")
			Setstrain=( 0.0 0.0 0.0 0.0 0.0 $strain )
			col=8
			;;
		*)
			echo "please rightly set your orientation value!!!"
			exit 1
			;;
	esac
	cp POSCAR.orig POSCAR
	args=`echo ${Setstrain[*]}`
	deform ${args}
	echo "iretation=$i, present strain is $strain"
	$mpiexec > /dev/null 2>&1
    grep " in kB" OUTCAR | tail -1 | awk -v col=$col -v strain=$strain '{print strain, $col}' >> stress_strain.all
    conver=`grep "reached required" OUTCAR `
    #check whether present iretation converge
	if [ -n "$conver" ]; then
		echo "present iretation converged"
	else
		echo "FBI warning: present iretation not converge, you may adjust your INCAR!"
	fi
done
	
echo "all the iretation have finished"
#plot
gnuplot -e "set term dumb; plot 'stress_strain.all' u 1:2 w l"

