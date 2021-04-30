#!/bin/bash
#20200414 final version.
#a simple bash shell to perform ideal tansile or shear process on the direct coordination,
#and get the stress-strain curve
#usage:prepare needed file for VASP, set ISIF=2, run this script, use stress_strain to plot curve
#warning: you should use the direct coordition 
#warning: your cell should better to be orthogonallty
#ponychen 2019/05/23
#email: 18709821294@outlook.com

#submission scripts related parameters

##PBS related parameters
##PBS -L nodes=6:ppn=12
##PBS -L walltime=3:00:00
##PBS -V
##PBS -N test

##SBATCH related parameters
##SBATCH -p v3_64
##SBATCH -N 2
##SBATCH -n 48

#change following parameters as you like
orientation="XX"  #set which type of strain, you should set mystrain blank!!!
initial=0.00       #set the initial strain
step=0.01         #set the step size
num=10             #how much strain to apply
mpiexec="vasp_std" #command to run vasp, modified by yourself
#mpiexec="mpirun -np 12 vasp_std"
#mpiexec="yhrun vasp_std"
mystrain=() #set all elements the strain tensor XX YY ZZ XY XZ YZ, if specified, orientation value is ommtied.

#do not change following codes unless you know what you are doing

#some functions

deform(){
#this function apply the specific strain on cell
local newstrain
newstrain=(` echo "$@" `)

eval $(awk -v arr1="${newstrain[*]}" '
           BEGIN{split(arr1, adds, " ");
	             adds[1]+=1;
			     adds[2]+=1;
			     adds[3]+=1;
			     adds[4]/=2;
			     adds[5]/=2;
			     adds[6]/=2;}
		   NR>=3 && NR<=5 {r0[NR]=$1;s0[NR]=$2;t0[NR]=$3}
		   END{r1[3]=adds[1]*r0[3]+adds[4]*r0[4]+adds[5]*r0[5];
		       s1[3]=adds[1]*s0[3]+adds[4]*s0[4]+adds[5]*s0[5];
			   t1[3]=adds[1]*t0[3]+adds[4]*t0[4]+adds[5]*t0[5];
			   r1[4]=adds[4]*r0[3]+adds[2]*r0[4]+adds[6]*r0[5];
			   s1[4]=adds[4]*s0[3]+adds[2]*s0[4]+adds[6]*s0[5];
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

unify(){
	#this function unify the prefactor of POSCAR to 1.000
	eval $(awk '
	    NR==2 {scal=$1}
		NR>=3 && NR<=5 {r[NR]=$1*scal;s[NR]=$2*scal;t[NR]=$3*scal}
	END{for(i=3;i<=5;i++){
	    printf("newaxis[%s]=\"%9.6f\t%9.6f\t%9.6f\n\";",i,r[i],s[i],t[i])}
	}' POSCAR)
    sed -i "2c 1.0000"      POSCAR
	sed -i "3c${newaxis[3]}" POSCAR
	sed -i "4c${newaxis[4]}" POSCAR
	sed -i "5c${newaxis[5]}" POSCAR
}

# check whether new starting or restarting, added by xuefei Liu
if [ -f "./POSCAR.orig" ];then
	echo "Restarting !!"
	recal=true
else
        cp POSCAR POSCAR.orig
	echo "New starting !!"
	#unify the prefactor of POSCAR
	echo "unify the prefactor of POSCAR!"
	unify
	recal=false

fi

#delete old output files
if $recal;then
	echo "restarting..we will not delete the old stress_strain.all files."
else
        if [ ! -f "./engineering_stress_strain.all" ];then
	        echo "engineering_stress_strain.all is not exsited"
        else
	        echo "rm the old engineering_stress_strain.all"
	        rm ./engineering_stress_strain.all
        fi
        if [ ! -f "./true_stress_strain.all" ];then
	        echo "true_stress_strain.all is not exsited"
        else
	        echo "rm the old true_stress_strain.all"
	        rm ./true_stress_strain.all
        fi
fi

#start calculation

echo "starting calculation, hold on, drink coffee and sleeping..."

#get the strain type

for((i=0;i<=num;i++))
do
    if  test -z "${mystrain[*]}" ;then
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
			col=8
			;;
		"YZ")
			Setstrain=( 0.0 0.0 0.0 0.0 0.0 $strain )
			col=7
			;;
		*)
			echo "please rightly set your orientation value!!!"
			exit 1
			;;
	esac
    else
    eval $(awk -v arr1="${mystrain[*]}" -v i=$i -v s=$step -v ini=$initial'
	BEGIN {split(arr1, adds, " ");
          for(j=1;j<=6;j++){
          adds[j]=ini+adds[j]*i*s;}
          printf("Setstrain=( %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f );",adds[1],adds[2],adds[3],adds[4],adds[5],adds[6])}' )
    fi
	cp POSCAR.orig POSCAR
	args=`echo ${Setstrain[*]}`
	deform ${args}
	echo "iretation=$i, present applied strain is ${Setstrain[*]}"
	$mpiexec > /dev/null 2>&1
	#get the result
	if test -z "${mystrain[*]}";then
    grep " in kB" OUTCAR | tail -1 | awk -v col=$col -v strain=$strain '{print strain, -$col/10}' >> engineering_stress_strain.all
	grep " in kB" OUTCAR | tail -1 | awk -v col=$col -v strain=$strain '
	    {estress=-$col/10;
			tstress=estress*(1+strain);
			tstrain=log(1+strain);
			print tstrain, tstress}' >> true_stress_strain.all
		else
        grep " in kB" OUTCAR | tail -1 | awk -v num=$i -v arr1="${mystrain[*]}" -v step=$step '
		BEGIN{split(arr1, strain, " ");}
		{for(i=3;i<=8;i++){
			stress[i-2]=-$i/10}
		m=0;n=0;s=0;
		for(i=1;i<=6;i++){
			m+=stress[i]*strain[i]*step;
			n+=strain[i]*step;
			s+=(strain[i]*step)^2;}
		stresstol=m/n;
		straintol=sqrt(s)*num;
        print straintol, stresstol;
	    }' >> engineering_stress_strain.all
        grep " in kB" OUTCAR | tail -1 | awk -v num=$i -v arr1="${mystrain[*]}" -v step=$step '
		BEGIN{split(arr1, strain, " ");}
		{for(i=3;i<=8;i++){
			stress[i-2]=-$i/10}
		m=0;n=0;s=0;
		for(i=1;i<=6;i++){
			m+=stress[i]*strain[i]*step;
			n+=strain[i]*step;
			s+=(strain[i]*step)^2;}
		stresstol=m/n;
		straintol=sqrt(s)*num;
		tstress=stresstol*(1+straintol);
		tstrain=log(1+straintol);
        print tstrain, tstress;
	    }' >> true_stress_strain.all
	fi
    conver=`grep "reached required" OUTCAR `
    #check whether present iretation converge
	if [ -n "$conver" ]; then
		echo "present iretation converged"
	else
		echo "FBI warning: present iretation not converge, you may adjust your INCAR!"
		break
	fi
done
	
echo "all the iretation have finished, the unit of stress is GPa!"

#plot the curve, of course, just have a roughly look

gnuplot -e "set term dumb; plot 'engineering_stress_strain.all'  w l"
#gnuplot -e "set term dumb; plot 'true_stress_strain.all' w l"
