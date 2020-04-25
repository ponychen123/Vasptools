#!/bin/bash
#20190528 fix a bug in prefactor of POSCAR
#20190524 fix a bug
#a simple bash shell to add a series of specific strain on the cell
#warning: you should use the direct coordition 
#warning: your cell should better to be orthogonallty
#ponychen 2019/05/23

orientation="XY" #set which type of strain
initial=0.0      #set the initial strain
step=0.1         #set the step size
num=5            #how much strain to apply

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

#unify the strain type
unify

cp POSCAR POSCAR.orig

for((i=0;i<=num;i++))
do
	eval $(awk -v i=$i -v step=$step -v initial=$initial '
	BEGIN {strain=initial+i*step;
	       printf("strain=%f;",strain)}'
	)
	case $orientation in
		"XX")
			Setstrain=( $strain 0.0 0.0 0.0 0.0 0.0 )
			;;
		"YY")
			Setstrain=( 0.0 $strain 0.0 0.0 0.0 0.0 )
			;;
		"ZZ")
			Setstrain=( 0.0 0.0 $strain 0.0 0.0 0.0 )
			;;
		"XY")
			Setstrain=( 0.0 0.0 0.0 $strain 0.0 0.0 )
			;;
		"XZ")
			Setstrain=( 0.0 0.0 0.0 0.0 $strain 0.0 )
			;;
		"YZ")
			Setstrain=( 0.0 0.0 0.0 0.0 0.0 $strain )
			;;
		*)
			echo "please rightly set your orientation value!!!"
			exit 1
			;;
	esac
	cp POSCAR.orig POSCAR
	args=`echo ${Setstrain[*]}`
	deform ${args}
	mkdir $strain
	cp POSCAR INCAR POTCAR KPOINTS ./$strain
done
	


