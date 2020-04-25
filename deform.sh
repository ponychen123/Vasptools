#!/bin/bash
#2020/04/25 fix a definition of shear strain
#20190524 fix a bug
#a simple bash shell to add specific strain on the cell
#warning: you should use the direct coordition 
#ponychen 2019/05/22

#the XX, YY, ZZ, XY, XZ, YZ portion of the strain tensor
Setstrain=( 0.0 0.0 0.0 0.5 0.0 0.0 )

eval $(awk -v arr1="${Setstrain[*]}" '
           BEGIN{split(arr1, adds, " ");
	             adds[1]+=1;
			     adds[2]+=1;
			     adds[3]+=1;
			     adds[4]/=2;
			     adds[5]/=2;
			     adds[6]/=2}
		   NR==2 {scal=$1}
		   NR>=3 && NR<=5 {r0[NR]=$1*scal;s0[NR]=$2*scal;t0[NR]=$3*scal}
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

