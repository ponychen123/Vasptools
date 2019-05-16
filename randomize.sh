#!/bin/bash
#this script randomize the position in POSCAR, you can adjust the step size 
#usage: ./randomize.sh   only support POSCAR
#ponychen 20190516
#email:18709821294@outlook.com

step=0.05 #step size for randomization, 0.1A is surfficien, maybe smaller
end=`awk 'NR == 7 {print $1+$2+$3+$4+$5+8}' POSCAR`;

grep -in "Sel" POSCAR > /dev/null 2>&1 && constr=0 || constr=1;
grep -in "direct" POSCAR > /dev/null 2>&1 && dire=0 || dire=1;

if [ $dire -eq 0 ];then
	echo "Please transform your POSCAR into Cartesian !!!";
else
	if [ $constr -eq 1 ];then
		begin=9
	else
		begin=10
	fi
	hed=$(($begin-1))
	head -n $hed POSCAR | cat - > POSCAR.new
	awk -v begin=$begin -v end=$end -v step=$step '
	    NR>=begin && NR<=end {x0[NR]=$1;y0[NR]=$2;z0[NR]=$3;
		    u[NR]=$4;v[NR]=$5;w[NR]=$6}
		END{for(i=begin;i<=end;i++){
		        r=2*rand()-1;
				s=2*rand()-1;
				t=2*rand()-1;
				d=sqrt(r^2+s^2+t^2);
				r=r/d*step;
				s=s/d*step;
				t=t/d*step;
				x1[i]=x0[i]*(1+r);
				y1[i]=y0[i]*(1+s);
				z1[i]=z0[i]*(1+t);
				printf("%9.6f\t%9.6f\t%9.6f\t%s\t%s\t%s\t\n",x1[i],y1[i],z1[i],u[i],v[i],w[i])
			}}' POSCAR >> POSCAR.new
	fi

