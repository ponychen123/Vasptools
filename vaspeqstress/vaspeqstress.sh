#!/bin/bash
#Instructions:
#1.Input relative parameters in the script (see below)
#2.remeber setting ISIF=2 in the INCAR
#3.put the script in the working directory, run it
#
#General algorithm:
#1.run vasp with ISIF=2 to get current pressure tensor
#2.using generalized Hooke's law and eleastic modulus you input, generating 
#  new POSCAR
#3.repeate 1-2 unitl the current pressure - target pressure in within convegency
#
#this script was first write by xiaohoufeichuan (sorry i cannot write Chinese
#here) with python2, pony chen update this code with python3 and polish the
#code. Futher pony chen will develop this code in bash 
#author: ponychen  email:18709821294@outlook.com
#2019/04/18
#2019/05/13 now ponychen write this in bash shell
#2019/05/17 add PBS columns
##PBS -L nodes=6:ppn=12
##PBS -L walltime=3:00:00
##PBS -V
##PBS -N test
#Set the target pressure tensor(xx, yy, zz, xy, yz, zx) in Kbar. Note pressure
#= -stress in VASP
Setpress=( 0.0 100.0 0.0 0.0 0.0 0.0 )
#convergency criteria for pressure,
presscirt=0.1
#Young's modulus in Kbar
E=2010
#Possion ratio
v=0.29
#Shear modulus in Kbar
G=$(echo "scale=8;$E/(2+2*$v)"|bc -l)
#maximum iteration cycles
imax=100
#damping factor for deformation
P=0.9
#Set your path effectively to run vasp
mpiexe="vasp_std"
#mpiexe='mpirun -np 12 vasp_std"

###do not change following codes unless you know what you are doing###
#######initial run to get the pressure tensor present####
cp POSCAR poscar.0
echo "starting calculation" > pressure.all
$mpiexe
cat OSZICAR > oszicar.all
cat OUTCAR > outcar.all

itr=0
num=`expr $imax - 1`
while [ $itr -lt $imax ]; do
	itr=`expr $itr + 1`
    echo "iteration = " $itr
	####scan whether the previous calculation is done######
	noempty=`grep "Total CPU time used (sec):" OUTCAR | tail -1`
    if [ -n "$noempty" ]; then
		pressstr=`grep "in kB" OUTCAR | tail -1 | awk '{print $3, $4,$5,$6,$7,$8}'`
		echo $pressstr >> pressure.all
		##calculate additional pressure tensor needed to achieve Setpress##
		eval $(awk -v arr1="${Setpress[*]}" -v arr2="${pressstr[*]}" -v cirt=$presscirt '
		          BEGIN{split(arr1,setp," ");split(arr2,nowp," ");
			      for(i=0;i<=5;i++){addp[i]=strtonum(setp[i+1])-strtonum(nowp[i+1])};
				      printf("addpress[0]=%f;addpress[1]=%f;addpress[2]=%f;\
				          addpress[3]=%f;addpress[4]=%f;addpress[5]=%f;",\
					      addp[0],addp[1],addp[2],addp[3],addp[4],addp[5])
				  max_p=0; 
				  for(j=0;j<=5;j++){
					  if( addp[j] < 0 ){addp[j]=-addp[j]};
					  if( addp[j] > max_p ){max_p=addp[j]}};
						  if( max_p > cirt ){goon=1}
						  else{goon=0};
				  printf("goon=%s;",goon)}
				  ')
	###stop if reaching the target pressure####
	if [ $goon == 0 ]; then
		echo "aborting calculations for the convergency was reaching!";
		break;
	else
		#add the strain needed to the new POSCAR
        cp CONTCAR POSCAR;
        eval $(awk -v arr1="${addpress[*]}" -v E=$E -v v=$v -v G=$G -v P=$P '
		    BEGIN{split(arr1, addp, " ");
			addM[0]=(addp[1]-v*(addp[2]+addp[3]))/(-E)*P+1;
			addM[1]=(addp[2]-v*(addp[1]+addp[3]))/(-E)*P+1;
			addM[2]=(addp[3]-v*(addp[2]+addp[1]))/(-E)*P+1;
            addM[3]=addp[4]/(-G)/2*P;
			addM[4]=addp[5]/(-G)/2*P;
			addM[5]=addp[6]/(-G)/2*P;}
            NR>=3 && NR <= 5 {r[NR-3]=$1;s[NR-3]=$2;t[NR-3]=$3;}
		END{for(i=0;i<=2;i++){
		    x[i]=r[i]*addM[0]+s[i]*addM[3]+t[i]*addM[5];
			y[i]=r[i]*addM[3]+s[i]*addM[1]+t[i]*addM[4];
			z[i]=r[i]*addM[5]+s[i]*addM[4]+t[i]*addM[2];
			printf("newM[%s]=\" %9.6f %9.6f %9.6f\n\";",i,x[i],y[i],z[i])}
            }
			' POSCAR);
         sed -i "3c${newM[0]}" POSCAR
		 sed -i "4c${newM[1]}" POSCAR
		 sed -i "5c${newM[2]}" POSCAR
		 cp POSCAR poscar.$itr
		 #rerun vasp with new POSCAR
         $mpiexe
		 cat OSZICAR >> oszicar.all
		 if [ $itr == $num ]; then
			 echo "reaching maximum steps, you may turn up the damping factor or lower the presscirt"
		 fi
	fi

fi
done

