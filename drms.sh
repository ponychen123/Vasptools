#!/bin/bash
#this script tells you the root mean squre bonds length
#Usage: dist.sh inputfile
#2019/04/25 pony chen
#email:cjchen16s@imr.ac.cn

end=`awk 'NR == 7 {print $1+$2+$3+$4+$5+8}' $1`
#end    
r0=2.884
#r0 is the cutoff for 1NN, you shold set a proper value
grep -in "Sel" $1 > /dev/null 2>&1 && constr=0 || constr=1;
grep -in "direct" $1 > /dev/null 2>&1 && dire=0 || dire=1;

if [ $constr != 0 ]; then
		begin=9
else
		begin=10
	end=`expr $end + 1`
fi

if [ $dire -eq 0 ];then


	     awk -v begin=$begin -v end=$end -v r0=$r0 '
		   NR>=3 && NR<=5 {r[NR]=$1;s[NR]=$2;t[NR]=$3}
		   NR>=begin && NR<=end {x0[NR]=$1;y0[NR]=$2;z0[NR]=$3}
	       END{for(i=begin;i<=end;i=i+1){
		        x1[i]=x0[i]*r[3] + y0[i]*r[4] + z0[i]*r[5];
				y1[i]=x0[i]*s[3] + y0[i]*s[4] + z0[i]*s[5];
				z1[i]=x0[i]*t[3] + y0[i]*t[4] + z0[i]*t[5];
		    	};
				u=0;
				v=0;
               for(j=begin;j<end;j=j+1){
				   for(k=j+1;k<=end;k=k+1){
					   m=(x1[k]-x1[j])^2 + (y1[k]-y1[j])^2 + (z1[k]-z1[j])^2;
					   if (sqrt(m) < r0 ){
                           u=u + m;
						   v=v + 1;
                       };				
				   };
			   };
			   q=sqrt(u/v);
			   printf("the rms bond distance is %f \n", q);
		}' $1 
else
	     awk -v begin=$begin -v end=$end -v r0=$r0 '
		         NR>=begin && NR<=end {x1[NR]=$1;y1[NR]=$2;z1[NR]=$3}
			     END{u=0;
				     v=0;
				     for(j=begin;j<end;j=j+1){
				         for(k=j+1;k<=end;k=k+1){
					        m=(x1[k]-x1[j])^2 + (y1[k]-y1[j])^2 + (z1[k]-z1[j])^2;
							if (sqrt(m) < r0 ){
						          u=u + m ;
								  v=v + 1 ;
							  };
				         };
			            };
						q=sqrt(u/v);
						printf("the rms bond distance is %f \n", q);
				 }' $1 ;
fi



