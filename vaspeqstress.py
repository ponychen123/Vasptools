#!/usr/bin/python3
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

import os
import subprocess
import shutil
import numpy as np
import linecache

##############you need to configure the following variable#################
#Set the target pressure tensor(xx, xy, zz, xy, yz, zx) in kB. Note pressure
#-stress in VASP
Setpress = np.array([60, 0, 0, 0, 0, 0])
#convergency criteria for pressure, unit in kB
presscirt = 0.1
#Young's modulus in Kb
E = 2010
#Possion ratio
v = 0.29
#Shear modulus in kB
G = E / ( 2 + 2*v)
#maximum iteration cycles
imax = 30
#damping factor for deformation
P = 0.9
#set your path effective to run vasp
mpiexe = 'vasp_std'

###do not change following codes unless you know what you are doing####
#######initial run to get the pressure tensor present#######
shutil.copy('POSCAR', 'poscar.0')
subprocess.run("echo 'starting calculation' > pressure.all", shell=True)
subprocess.run(mpiexe, shell=True)

subprocess.run("cat OSZICAR > oszicar.all", shell=True)
subprocess.run("cat OUTCAR > outcar.all", shell=True)
subprocess.run("grep 'Total CPU time used (sec):' OUTCAR | tail -n 1 > end.txt", 
        shell=True)

itr = 0
while itr < imax:
    itr += 1
    print("iteration = ", itr)
######scan whether the previous calculation is complicated########    
    subprocess.run("grep 'Total CPU time used (sec):' OUTCAR | \
        tail -n 1 > end.txt", shell = True)
    noempty = os.path.getsize('end.txt')
    os.remove('end.txt')

    if noempty:
######read the pressure tensor from OUTCAR #######################        
        subprocess.run("grep 'in kB' OUTCAR | tail -n 1 > pressure.txt",
                 shell=True)
        subprocess.run("cat pressure.txt >> pressure.all", shell=True)
        for line in open("pressure.txt"):
            press = np.array([float(x) for x in line.split()[2:]])
            print("pressure present = ", press)
        os.remove("pressure.txt")
#########calculate the additional pressure tensor meede to achieve Setpress
        addpress = Setpress - press
        print("adding pressure = ", addpress)

        abs_P = max([abs(y) for y in addpress])
#########stop if the additional pressure is small enough#############        
        if abs_P < presscirt:
            subprocess.run("echo 'aborting calculation for thr convergence is reached' >> pressure.all",
                    shell=True)
            break
        else:
            shutil.copy('CONTCAR', 'POSCAR')
########read current POSCAR matrix            
            pos = linecache.getlines('POSCAR')
            M = np.zeros((3, 3))
            addM = np.zeros((3, 3))
            M[0, :] = np.array([float(x) for x in pos[2].split()])
            M[1, :] = np.array([float(x) for x in pos[3].split()])
            M[2, :] = np.array([float(x) for x in pos[4].split()])
#calculate the strain we nedd to achieve addpress according to general
#Hook's law.Beaware that VASP defines compressive pressure as positive,
#hence /E becomes /(-E)
            addM[0, 0] = (addpress[0] - v*(addpress[1] + addpress[2]))/(-E)
            addM[1, 1] = (addpress[1] - v*(addpress[0] + addpress[2]))/(-E)
            addM[2, 2] = (addpress[2] - v*(addpress[1] + addpress[0]))/(-E)
            addM[0, 1] = addM[1, 0] = addpress[3]/(-G)/2
            addM[1, 2] = addM[2, 1] = addpress[4]/(-G)/2
            addM[0, 2] = addM[2, 0] = addpress[5]/(-G)/2

            addM = np.diag([1, 1, 1]) + addM*P
            M = np.dot(M, addM)
            print("adjusting matrix to : ", M)
            np.savetxt("M.txt", M)
            linecache.clearcache()
            for i in range(3):
                pos[i+2] = str(M[i,0]) + "    " + str(M[i,1]) + "    " \
                        + str(M[i,2]) + "\n"
            fo = open("POSCAR", "w+")
            fo.writelines(pos)
            fo.close() 
            posname = 'poscar.' + str(itr)
            shutil.copy('POSCAR', posname)
####run vasp to calculate pressure for this POSCAR
            subprocess.run(mpiexe, shell=True)
            subprocess.run("cat OSZICAR >> oszicar.all", shell=True)
   
    else:
        print("calculation aborted")



