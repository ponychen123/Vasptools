#!/usr/bin/python3
#minimum atom displacement vector (madv) method by ponychen
#just use!
#author:ponychen
#20190819
#email:18709821294@outlook.com

import numpy as np
import os
import re
import sys

#some default values, you may change it depend on you condition
step_init = 0.3 #step size
dmin = 2.3  #below which the atoms are thought to be close, unit angstrom
maxitr = 100 #the maximum iteration for each cycle in optimization
amp = 2 #the amptitude of weight function A/d^-N
power = 4 #the power of weight function A/d^-N
uniform = False #experimental use

#read from user input
images = int(input("please input number of images: "))
ininame = input("please input name of initial structure: ")
finname = input("please input name of final structure: ")

#read initial structures
fileopen = open(ininame,'r')
ini_data = fileopen.readlines()
fileopen.close()

#check whether atom are being frozen
if (re.search('sel', ini_data[7], re.I)):
    head = ini_data[:9]
    atom_num = sum(map(int, head[6].split()))
    ini_data = ini_data[9:9+atom_num]
    frozen = 1
    head[8] = "Cartesian \n"
else:
    head = ini_data[:8]
    atom_num = sum(map(int, head[6].split()))
    ini_data = ini_data[8:8+atom_num]
    frozen = 0
    head[7] = "Cartesian \n"

tmp = []
fix = []
fixx = [0, 0, 0]
for i in range(atom_num):
    tmp.append(list(map(float, ini_data[i].split()[0:3])))
    if frozen == 1:
        for i in range(3):
            if ini_data[i].split()[j+3] == "F":
                fixx[j] = "F"
            else:
                fixx[j] = "T"
        fix.append([fixx[0], fixx[1], fixx[2]])

pos_a = np.array(tmp)

#read the coordition matrix of three bias axis, not support for the case of ssNEB
tmp = []
for i in range(2,5):
    tmp.append(list(map(float, head[i].split())))
axis = np.array(tmp)

#read final structure
fileopen = open(finname, "r")
fin_data = fileopen.readlines()

#keep frozen condition same with initial structure
if frozen == 1:
    fin_data = fin_data[9:9+atom_num]
else:
    fin_data = fin_data[8:8+atom_num]

tmp = []
for i in fin_data:
    tmp.append(list(map(float, i.split()[0:3])))
pos_b = np.array(tmp)

#correction of periodic boundary condition only support direct format
for i in range(atom_num):
    for j in range(3):
        if pos_a[i,j] - pos_b[i,j] > 0.5:
            pos_a[i,j] -= 1
        if pos_a[i,j] - pos_b[i,j] < -0.5:
            pos_b[i,j] -= 1

#get linear interpolation between initial and final structure
pos_im = np.zeros([images, atom_num, 3]) #3D position matrix
for i in range(images):
    pos_im[i] = pos_a+(i+1)*(pos_b-pos_a)/(images+1.0)

#transfer the direct coordination to cartesian coordination
for i in range(images):
    for j in range(atom_num):
        pos_im[i,j] = np.dot(pos_im[i,j],axis)

for i in range(atom_num):
    pos_a[i] = np.dot(pos_a[i],axis)
    pos_b[i] = np.dot(pos_b[i],axis)

#optimize the atoms that are too close based on hard sphere model
for i in range(images):
    flag = 1 #start the following loop
    itr = 0 #initilize the step number
    while flag:
        #initialize
        flag = 0 #agin initialize, this means this program default was there no atoms being closely
        advec = np.zeros([atom_num,3]) #a matrix to store the displacement vector of each atom
        itr += 1
        print("now calculating the "+str(i+1)+" image, iteration "+str(i))

        for j in range(atom_num):
            for k in range(atom_num):
                if j != k:
                    tmp = pos_im[i,j]-pos_im[i,k]
                    ds = np.sqrt(sum(tmp**2))
                    #check whther this two atom meet too cloaely
                    if ds < dmin:
                        #apply weighting function
                        flag += 1
                        advec[j] += amp/ds**power*tmp
        #displace atoms by displacement vector
        if frozen == 1:
            for j in range(atom_num):
                for k in range(3):
                    if fix[j,k] == "T":
                        pos_im[i,j,k] += step_init*advec[j.k]
        else:
            for j in range(atom_num):
                pos_im[i,j] += step_init*advec[j]

        if flag == 0:
            print("the "+str(i+1)+" image has converged!")

        #if iteration steps reaching maxitr, break and you should check relative parameters
        if itr > maxitr:
            sys.exit("buddy, please check default parameters, maybe you should alter them, or you can try idpp.py write by me.")


#mkdir and generate poscar file for neb
if images + 1 < 10:
    num = "0" + str(images+1)
else:
    num = str(images+1)
os.system("mkdir 00")
os.system("cp "+ininame+" 00/POSCAR ")
os.system("mkdir "+num)
os.system("cp "+finname+" "+num+"/POSCAR")
for i in range(images):
    if i+1<10:
        num = "0"+str(i+1)
    else:
        num = str(i+1)
    os.system("mkdir "+num)
    data = pos_im[i].tolist()
    filename = num + "/POSCAR"
    f = open(filename, "a+")
    f.writelines(head)
    for j in range(atom_num):
        line = map(str, data[j])
        line = " ".join(line)
        if frozen == 1:
            line = line + "    " + fix[j][0] +fix[j][1] +fix[j][2] + "\n"
        else:
            line += "\n"
        f.write(line)
    f.close()
