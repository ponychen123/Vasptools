#!/usr/bin/python3
#usage: a simple script to make linear path for NEB using IDPP(Hannes Jonsson 2014)
#       just run 
#first write by lipai@mail.ustc.edu.cn 2016/06/22
#ponychen write this in python3, add support for no frozen atom type ,fix
#bugs in three axiss of select frozen type
#2019/06/15
#email:18709821294@outlook.com

import numpy as np
import os
import re

step_init = 0.0001
images = int(input("please input number of images: "))
ininame = input("please input name of initial structure: ")
finname = input("please input name of final structure: ")

#read initial structures
fileopen = open(ininame,'r')
ini_data = fileopen.readlines()

#check whether atom are being fronzen
if (re.search('sel', ini_data[7], re.I)):
    head = ini_data[:9]
    atom_num = sum(map(int, head[6].split()))
    ini_data = ini_data[9:9+atom_num]
    frozen = 1
else:
    head = ini_data[:8]
    atom_num = sum(map(int, head[6].split()))
    ini_data = ini_data[8:8+atom_num]
    frozen = 0

tmp = []
fix = []
fixx = [0, 0, 0]
for i in range(atom_num):
    tmp.append(list(map(float, ini_data[i].split()[0:3])))
    if ( frozen == 1):
        for j in range(3):
            if(ini_data[i].split()[j+3] == 'F'):
                fixx[j] = "F"
            else:
                fixx[j] = "T"
        fix.append([fixx[0], fixx[1], fixx[2]])
pos_a = np.array(tmp)

#read final structure
fileopen = open(finname, 'r')
fin_data = fileopen.readlines()

#keep frozen condition same with initial structure
if (frozen == 1):
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
        if(pos_a[i][j]-pos_b[i][j] > 0.5):
            pos_a[i][j] -= 1
            if(pos_a[i][j]-pos_b[i][j] < -0.5):
                pos_b[i][j] -= 1

#get distance matrix and linear interpolation
dist_a = np.zeros([atom_num, atom_num])
dist_b = np.zeros([atom_num, atom_num])
for i in range(atom_num):
    for j in range(atom_num):
        tmp_a = 0
        tmp_b = 0
        for k in range(3):
            tmp_a += (pos_a[i][k]-pos_a[j][k])**2
            tmp_b += (pos_b[i][k]-pos_b[j][k])**2
        dist_a[i,j] = np.sqrt(tmp_a) #distance between i and j atoms in pos_a
        dist_b[i,j] = np.sqrt(tmp_b)

dist_im = np.zeros([images, atom_num, atom_num]) #3D distance matrix
pos_im = np.zeros([images, atom_num, 3])  #3D position matrix

for i in range(images):           #linear interpolation
    dist_im[i] = dist_a+(i+1.0)*(dist_b-dist_a)/(images+1.0)
    pos_im[i] = pos_a+(i+1.0)*(pos_b-pos_a)/(images+1.0)

#optimization using steepest descent method
pos_tmp = np.zeros([atom_num, 3])
dist_tmp = np.zeros([atom_num, atom_num])
s0 = np.zeros(images)
s1 = np.zeros(images)
flag = np.zeros(images)

for im in range(images):
    if (flag[im] == 1): #avoid repetition
        continue
    step = step_init
    print("generate image " + str(im+1))
    loop = 0
    while(1):
        for i in range(atom_num): #get the distant matrix for each image
            for j in range(atom_num):
                if (i == j):
                    dist_tmp[i, j] = 10
                else:
                    tmp = 0
                    for k in range(3):
                        tmp += (pos_im[im][i][k]-pos_im[im][j][k])**2
                    dist_tmp[i,j] = np.sqrt(tmp)

        for i in range(atom_num):
            for sigma in range(3):
                grad = 0
                if (frozen == 1 and fix[i][sigma] == "T") or frozen == 0:
                    for j in range(atom_num): #get the partial differencial of Sidpp
                        if (j != i): #get the vitual force
                            grad += 2*(dist_im[im][i][j]-dist_tmp[i][j])*(pos_im[im][i][sigma]-pos_im[im][j][sigma])\
                                    *(2*dist_im[im][i][j]-dist_tmp[i][j])/dist_tmp[i, j]**5
                pos_tmp[i][sigma] = pos_im[im][i][sigma] + step*grad
        pos_im[im] = pos_tmp

        #judge convergence
        s0[im] = s1[im]
        s1[im] = 0
        for i in range(atom_num):
            for j in range(i):
                s1[im] += (dist_im[im][i][j]-dist_tmp[i][j])**2/dist_tmp[i][j]**4
        loop += 1
        print("loop: " + str(loop))
        if (abs(s0[im]-s1[im]) < 0.01):
            print("image "+ str(im+1) +" converge!!!")
            flag[im] = 1
            break
        if (loop > 1 and s1[im] > s0[im]): #decrease step size when cross hollow
            step = step/3

#mkdir and generate poscar file for neb
if (images + 1 < 10):
    num = '0' + str(images +1)
else:
    num = str(images +1)
os.system(" mkdir 00 ")
os.system(" cp " + ininame + " 00/POSCAR ")
os.system(" mkdir " + num)
os.system(" cp  " + finname + " " + num + "/POSCAR")
for i in range(images):
    if (i + 1 < 10):
        num = '0' + str(i+1)
    else:
        num = str(i+1)
    os.system("mkdir " + num)
    data = pos_im[i].tolist()
    filename = num + "/POSCAR"
    f = open(filename, "a+")
    f.writelines(head)
    for j in range(atom_num):
        line = map(str, data[j])
        line = " ".join(line)
        if frozen == 1:
            line = line + "    " + str(fix[j][0]) + "    " + str(fix[j][1]) + "    " + str(fix[j][2]) + "\n"       
        else:
            line = line + "\n"
        f.write(line)
    f.close()

