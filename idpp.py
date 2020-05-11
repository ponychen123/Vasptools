#!/usr/bin/python3
#usage: a simple script to make linear path for NEB using IDPP(Hannes Jonsson 2014)
##################################################################################
#there are two way to run this script:
#1. just run ./idpp.py and follow the instrucment
#2. ./idpp.py numberofimages initialfilename finalfilename
#################################################################################
#first write by lipai@mail.ustc.edu.cn 2016/06/22
#ponychen write this in python3, add support for no frozen atom type ,fix
#bugs in three axiss of select frozen type, add NEB method
#2019/06/15
#email:18709821294@outlook.com
#20190820:add support for POSCAR in cartesian format, and add output for a XDACAR watching movie of
#transition path.All the output images are in cartesian format. by ponychen
#20190822:add support reading user designed initial transition path, it should be careful that all
#20200511:add support NEB method, and another mode for reding argv
#the images should in same format. by ponychen

import numpy as np
import os
import re
import sys

#some default values, you may change it depend on your condition
<<<<<<< HEAD
step_init = 0.01  # step size, unit, Angstrom
readfromexits = False # read transition path from user built, default False
conver = 0.1 # converge thershould
linear = False #just to do linear interpolation only
lneb = True # whther to open NEB
spring = 1 # spring constant, indead it means the fraction of spring force
scale = 3 # scale constant to decrease step size when crossing hollow
=======
step_init = 0.0001  #step size, small in the case od direct format
readfromexits = False #read transition path from user built, default False
conver = 10 #converge thershould
linear = False #just to do linear interpolation only
>>>>>>> refs/remotes/origin/master

#input filename and number of images
if sys.argv[1]:
    images = int(sys.argv[1])
    if not readfromexits:
        ininame = sys.argv[2]
        finname = sys.argv[3]
else:
    images = int(input("please input number of images: "))
    if not readfromexits:
        ininame = input("please input name of initial structure: ")
        finname = input("please input name of final structure: ")

#read initial structures
if readfromexits:
    fileopen = open("00/POSCAR", 'r') #read from the initial structure of exsting transition path
else:
    fileopen = open(ininame,'r')
ini_data = fileopen.readlines()
fileopen.close()

#check whether atom are being fronzen
if (re.search('sel', ini_data[7], re.I)):
    head = ini_data[:9]
    atom_num = sum(map(int, head[6].split()))
    ini_data = ini_data[9:9+atom_num]
    frozen = 1
    #check whether the coordination are in cartesian format or direct format
    if re.search('dir', head[8], re.I):
        direct = 1
        head[8] = "Cartesian\n"
    else:
        direct = 0
else:
    head = ini_data[:8]
    atom_num = sum(map(int, head[6].split()))
    ini_data = ini_data[8:8+atom_num]
    frozen = 0
    #check whether the coordination are in Cartesian format or direct format
    if re.search('dir', head[7], re.I):
        direct = 1
        head[7] = "Cartesian\n"
    else:
        direct = 0

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
if readfromexits:
    if images < 9:
        filename = "0"+str(images+1)+"/POSCAR"
    else:
        filename = str(images+1)+"/POSCAR"
    fileopen = open(filename, "r")
else:
    fileopen = open(finname, "r")
fin_data = fileopen.readlines()
fileopen.close()

#keep frozen condition same with initial structure
if (frozen == 1):
    fin_data = fin_data[9:9+atom_num]
else:
    fin_data = fin_data[8:8+atom_num]

tmp = []
for i in fin_data:
    tmp.append(list(map(float, i.split()[0:3])))
pos_b = np.array(tmp)

#if read from existing path, then read all the images to pos_im
if readfromexits:
    pos_im = np.zeros([images, atom_num, 3])
    for i in range(images):
        if i+1 < 10:
            filename = "0"+str(i+1)+"/POSCAR"
        else:
            filename = str(i+1)+"/POSCAR"
        fileopen = open(filename, "r")
        image_data = fileopen.readlines()
        fileopen.close()
        if frozen == 1:
            image_data = image_data[9:9+atom_num]
        else:
            image_data = image_data[8:8+atom_num]
        tmp = []
        for j in image_data:
            tmp.append(list(map(float, j.split()[0:3])))
        pos_im[i] = np.array(tmp)

#if input POSCARs are in cartesian format, then transfer them into direct format
#firstly read the coordination matrix of three bias axis, not support for the case of ssNEB
if direct:
    tmp = []
    for i in range(2,5):
        tmp.append(list(map(float, head[i].split())))
    axis = np.array(tmp)
    
    if readfromexits:
        for i in range(atom_num):
            pos_a[i] = np.dot(pos_a[i], axis)
            pos_b[i] = np.dot(pos_b[i], axis)
        for i in range(images):
            for j in range(atom_num):
                pos_im[i,j] = np.dot(pos_im[i,j], axis)
    else:
        for i in range(atom_num):
            pos_a[i] = np.dot(pos_a[i], axis)
            pos_b[i] = np.dot(pos_b[i], axis)

#get distance matrix and linear interpolation
dist_a = np.zeros([atom_num, atom_num])
dist_b = np.zeros([atom_num, atom_num])
for i in range(atom_num):
    for j in range(atom_num):
        dist_a[i,j] = np.linalg.norm(pos_a[i]-pos_b[j]) #distance between i and j atoms in pos_a
        dist_b[i,j] = np.linalg.norm(pos_a[i]-pos_b[j])

dist_im = np.zeros([images, atom_num, atom_num]) #3D distance matrix
if readfromexits:
    for  i in range(images):
        for j in range(atom_num):
            for k in range(atom_num):
                dist_im[i,j,k] = np.linalg.norm(pos_im[i,j]-pos_im[i,k])
else:
    pos_im = np.zeros([images, atom_num, 3])  #3D position matrix
    for i in range(images):           #linear interpolation
        dist_im[i] = dist_a+(i+1.0)*(dist_b-dist_a)/(images+1.0)
        pos_im[i] = pos_a+(i+1.0)*(pos_b-pos_a)/(images+1.0)
if not linear:
    #optimization using steepest descent method
    pos_tmp = np.zeros([atom_num, 3])
<<<<<<< HEAD
    dist_tmp = np.zeros([images, atom_num, atom_num])
    s0 = np.zeros(images)
    s1 = np.zeros(images)
    flag = 0
    loop = 0
    grad1 = np.zeros([images, atom_num, 3]) # atom force
    grad2 = np.zeros([images, atom_num, 3]) # spring force
    grad3 = np.zeros([images, atom_num, 3]) # true force
    nn = 0

    while(1):
        #decrease step size when crossing hollow
        if (nn > 0):
            step = step_init/scale
        else:
            step = step_init
        #atomic force
        for im in range(images):
            for i in range(atom_num): #get the distant matrix for each image
                for j in range(atom_num):
                    if (i == j):
                        dist_tmp[im,i, j] = 10
                    else:
                        tmp = np.linalg.norm(pos_im[im,i]-pos_im[im,j])
                        dist_tmp[im,i,j] = tmp

            for i in range(atom_num):
                for sigma in range(3):
                    gradtmp1 = 0
                    if (frozen == 1 and fix[i,sigma] == "T") or frozen == 0:
                        for j in range(atom_num): #get the partial differencial of Sidpp
                            if (j != i): #get the vitual force
                                gradtmp1 += 2*(dist_im[im,i,j]-dist_tmp[im,i,j])*(pos_im[im,i,sigma]-pos_im[im,j,sigma])\
                                        *(2*dist_im[im,i,j]-dist_tmp[im,i,j])/dist_tmp[im,i, j]**5
                    grad1[im,i,sigma] = gradtmp1
                # unify the vitual force
                length = np.linalg.norm(grad1[im,i])
                grad1[im,i] /= length
        
        # spring force
        if lneb:
            tau1 = np.zeros([atom_num,3])
            tau2 = np.zeros([atom_num,3])
            tau = np.zeros([images, atom_num, 3])
            #normalized line segment
            for im in range(images):
                if (im == 0):
                    tmp1 = pos_im[im] - pos_a
                    tmp2 = pos_im[im+1] - pos_im[im]
                    for i in range(atom_num):
                        tau1[i] = tmp1[i]/np.linalg.norm(tmp1[i])
                        tau2[i] = tmp2[i]/np.linalg.norm(tmp2[i])
                        tau[im,i] = (tau1[i]+tau2[i])/np.linalg.norm(tau1[i]+tau2[i])
                elif (im == images - 1):
                    tmp1 = pos_im[im] - pos_im[im-1]
                    tmp2 = pos_b - pos_im[im]
                    for i in range(atom_num):
                        tau1[i] = tmp1[i]/np.linalg.norm(tmp1[i])
                        tau2[i] = tmp2[i]/np.linalg.norm(tmp2[i])
                        tau[im,i] = (tau1[i]+tau2[i])/np.linalg.norm(tau1[i]+tau2[i])
                else: 
                    tmp1 = pos_im[im] - pos_im[im-1]
                    tmp2 = pos_im[im+1] - pos_im[im]
                    for i in range(atom_num):
                        tau1[i] = tmp1[i]/np.linalg.norm(tmp1[i])
                        tau2[i] = tmp2[i]/np.linalg.norm(tmp2[i])
                        tau[im,i] = (tau1[i]+tau2[i])/np.linalg.norm(tau1[i]+tau2[i])
            # spring force
            for im in range(images):
                if (im == 0):
                    for i in range(atom_num):
                        tmp = pos_im[im+1,i]-pos_a[i]
                        grad2[im,i] = np.dot(tmp,tau[im,i])*tau[im,i]
                elif (im == images-1):
                    for i in range(atom_num):
                        tmp = pos_b[i]-pos_im[im-1,i]
                        grad2[im,i] = np.dot(tmp,tau[im,i])*tau[im,i]
                else:
                    for i in range(atom_num):
                        tmp = pos_im[im+1,i]-pos_im[im-1,i]
                        grad2[im,i] = np.dot(tmp,tau[im,i])*tau[im,i]
            # true force perpendicular to the tangent
            for im in range(images):
                for i in range(atom_num):
                    tmp = np.dot(grad1[im,i],tau[im,i])*tau[im,i]
                    grad3[im,i] = tmp
            grad3 = grad1-grad3
        # update atomic coordination
        if lneb:
            grad = spring*grad2+grad3
        else:
            grad = grad1
        
        if frozen == 0:
            pos_im += step*grad
        else:
            #do not relax fixed ions
            for i in range(atom_num):
                for j in range(3):
                    if (fix[i,j] == "T"):
                        grad[:,i,j] = 0.0
            pos_im += step*grad

        #judge convergence
        loop += 1
        print("loop: "+str(loop))
    
        flag = 0 #reset converge flag
        nn = 0
        for im in range(images):
=======
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
>>>>>>> refs/remotes/origin/master
            s0[im] = s1[im]
            s1[im] = 0
            for i in range(atom_num):
                for j in range(i):
<<<<<<< HEAD
                    s1[im] += (dist_im[im,i,j]-dist_tmp[im,i,j])**2/dist_tmp[im,i,j]**4
            if (abs(s0[im]-s1[im]) < conver):
                print("loop: "+str(loop)+" image "+ str(im+1) +" converge!!!")
                flag += 1
            if (loop > 0 and s1[im] > s0[im]):
                nn += 1

        if (flag == images):
            print("Hey! IDPP path has converged!!!")
            break
        
=======
                    s1[im] += (dist_im[im][i][j]-dist_tmp[i][j])**2/dist_tmp[i][j]**4
            loop += 1
            print("loop: " + str(loop))
            if (abs(s0[im]-s1[im]) < conver):
                print("image "+ str(im+1) +" converge!!!")
                flag[im] = 1
                break
            if (loop > 1 and s1[im] > s0[im]): #decrease step size when cross hollow
                step = step/3
>>>>>>> refs/remotes/origin/master

#mkdir and generate poscar file for neb
if (images + 1 < 10):
    num = '0' + str(images +1)
else:
    num = str(images +1)

if not readfromexits:
    os.system("mkdir -p 00")
    f = open("00/POSCAR", "w")
    f.writelines(head)
    data = pos_a.tolist()
    for i in range(atom_num):
        line = map(str, data[i])
        line = " ".join(line)
        if frozen == 1:
            line = line + "    " +fix[i,0]+fix[i,1]+fix[i,2]+"\n"
        else:
            line += "\n"
        f.write(line)
    f.close()
    os.system("mkdir -p "+num)
    filename = str(num)+"/POSCAR"
    f = open(filename, "w")
    f.writelines(head)
    data = pos_b.tolist()
    for i in range(atom_num):
        line = map(str, data[i])
        line = " ".join(line)
        if frozen == 1:
            line = line + "    " +fix[i,0]+fix[i,1]+fix[i,2]+"\n"
        else:
            line += "\n"
        f.write(line)
    f.close()
    for i in range(images):
        if i+1<10:
            num = "0"+str(i+1)
        else:
            num = str(i+1)
        os.system("mkdir -p "+num)
        data = pos_im[i].tolist()
        filename = num + "/POSCAR"
        f = open(filename, "w")
        f.writelines(head)
        for j in range(atom_num):
            line = map(str, data[j])
            line = " ".join(line)
            if frozen == 1:
                line = line + "    " + fix[j,0] +fix[j,1] +fix[j,2] + "\n"
            else:
                line += "\n"
            f.write(line)
        f.close()
else:
    os.system("mkdir new")
    os.system("mkdir new/00")
    os.system("cp 00/POSCAR  new/00/POSCAR ")
    os.system("mkdir new/"+num)
    os.system("cp "+num+"/POSCAR"+" new/"+num+"/POSCAR")
    for i in range(images):
        if i+1<10:
            num = "new/"+"0"+str(i+1)
        else:
            num = "new/"+str(i+1)
        os.system("mkdir "+num)
        data = pos_im[i].tolist()
        filename = num + "/POSCAR"
        f = open(filename, "a+")
        f.writelines(head)
        for j in range(atom_num):
            line = map(str, data[j])
            line = " ".join(line)
            if frozen == 1:
                line = line + "    " + fix[j,0] +fix[j,1] +fix[j,2] + "\n"
            else:
                line += "\n"
            f.write(line)
        f.close()

#generate a XDATCAR watching the movie
f = open("XDATCAR", "w")
f.writelines(head[:7])
f.write("Cartesian configuration=     1\n")
for i in range(atom_num):
    line = map(str, pos_a[i])
    line = " ".join(line)
    line += "\n"
    f.write(line)
for i in range(images):
    f.write("Cartesian configuration=     "+str(i+2)+"\n")
    for j in range(atom_num):
        line = map(str, pos_im[i,j])
        line = " ".join(line)
        line += "\n"
        f.write(line)
f.write("Cartesian configuration=     "+str(images+2)+"\n")
for i in range(atom_num):
    line = map(str, pos_b[i])
    line = " ".join(line)
    line += "\n"
    f.write(line)
f.close()
