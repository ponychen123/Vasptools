#!/usr/bin/python3
#this script sort the atom coordinations along specific crystal axis
#ponychen
#20190620
#email:18709821294@outlook.com

import re

#read in the enssential parameters
axis = int(input("please input which crystal axis do you want to sort :(1; (2; (3: "))
reverse = int(input("sort in ascending order(0) or descending order(1): "))

#open POSCAR and read in the atom coordinations and atom constraints information
fileopen = open("POSCAR",'r')
ini_data = fileopen.readlines()
fileopen.close()

if re.search('sel', ini_data[7], re.I):
    head = ini_data[:9]
    species_num = list(map(int, head[6].split()))
    atom_num = sum(species_num)
    ini_data = ini_data[9:9+atom_num]
    frozen = 1
    fomer = -1

else:
    head = ini_data[:8]
    species_num = list(map(int, head[6].split()))
    atom_num = sum(species_num)
    ini_data = ini_data[8:8+atom_num]
    frozen = 0

species = []
for i in range(len(species_num)):
    struc = []
    for j in range(species_num[i]):
        fomer += 1
        struc.append(ini_data[fomer])
            
    species.append(struc)    

for i in range(len(species_num)):
    tmp = []
    for j in range(species_num[i]):
        tmp = list(map(float, species[i][j].split()[0:3]))
        if frozen == 1:
            for k in range(3):
                tmp.append(species[i][j].split()[k+3])
        species[i][j] = tuple(tmp)

#sort the atoms
for i in range(len(species_num)):
    if axis == 3:
        if reverse == 0:
            species[i].sort(key=lambda x: x[2])
        else:
            species[i].sort(key=lambda x: x[2], reverse=True)
    elif axis == 2:
        if reverse == 0:
            species[i].sort(key=lambda x: x[1])
        else:
            species[i].sort(key=lambda x: x[1], reverse=True)
    elif axis == 1:
        if reverse:
            species[i].sort(key=lambda x: x[0])
        else:
            species[i].sort(key=lambda x: x[0], reverse=True)
    else:
        print("please choose axis direction in 1 2 or 3 !!!")

#write the sorted atoms coordination to POSCAR.sort
outname = "POSCAR.sort"
f = open(outname, "w+")
f.writelines(head)
for i in range(len(species_num)):
    for j in range(species_num[i]):
        line = map(str, species[i][j])
        line = " ".join(line)
        line += "\n"
        f.write(line)
f.close()
