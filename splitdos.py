#!/usr/bin/python3
#this script split dos into specific format, basely need DOSCAR and sometimes need POSCAR
#just run this script and it will tell you what to choose
#ponychen
#20190624
#email:18709821294@outlook.com

#read DOSCAR only need DOSCAR
fileopen = open("DOSCAR", 'r')
ini_data = fileopen.readlines()
fileopen.close()

#get some base information of DOS
dos_start = float(ini_data[5].split()[0]) #start energy of DOS mesh
dos_end = float(ini_data[5].split()[1])   #end energy of DOS mesh
dos_num = int(ini_data[5].split()[2])    # number of points in DOS mesh
fermi_level = float(ini_data[5].split()[3]) 
atom_num = int(ini_data[0].split()[0])   #total atom numbers

def totalDOS():
    #get the total DOS of system, both suitable for no magnetic and no colineair
    totalDOS = []
    for i in range(dos_num):
        tmp = list(map(float, ini_data[6+i].split()[:2]))
        tmp[0] = round(tmp[0]-fermi_level, 4)
        tmp[1] = round(tmp[1], 4)
        totalDOS.append(tmp)
    
    #write total DOS to total_dos.dat
    outname = "total_dos.dat"
    f = open(outname, "w+")
    f.write("#energy   totalDOS\n")  
    for i in range(dos_num):
        line = map(str,totalDOS[i])
        line = "    ".join(line)
        line += "\n"
        f.write(line)
    f.close()

def LDOS(ini_data=ini_data):
    #get the atom projected dos 
    LDOS = []

    def specific_ldos(atom_no, ini_data=ini_data, flag=0):
        #get LDOS for one specifc atom
        #get relative block in DOSCAR for specific atom
        atom_no = int(atom_no)
        dos_begin = 5 + (atom_no+1)*dos_num + atom_no + 2
        ini_data = ini_data[dos_begin:dos_begin+dos_num] #be aware of the section
        if flag == 0: #get total LDOS for one specific atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                counter = 0
                for j in tmp[1:]:
                    counter += j
                counter = round(counter, 4)
                LDOS.append([tmp[0], counter])
        
            #write ldos of specific atom to output file
            outname = "atom" + str(atom_no) + "ldos.dat"
            f = open(outname, "w+")
            f.write("#energy    totalDOS\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                f.write(line)
            f.close()

        elif flag == 1: #get shell projected DOS for one specific atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                shell_type = len(tmp)-1
                if shell_type == 1:
                    s = round(tmp[1], 4)
                    LDOS.append([tmp[0], s])
                elif shell_type == 4:
                    s = round(tmp[1], 4)
                    p = round(sum(tmp[2:5]), 4)
                    LDOS.append([tmp[0], s, p])
                elif shell_type == 9:
                    s = round(tmp[1], 4)
                    p = round(sum(tmp[2:5]), 4)
                    d = round(sum(tmp[5:10]), 4)
                    LDOS.append([tmp[0], s, p, d])
                elif shell_type == 16:
                    s = round(tmp[1], 4)
                    p = round(sum(tmp[2:5]), 4)
                    d = round(sum(tmp[5:10]), 4)
                    f = round(sum(tmp[10:17]), 4)
                    LDOS.append([tmp[0], s, p, d, f])
                else:
                    print("shell nums out of range! impossible!!!")

            outname = "atom" + str(atom_no) + "ldos.dat"
            ff = open(outname, "w+")
            ff.write("#energy   s   p   d   f\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()    
        elif flag == 2: #get orbital projected LDOS for one specifc atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                for j in tmp[1:]:
                    j = round(j, 4)
                LDOS.append(tmp)

            outname = "atom" + str(atom_no) + "ldos.dat"
            ff = open(outname, "w+")
            ff.write("#energy   s py pz px dxy dyz dz2 dxz dx2\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()
        else:
            print("write set your flag!!!!")

    flag = int(input("choose type of atom projection: 0)total 1)shell 2)orbital:  "))
    if flag == 0:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n)
                else:
                    print("your input atom number was out of range!!check!!")
    elif flag == 1:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k,flag=1)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n,flag=1)
                else:
                    print("your input atom number was out of range!!check!!")
    elif flag == 2:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k,flag=2)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n,flag=2)
                else:
                    print("your input atom number was out of range!!check!!")
    else:
        print("please write input flag parameters!!!!")

def PDOS(ini_data=ini_data):
    #get the shell or orbital projected PDOS
    PDOS = []
    fileopen = open("POSCAR", 'r+')
    pos_data = fileopen.readlines()
    fileopen.close()
    ele_name = pos_data[5].split()
    ele_num = list(map(int, pos_data[6].split()))
    def specific_pdos(ini_data=ini_data, flag=0):
        #get PDOS for all elements species
        if flag == 0: #shell projection for one specific element
            counter = -1
            tol_dos = [[0,0,0,0,0] for row in range(dos_num)] #shell projection for all atoms
            for i in range(len(ele_num)):
                ele_dos = [[0, 0, 0, 0, 0] for row in range(dos_num)]
                for j in range(ele_num[i]):
                    counter += 1
                    dos_begin = 5 + (counter+1)*dos_num + counter + 2
                    ini_data2 = ini_data[dos_begin:dos_begin+dos_num]
                    for k in range(dos_num):
                        tmp = list(map(float, ini_data2[k].split()))
                        ele_dos[k][0] = tmp[0]-fermi_level
                        tol_dos[k][0] = tmp[0]-fermi_level
                        shell_type = len(tmp)-1
                        if shell_type == 1:
                            ele_dos[k][1] += tmp[1]
                            tol_dos[k][1] += tmp[1]
                        elif shell_type == 4:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += sum(tmp[2:5])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += sum(tmp[2:5])
                        elif shell_type == 9:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += sum(tmp[2:5])
                            ele_dos[k][3] += sum(tmp[5:10])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += sum(tmp[2:5])
                            tol_dos[k][3] += sum(tmp[5:10])
                        elif shell_type == 16:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += sum(tmp[2:5])
                            ele_dos[k][3] += sum(tmp[5:10])
                            ele_dos[k][4] += sum(tmp[10:17])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += sum(tmp[2:5])
                            tol_dos[k][3] += sum(tmp[5:10])
                            tol_dos[k][4] += sum(tmp[10:17])
                        else:
                            print("shell nums out of range! impossible!!!")
                for i in range(dos_num):
                    for j in range(5):
                        ele_dos[i][j] = round(ele_dos[i][j], 4)
                        tol_dos[i][j] = round(tol_dos[i][j], 4)
                PDOS.append(ele_dos)
            

        elif flag == 1: #orbital projection for one specifc element
            counter = -1
            tol_dos = [[0 for column in range(17)] for row in range(dos_num)]
            for i in range(len(ele_num)):
                ele_dos = [[0 for column in range(17)] for row in range(dos_num)]
                for j in range(ele_num[i]):
                    counter += 1
                    dos_begin = 5 + (counter+1)*dos_num + counter + 2
                    ini_data2 = ini_data[dos_begin:dos_begin+dos_num]
                    for k in range(dos_num):
                        tmp = list(map(float, ini_data2[k].split()))
                        ele_dos[k][0] = tmp[0] - fermi_level
                        tol_dos[k][0] = tmp[0] - fermi_level
                        for m in range(1,len(tmp)):
                            ele_dos[k][m] += tmp[m]
                            tol_dos[k][m] += tmp[m]
                for i in range(dos_num):
                    for j in range(17):
                        ele_dos[i][j] = round(ele_dos[i][j], 4)
                        tol_dos[i][j] = round(tol_dos[i][j], 4)
                PDOS.append(ele_dos)
        else:
            print("your flag was out of range!!!")

        #output file
        if flag == 0:
            head = "#energy    s p d f\n"
        else:
            head = "#energy    s py pz px dxy dyz dz2 dxz dx2\n"
        for i in range(len(ele_name)):
            outname = "element" + ele_name[i] +"pdos.dat"
            ff = open(outname, "w+")
            ff.write(head)
            for j in range(dos_num):
                line = map(str, PDOS[i][j])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()
        ff = open("total_pdos.dat", "w+")
        ff.write(head)
        for i in range(dos_num):
            line = map(str, tol_dos[i])
            line = "    ".join(line)
            line += "\n"
            ff.write(line)
        ff.close()
    
    flag = int(input("choose type of element projection: 0)shell 1)orbital : "))
    if flag == 0:
        specific_pdos()
    elif flag == 1:
        specific_pdos(flag=1)

def spn_totalDOS():
    #get the spin total DOS of system
    totalDOS = []
    for i in range(dos_num):
        tmp = list(map(float, ini_data[6+i].split()[:4]))
        tmp[0] = round(tmp[0]-fermi_level, 4)
        tmp[1] = round(tmp[1], 4)
        tmp[2] = round(tmp[2], 4)
        tmp[3] = round(tmp[1] + tmp[2], 4)
        totalDOS.append(tmp)
    
    #write total DOS to total_dos.dat
    outname = "total_dos.dat"
    f = open(outname, "w+")
    f.write("#energy    up down total\n")
    for i in range(dos_num):
        line = map(str,totalDOS[i])
        line = "    ".join(line)
        line += "\n"
        f.write(line)
    f.close()

def spn_LDOS(ini_data=ini_data):
    #get the atom projected dos 
    LDOS = []

    def specific_ldos(atom_no, ini_data=ini_data, flag=0):
        #get spin LDOS for one specifc atom
        #get relative block in DOSCAR for specific atom
        atom_no = int(atom_no)
        dos_begin = 5 + (atom_no+1)*dos_num + atom_no + 2
        ini_data = ini_data[dos_begin:dos_begin+dos_num] #be aware of the section
        if flag == 0: #get total LDOS for one specific atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                counter,upcounter,downcounter = 0,0,0
                for j in range(1,len(tmp)):
                    if j%2 == 1:
                        upcounter += tmp[j]
                    else:
                        downcounter += tmp[j] 
                counter = round(upcounter+downcounter, 4)
                upcounter = round(upcounter, 4)
                downcounter = round(downcounter, 4)
                LDOS.append([tmp[0], upcounter, downcounter, counter])
        
            #write ldos of specific atom to output file
            outname = "atom" + str(atom_no) + "ldos.dat"
            f = open(outname, "w+")
            f.write("#energy    up down total\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                f.write(line)
            f.close()

        elif flag == 1: #get spin shell projected DOS for one specific atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                shell_type = len(tmp)-1
                if shell_type == 2:
                    ups = round(tmp[1], 4)
                    dws = round(tmp[2], 4)
                    s = round(sum(tmp[1:3]), 4)
                    LDOS.append([tmp[0], ups, dws, s])
                elif shell_type == 8:
                    ups = round(tmp[1], 4)
                    dws = round(tmp[2], 4)
                    s = round(sum(tmp[1:3]), 4)
                    upp = round(sum(tmp[3:9:2]), 4)
                    dwp = round(sum(tmp[4:10:2]), 4)
                    p = round(sum(tmp[3:9]), 4)
                    LDOS.append([tmp[0],ups,dws,s,upp,dwp,p])
                elif shell_type == 18:
                    ups = round(tmp[1], 4)
                    dws = round(tmp[2], 4)
                    s = round(sum(tmp[1:3]), 4)
                    upp = round(sum(tmp[3:9:2]), 4)
                    dwp = round(sum(tmp[4:10:2]), 4)
                    p = round(sum(tmp[3:9]), 4)
                    upd = round(sum(tmp[9:19:2]), 4)
                    dwd = round(sum(tmp[10:20:2]), 4)
                    d = round(sum(tmp[9:19]), 4)
                    LDOS.append([tmp[0],ups,dws,s,upp,dwp,p,upd,dwd,d])
                elif shell_type == 32:
                    ups = round(tmp[1], 4)
                    dws = round(tmp[2], 4)
                    s = round(sum(tmp[1:3]), 4)
                    upp = round(sum(tmp[3:9:2]), 4)
                    dwp = round(sum(tmp[4:10:2]), 4)
                    p = round(sum(tmp[3:9]), 4)
                    upd = round(sum(tmp[9:19:2]), 4)
                    dwd = round(sum(tmp[10:20:2]), 4)
                    d = round(sum(tmp[9:19]), 4)
                    upf = round(sum(tmp[19:33:2]), 4)
                    dwf = round(sum(tmp[20:34:2]), 4)
                    f = round(sum(tmp[19:33]), 4)
                    LDOS.append([tmp[0],ups,dws,s,upp,dwp,p,upd,dwd,d,upf,dwf,f])
                else:
                    print("shell nums out of range! impossible!!!")

            outname = "atom" + str(atom_no) + "ldos.dat"
            ff = open(outname, "w+")
            ff.write("#energy    ups downs s upp downp p upd downd d upf downf f\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()    
        elif flag == 2: #get orbital projected LDOS for one specifc atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                for j in tmp[1:]:
                    j = round(j, 4)
                LDOS.append(tmp)

            outname = "atom" + str(atom_no) + "ldos.dat"
            ff = open(outname, "w+")
            ff.write("#energy    ups downs uppy downpy uppz downpz uppx downpx updxy downdxy updyz downdyz updz2 downdz2 updxz downdxz updx2 downdx2\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()
        else:
            print("write set your flag!!!!")

    flag = int(input("choose type of atom projection: 0)total 1)shell 2)orbital:  "))
    if flag == 0:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n)
                else:
                    print("your input atom number was out of range!!check!!")
    elif flag == 1:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k,flag=1)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n,flag=1)
                else:
                    print("your input atom number was out of range!!check!!")
    elif flag == 2:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k,flag=2)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n,flag=2)
                else:
                    print("your input atom number was out of range!!check!!")
    else:
        print("please write input flag parameters!!!!")

def spn_PDOS(ini_data=ini_data):
    #get the spin shell or orbital projected PDOS
    PDOS = []
    fileopen = open("POSCAR", 'r+')
    pos_data = fileopen.readlines()
    fileopen.close()
    ele_name = pos_data[5].split()
    ele_num = list(map(int, pos_data[6].split()))
    def specific_pdos(ini_data=ini_data, flag=0):
        #get spin PDOS for all elements species
        if flag == 0: #shell projection for one specific element
            counter = -1
            tol_dos = [[0 for column in range(13)] for row in range(dos_num)] #shell projection for all atoms
            for i in range(len(ele_num)):
                ele_dos = [[0 for column in range(13)] for row in range(dos_num)]
                for j in range(ele_num[i]):
                    counter += 1
                    dos_begin = 5 + (counter+1)*dos_num + counter + 2
                    ini_data2 = ini_data[dos_begin:dos_begin+dos_num]
                    for k in range(dos_num):
                        tmp = list(map(float, ini_data2[k].split()))
                        ele_dos[k][0] = tmp[0]-fermi_level
                        tol_dos[k][0] = tmp[0]-fermi_level
                        shell_type = len(tmp)-1
                        if shell_type == 2:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += tmp[2]
                            ele_dos[k][3] += sum(tmp[1:3])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += tmp[2]
                            tol_dos[k][3] += sum(tmp[1:3])
                        elif shell_type == 8:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += tmp[2]
                            ele_dos[k][3] += sum(tmp[1:3])
                            ele_dos[k][4] += sum(tmp[3:9:2])
                            ele_dos[k][5] += sum(tmp[4:10:2])
                            ele_dos[k][6] += sum(tmp[3:9])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += tmp[2]
                            tol_dos[k][3] += sum(tmp[1:3])
                            tol_dos[k][4] += sum(tmp[3:9:2])
                            tol_dos[k][5] += sum(tmp[4:10:2])
                            tol_dos[k][6] += sum(tmp[3:9])
                        elif shell_type == 18:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += tmp[2]
                            ele_dos[k][3] += sum(tmp[1:3])
                            ele_dos[k][4] += sum(tmp[3:9:2])
                            ele_dos[k][5] += sum(tmp[4:10:2])
                            ele_dos[k][6] += sum(tmp[3:9])
                            ele_dos[k][7] += sum(tmp[9:19:2])
                            ele_dos[k][8] += sum(tmp[10:20:2])
                            ele_dos[k][9] += sum(tmp[9:19])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += tmp[2]
                            tol_dos[k][3] += sum(tmp[1:3])
                            tol_dos[k][4] += sum(tmp[3:9:2])
                            tol_dos[k][5] += sum(tmp[4:10:2])
                            tol_dos[k][6] += sum(tmp[3:9])
                            tol_dos[k][7] += sum(tmp[9:19:2])
                            tol_dos[k][8] += sum(tmp[10:20:2])
                            tol_dos[k][9] += sum(tmp[9:19])
                        elif shell_type == 32:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += tmp[2]
                            ele_dos[k][3] += sum(tmp[1:3])
                            ele_dos[k][4] += sum(tmp[3:9:2])
                            ele_dos[k][5] += sum(tmp[4:10:2])
                            ele_dos[k][6] += sum(tmp[3:9])
                            ele_dos[k][7] += sum(tmp[9:19:2])
                            ele_dos[k][8] += sum(tmp[10:20:2])
                            ele_dos[k][9] += sum(tmp[9:19])
                            ele_dos[k][10] += sum(tmp[19:33:2])
                            ele_dos[k][11] += sum(tmp[20:34:2])
                            ele_dos[k][12] += sum(tmp[19:33])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += tmp[2]
                            tol_dos[k][3] += sum(tmp[1:3])
                            tol_dos[k][4] += sum(tmp[3:9:2])
                            tol_dos[k][5] += sum(tmp[4:10:2])
                            tol_dos[k][6] += sum(tmp[3:9])
                            tol_dos[k][7] += sum(tmp[9:19:2])
                            tol_dos[k][8] += sum(tmp[10:20:2])
                            tol_dos[k][9] += sum(tmp[9:19])
                            tol_dos[k][10] += sum(tmp[19:33:2])
                            tol_dos[k][11] += sum(tmp[20:34:2])
                            tol_dos[k][12] += sum(tmp[19:33])
                        else:
                            print("shell nums out of range! impossible!!!")
                for i in range(dos_num):
                    for j in range(13):
                        ele_dos[i][j] = round(ele_dos[i][j], 4)
                        tol_dos[i][j] = round(tol_dos[i][j], 4)
                PDOS.append(ele_dos)
            

        elif flag == 1: #orbital projection for one specifc element
            counter = -1
            tol_dos = [[0 for column in range(33)] for row in range(dos_num)]
            for i in range(len(ele_num)):
                ele_dos = [[0 for column in range(33)] for row in range(dos_num)]
                for j in range(ele_num[i]):
                    counter += 1
                    dos_begin = 5 + (counter+1)*dos_num + counter + 2
                    ini_data2 = ini_data[dos_begin:dos_begin+dos_num]
                    for k in range(dos_num):
                        tmp = list(map(float, ini_data2[k].split()))
                        ele_dos[k][0] = tmp[0] - fermi_level
                        tol_dos[k][0] = tmp[0] - fermi_level
                        for m in range(1,len(tmp)):
                            ele_dos[k][m] += tmp[m]
                            tol_dos[k][m] += tmp[m]
                for i in range(dos_num):
                    for j in range(33):
                        ele_dos[i][j] = round(ele_dos[i][j], 4)
                        tol_dos[i][j] = round(tol_dos[i][j], 4)
                PDOS.append(ele_dos)
        else:
            print("your flag was out of range!!!")

        #output file 
        if flag == 0:
            head = "#energy    ups downs s upp downp p upd downd d upf downf f\n"
        else:
            head = "#energy    ups downs uppy downpy uppz downpz uppx downpx updxy downdxy updyz downdyz updz2 downdz2 updxz downdxz updx2 downdx2\n"
        for i in range(len(ele_name)):
            outname = "element" + ele_name[i] +"pdos.dat"
            ff = open(outname, "w+")
            ff.write(head)
            for j in range(dos_num):
                line = map(str, PDOS[i][j])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()
        ff = open("total_pdos.dat", "w+")
        ff.write(head)
        for i in range(dos_num):
            line = map(str, tol_dos[i])
            line = "    ".join(line)
            line += "\n"
            ff.write(line)
        ff.close()
    
    flag = int(input("choose type of element projection: 0)shell 1)orbital : "))
    if flag == 0:
        specific_pdos()
    elif flag == 1:
        specific_pdos(flag=1)

def ncs_LDOS(ini_data=ini_data):
    #get the no collinear atom projected dos 
    LDOS = []

    def specific_ldos(atom_no, ini_data=ini_data, flag=0):
        #get LDOS for one specifc atom
        #get relative block in DOSCAR for specific atom
        atom_no = int(atom_no)
        dos_begin = 5 + (atom_no+1)*dos_num + atom_no + 2
        ini_data = ini_data[dos_begin:dos_begin+dos_num] #be aware of the section
        if flag == 0: #get total LDOS for one specific atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                counter = 0
                for j in tmp[1:]:
                    counter += j
                counter = round(counter/2, 4)
                LDOS.append([tmp[0], counter])
        
            #write ldos of specific atom to output file
            outname = "atom" + str(atom_no) + "ldos.dat"
            f = open(outname, "w+")
            f.write("#energy    totaldos\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                f.write(line)
            f.close()

        elif flag == 1: #get shell projected DOS for one specific atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                shell_type = len(tmp)-1
                if shell_type == 4:
                    s = round(tmp[1], 4)
                    LDOS.append([tmp[0], s])
                elif shell_type == 16:
                    s = round(tmp[1], 4)
                    p = round(sum(tmp[5:14:4]), 4)
                    LDOS.append([tmp[0], s, p])
                elif shell_type == 36:
                    s = round(tmp[1], 4)
                    p = round(sum(tmp[5:14:4]), 4)
                    d = round(sum(tmp[17:34:4]), 4)
                    LDOS.append([tmp[0], s, p, d])
                elif shell_type == 64:
                    s = round(tmp[1], 4)
                    p = round(sum(tmp[5:14:4]), 4)
                    d = round(sum(tmp[17:34:4]), 4)
                    f = round(sum(tmp[37:62:4]), 4)
                    LDOS.append([tmp[0], s, p, d, f])
                else:
                    print("shell nums out of range! impossible!!!")

            outname = "atom" + str(atom_no) + "ldos.dat"
            ff = open(outname, "w+")
            ff.write("#energy    s p d f\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()    
        elif flag == 2: #get cartersian space projected LDOS for one specifc atom
            for i in range(dos_num):
                tmp = list(map(float, ini_data[i].split()))
                tmp[0] = round(tmp[0]-fermi_level, 4)
                for j in range(1,len(tmp)):
                    if (j-1)%4 == 0:
                        tmp[j] = ' '
                tmp = [ x for x in tmp if x != ' ']        
                for j in tmp[1:]:        
                    j = round(j, 4)
                LDOS.append(tmp)

            outname = "atom" + str(atom_no) + "ldos.dat"
            ff = open(outname, "w+")
            ff.write("#energy    x y z cartesian component of s py pz px dxy dyz dz2 dxz dx2\n")
            for i in range(dos_num):
                line = map(str, LDOS[i])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()
        else:
            print("write set your flag!!!!")

    flag = int(input("choose type of atom projection: 0)total 1)shell 2)cartesian:  "))
    if flag == 0:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n)
                else:
                    print("your input atom number was out of range!!check!!")
    elif flag == 1:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k,flag=1)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n,flag=1)
                else:
                    print("your input atom number was out of range!!check!!")
    elif flag == 2:
        atom_nos = input("choose which atom to project DOS: )atom number or all    ")
        if atom_nos == "all":
            for k in range(atom_num):
                specific_ldos(atom_no=k,flag=2)
        else:
            for n in atom_nos.strip().split():
                if int(n)<=atom_num-1 and int(n) >= 0:
                    specific_ldos(atom_no=n,flag=2)
                else:
                    print("your input atom number was out of range!!check!!")
    else:
        print("please write input flag parameters!!!!")

def ncs_PDOS(ini_data=ini_data):
    #get the shell or cartesian projected PDOS
    PDOS = []
    fileopen = open("POSCAR", 'r+')
    pos_data = fileopen.readlines()
    fileopen.close()
    ele_name = pos_data[5].split()
    ele_num = list(map(int, pos_data[6].split()))
    def specific_pdos(ini_data=ini_data, flag=0):
        #get no colllniear PDOS for all elements species
        if flag == 0: #shell no collilnear projection for one specific element
            counter = -1
            tol_dos = [[0,0,0,0,0] for row in range(dos_num)] #shell projection for all atoms
            for i in range(len(ele_num)):
                ele_dos = [[0, 0, 0, 0, 0] for row in range(dos_num)]
                for j in range(ele_num[i]):
                    counter += 1
                    dos_begin = 5 + (counter+1)*dos_num + counter + 2
                    ini_data2 = ini_data[dos_begin:dos_begin+dos_num]
                    for k in range(dos_num):
                        tmp = list(map(float, ini_data2[k].split()))
                        ele_dos[k][0] = tmp[0]-fermi_level
                        tol_dos[k][0] = tmp[0]-fermi_level
                        shell_type = len(tmp)-1
                        if shell_type == 4:
                            ele_dos[k][1] += tmp[1]
                            tol_dos[k][1] += tmp[1]
                        elif shell_type == 16:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += sum(tmp[5:14:4])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += sum(tmp[5:14:4])
                        elif shell_type == 36:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += sum(tmp[5:14:4])
                            ele_dos[k][3] += sum(tmp[17:34:4])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += sum(tmp[5:14:4])
                            tol_dos[k][3] += sum(tmp[17:34:4])
                        elif shell_type == 64:
                            ele_dos[k][1] += tmp[1]
                            ele_dos[k][2] += sum(tmp[5:14:4])
                            ele_dos[k][3] += sum(tmp[17:34:4])
                            ele_dos[k][4] += sum(tmp[37:62:4])
                            tol_dos[k][1] += tmp[1]
                            tol_dos[k][2] += sum(tmp[5:14:4])
                            tol_dos[k][3] += sum(tmp[17:34:4])
                            tol_dos[k][4] += sum(tmp[37:62:4])
                        else:
                            print("shell nums out of range! impossible!!!")
                for i in range(dos_num):
                    for j in range(5):
                        ele_dos[i][j] = round(ele_dos[i][j], 4)
                        tol_dos[i][j] = round(tol_dos[i][j], 4)
                PDOS.append(ele_dos)
            

        elif flag == 1: #orbital projection for one specifc element
            counter = -1
            tol_dos = [[0 for column in range(65)] for row in range(dos_num)]
            for i in range(len(ele_num)):
                ele_dos = [[0 for column in range(65)] for row in range(dos_num)]
                for j in range(ele_num[i]):
                    counter += 1
                    dos_begin = 5 + (counter+1)*dos_num + counter + 2
                    ini_data2 = ini_data[dos_begin:dos_begin+dos_num]
                    for k in range(dos_num):
                        tmp = list(map(float, ini_data2[k].split()))
                        ele_dos[k][0] = tmp[0] - fermi_level
                        tol_dos[k][0] = tmp[0] - fermi_level
                        for m in range(1,len(tmp)):
                            ele_dos[k][m] += tmp[m]
                            tol_dos[k][m] += tmp[m]
                for i in range(dos_num):
                    for j in range(65):
                        if (j-1)%4 == 0:
                            ele_dos[i][j] = ' '
                        else:
                            ele_dos[i][j] = round(ele_dos[i][j], 4)
                            tol_dos[i][j] = round(tol_dos[i][j], 4)
                ele_dos = [[j for j in ele_dos[i] if j != ' '] for i in range(dos_num)]
                PDOS.append(ele_dos)
            for i in range(dos_num):
                for j in range(1,65):
                    if (j-1)%4 == 0:
                        tol_dos[i][j] = ' '
            tol_dos = [[j for j in tol_dos[i] if j != ' '] for i in range(dos_num)]
        else:
            print("your flag was out of range!!!")

        #output file 
        if flag == 0:
            head = "#energy    s p d f\n"
        else:
            head = "#energy x y z cartesian component of s py pz px dxy dyz dz2 dxz dx2\n"
        for i in range(len(ele_name)):
            outname = "element" + ele_name[i] +"pdos.dat"
            ff = open(outname, "w+")
            ff.write(head)
            for j in range(dos_num):
                line = map(str, PDOS[i][j])
                line = "    ".join(line)
                line += "\n"
                ff.write(line)
            ff.close()
        ff = open("total_pdos.dat", "w+")
        ff.write(head)
        for i in range(dos_num):
            line = map(str, tol_dos[i])
            line = "    ".join(line)
            line += "\n"
            ff.write(line)
        ff.close()
    
    flag = int(input("choose type of element projection: 0)shell 1)cartesian : "))
    if flag == 0:
        specific_pdos()
    elif flag == 1:
        specific_pdos(flag=1)

#main program start
switch1 = int(input("choose type of your calculation: 0)no magnetic 1)linear spin 2)nolinear or SOC : "))
switch2 = int(input("choose type of projection : 0)total DOS 1)LDOS 2)PDOS : "))
if switch1 == 0:
    if switch2 == 0:
        totalDOS()
    elif switch2 == 1:
        LDOS()
    elif switch2 == 2:
        PDOS()
    else:
        print("please properly input your case, boy!!!")
elif switch1 == 1:
    if switch2 == 0:
        spn_totalDOS()
    elif switch2 == 1:
        spn_LDOS()
    elif switch2 == 2:
        spn_PDOS()
    else:
        print("please properly input your case, boy!!!")
elif switch1 == 2:
    if switch2 == 0:
        totalDOS()
    elif switch2 == 1:
        ncs_LDOS()
    elif switch2 == 2:
        ncs_PDOS()
    else:
        print("please properly input your case, boy!!!")
else:
    print("please properly input your case, boy!!!")

