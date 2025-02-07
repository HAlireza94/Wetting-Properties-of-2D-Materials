import numpy as np
index, atom, x, y, z = [], [], [], [], []
with open('oh-sup2.gro') as f:
    for line in f:
        p = line.split()
        if len(p) == 9:
            index.append(p[2])
            atom.append(p[1])
            x.append(np.round(float(p[3])*10,3))
            y.append(np.round(float(p[4])*10,3))
            z.append(np.round(float(p[5])*10,3))

x_up, y_up, z_up, index_up, atoms_up =[], [], [], [], []
for i in range(len(x)):
    if z[i] > 11:
        x_up.append(x[i])
        y_up.append(y[i])
        z_up.append(z[i])
        index_up.append(index[i])
        atoms_up.append(atom[i])

x_down, y_down, z_down, index_down, atoms_down =[], [], [], [], []
for i in range(len(x)):
    if z[i] < 7.5:
        x_down.append(x[i])
        y_down.append(y[i])
        z_down.append(z[i])
        index_down.append(index[i])
        atoms_down.append(atom[i])




######################## Ti, O, and H on top ########################

## finding Ti ###
x_up_ti, y_up_ti, z_up_ti, index_up_ti, atoms_up_ti = [], [], [], [], []
for i in range(len(x_up)):
    if z_up[i] > 11 and z_up[i] < 13:
        x_up_ti.append(x_up[i])
        y_up_ti.append(y_up[i])
        z_up_ti.append(z_up[i])
        index_up_ti.append(index_up[i])
        atoms_up_ti.append(atoms_up[i])

## finding O ###
x_up_O, y_up_O, z_up_O, index_up_O, atoms_up_O = [], [], [], [], []
for i in range(len(x_up)):
    if z_up[i] > 13.9 and z_up[i] < 14.9:
        x_up_O.append(x_up[i])
        y_up_O.append(y_up[i])
        z_up_O.append(z_up[i])
        index_up_O.append(index_up[i])
        atoms_up_O.append(atoms_up[i])

x_up_h, y_up_h, z_up_h, index_up_h, atoms_up_h = [], [], [], [], []
for i in range(len(x_up)):
    if z_up[i] > 14.9:
        x_up_h.append(x_up[i])
        y_up_h.append(y_up[i])
        z_up_h.append(z_up[i])
        index_up_h.append(index_up[i])
        atoms_up_h.append(atoms_up[i])


######################## Ti, O, and H on bottom ########################

## finding Ti ###
x_down_ti, y_down_ti, z_down_ti, index_down_ti, atoms_down_ti = [], [], [], [], []
for i in range(len(x_down)):
    if z_down[i] > 6 and z_down[i] < 7:
        x_down_ti.append(x_down[i])
        y_down_ti.append(y_down[i])
        z_down_ti.append(z_down[i])
        index_down_ti.append(index_down[i])
        atoms_down_ti.append(atoms_down[i])

## finding O ###
x_down_O, y_down_O, z_down_O, index_down_O, atoms_down_O = [], [], [], [], []
for i in range(len(x_down)):
    if z_down[i] > 4 and z_down[i] < 5:
        x_down_O.append(x_down[i])
        y_down_O.append(y_down[i])
        z_down_O.append(z_down[i])
        index_down_O.append(index_down[i])
        atoms_down_O.append(atoms_down[i])

x_down_h, y_down_h, z_down_h, index_down_h, atoms_down_h = [], [], [], [], []
for i in range(len(x_down)):
    if z_down[i] < 4:
        x_down_h.append(x_down[i])
        y_down_h.append(y_down[i])
        z_down_h.append(z_down[i])
        index_down_h.append(index_down[i])
        atoms_down_h.append(atoms_down[i])


print("[ atoms ]")
print("; nr type  resnr    residue    atom     cgnr    charge       mass")

tota_charge = 0
for i in range(len(atom)):
    if atom[i] == 'C':
        q = -0.76
        Mw = 12.0107
        at = 'C'
    elif atom[i] == 'O':
        q = -0.79   
        Mw = 15.9994
        at = 'O'
    elif atom[i] == 'H':
        q = 0.35  
        Mw = 1.008
        at = 'H'
        
    else:
        if index[i] in index_up_ti:
            q = 0.88     
            Mw = 47.867
            at = 'Ti2'
        elif index[i] in index_down_ti:
            q = 0.88     
            Mw = 47.867
            at = 'Ti2'
        else:
            q  = 0.64     
            Mw = 47.867
            at = 'Ti1'
                
    tota_charge += q  
    res = 'MAA'
    print(f"{i+1:<7}{at:<7}{1:<7}{res:<7}{at:<7}{i+1:<7}{q:<7}{Mw:<4}")





print("")
print("")

print("[ bonds ]")
print("")


####### O-H & O-Ti bonds on Up ######
c_bond = 0
list1_top = []
for i in range(len(x_up)):
    for j in range(i,len(x_up)):
        if (np.round(np.sqrt((x_up[i]-x_up[j])**2 + (y_up[i]-y_up[j])**2 + (z_up[i]-z_up[j])**2),3)) < 1.1 and i != j:
            c_bond += 1
            sign = ";H-O"
            list1_top.append((index_up[i],index_up[j]))
            print(f"{index_up[i]:<5}{index_up[j]:<5}{1:<5}{0.0973:<10}{23185.0084:<15}{sign:<5}")



list2_top = []
for i in range(len(x_up_O)):
    for j in range(len(x_up_O)):
        if np.round(np.sqrt((x_up_ti[j] - x_up_O[i]) ** 2 + (y_up_ti[j] - y_up_O[i]) ** 2 + (z_up_ti[j] - z_up_O[i]) ** 2),3) == 2.9:
            sign = ";O - Ti"
            c_bond += 1
            list2_top.append((index_up_O[i],index_up_ti[j]))
            print(f"{index_up_O[i]:<5}{index_up_ti[j]:<5}{1:<5}{0.2189:<10}{5208.6616:<10}{sign:<5}")
            


####### O-H and O-Ti bonds on bottom ######
list1_bot = []
for i in range(len(x_down)):
    for j in range(i,len(x_down)):
        if (np.round(np.sqrt((x_down[i]-x_down[j])**2 + (y_down[i]-y_down[j])**2 + (z_down[i]-z_down[j])**2),3)) < 1 and i != j:
            # print((np.round(np.sqrt((x_down[i]-x_down[j])**2 + (y_down[i]-y_down[j])**2 + (z_down[i]-z_down[j])**2),3)))
            c_bond += 1
            sign = ";H-O"
            list1_bot.append((index_down[i],index_down[j]))
            print(f"{index_down[i]:<5}{index_down[j]:<5}{1:<5}{0.0973:<10}{23185.0084:<12}{sign:<5}")




list2_bot = []
for i in range(len(x_down_O)):
    for j in range(len(x_down_O)):
        distance = (np.round(np.sqrt((x_down_ti[j] - x_down_O[i]) ** 2 + (y_down_ti[j] - y_down_O[i]) ** 2 + (z_down_ti[j] - z_down_O[i]) ** 2),3))
        if distance == 2.907:
            c_bond += 1
            sign = ";O - Ti"
            list2_bot.append((index_down_O[i],index_down_ti[j]))
            print(f"{index_down_O[i]:<5}{index_down_ti[j]:<5}{1:<5}{0.2189:<10}{5208.6616:<10}{sign:<5}")
            

print("")
print("")
print("[ angles ]")
print("")              
list2_dict_top = {row[0]: row[1] for row in list2_top}
c_angle = 0
for col1, col2 in list1_top:
    if col2 in list2_dict_top:
        a1 = atom[int(col1)-1]
        a2 = atom[int(col2)-1]
        a3 = atom[int(list2_dict_top[col2])-1]
        c_angle += 1
        print(f"{col1:<10}{col2:<10}{list2_dict_top[col2]:<10}{1:<5}{125.87:<10}{193996.6709:<15}{';'+a1 + ' - ' + a2 +  ' - ' + a3:<5}")

list2_dict_bot = {row[0]: row[1] for row in list2_bot}
for col1, col2 in list1_bot:
    if col2 in list2_dict_bot:
        a1 = atom[int(col1)-1]
        a2 = atom[int(col2)-1]
        a3 = atom[int(list2_dict_bot[col2])-1]
        c_angle += 1
        print(f"{col1:<10}{col2:<10}{list2_dict_bot[col2]:<10}{1:<5}{125.87:<10}{193996.6709:<15}{';'+a1 + ' - ' + a2 +  ' - ' + a3:<5}")
