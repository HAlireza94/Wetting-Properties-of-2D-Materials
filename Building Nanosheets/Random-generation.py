import numpy as np
import pandas as pd

proportion = 0.25 # 25 % just O and 75 % OH 
atoms, index, x, y, z = [], [], [], [], []
with open('new-system.pdb') as f:
    for line in f:
        p = line.split()
        if len(p) == 8:
            index.append(p[1])
            atoms.append(p[2])
            x.append(p[5])
            y.append(p[6])
            z.append(p[7])
        elif len(p) == 10:
            Lx = float(p[3])
            Ly = float(p[6])
            Lz = float(p[9])
            
            
            

############### Hydrogens on Top ###############
index_h_top, atoms_h_top, x_h_top, y_h_top, z_h_top = [], [], [], [], []
for i in range(len(atoms)):
    if float(z[i]) > 17:
        index_h_top.append(index[i])
        atoms_h_top.append(atoms[i])
        x_h_top.append(x[i])
        y_h_top.append(y[i])
        z_h_top.append(z[i])
        
integer_index_h_top = []
for i in range(len(index_h_top)):
    integer_index_h_top.append(int(index_h_top[i]))

removing_h_top = np.random.choice(
    integer_index_h_top, size=int(len(integer_index_h_top) * proportion), replace=False
).tolist()
rht = [str(d) for d in removing_h_top]

############### Hydrogens on Bottom ###############
index_h_down, atoms_h_down, x_h_down, y_h_down, z_h_down = [], [], [], [], []
for i in range(len(atoms)):
    if float(z[i]) < 6:
        index_h_down.append(index[i])
        atoms_h_down.append(atoms[i])
        x_h_down.append(x[i])
        y_h_down.append(y[i])
        z_h_down.append(z[i])
        
integer_index_h_down = []
for i in range(len(index_h_down)):
    integer_index_h_down.append(int(index_h_down[i]))
removing_h_down = np.random.choice(
    integer_index_h_down, size=int(len(integer_index_h_down) * proportion), replace=False
).tolist()
rhd = [str(d) for d in removing_h_down]

############## Generating New Set ###############
df = pd.DataFrame({'index':index,'atoms':atoms,'X':x, 'Y':y,'Z':z})
f1 = df[~df['index'].isin(rht)]
f2 = f1[~f1['index'].isin(rhd)]


filename = "modified-system.pdb"
with open(filename, 'w') as file:
    file.write("REMARK    bx = {}  by = {}  bz = {}\n".format(Lx, Ly, Lz))
    file.write("REMARK    # of atoms = {}\n".format(len(f2)))

    # Loop through rows using zip for cleaner code
    for index, atom, x, y, z in zip(f2['index'], f2['atoms'], f2['X'], f2['Y'], f2['Z']):
        file.write(
            "ATOM  {:>5} {:<4} XYZ {:>4}    {:8.3f}{:8.3f}{:8.3f}\n".format(
                int(index), atom, 1, float(x), float(y), float(z)
            )
        )



