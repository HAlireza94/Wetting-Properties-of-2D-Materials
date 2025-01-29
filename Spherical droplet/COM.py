

x, y, z = [], [], []
with open('checking.gro') as f:
    for line in f:
        p = line.split()
        if len(p) == 6:
            x.append(float(p[3]))
            y.append(float(p[4]))
            z.append(float(p[5]))

Mw = 15.999

X, Y, Z = 0, 0, 0
for i in range(len(x)):
    X += x[i]*Mw
    Y += y[i]*Mw
    Z += z[i]*Mw

#print('X_com = '+str(X/(len(x)*Mw)))
#print('Y_com = '+str(Y/(len(y)*Mw)))
#print('Z_com = '+str(Z/(len(z)*Mw)))

print(f"{X/(len(x)*Mw):<25}{Y/(len(y)*Mw):<25}{Z/(len(z)*Mw):<25}")



