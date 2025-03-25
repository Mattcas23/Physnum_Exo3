import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/a-x-3/Desktop/Ex3_2024_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/a-x-3/Desktop/Ex3_2024_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file

#----------------------------------------- Valeurs du Configfile --------------------------------- # 

Values = np.genfromtxt("configuration.in.example" , comments = '//')

tfin = Values[0,-1]
ms = Values[1,-1]
mj = Values[2,-1]
rs= Values[3,-1]
rj= Values[4,-1]
G_grav = Values[5,-1]
a   =  Values[6,-1]
eps    = Values[7,-1]
alpha  = Values[8,-1]
maxit  = Values[9,-1]
output = Values[10,-1]
sampling = Values[11,-1]
vx0 = Values[12,-1]
vy0 = Values[13,-1]
x0 = Values[14,-1]
y0 = Values[15,-1]
#output = Values[16,-1]
#nsteps = Values[17,-1]
f = Values[18,-1]
adaptative = Values[19,-1]

# ---------------------------------------------------------------

nsteps = np.array([30e03]) # TODO change
neps = np.array([0])
# nsteps = np.array([ 20e3, 40e3, 50e3 , 60e3 , 80e3 , 200e3])
# nsteps = np.array([4000,40000])
#nsteps = np.array([200e3])
# nsteps = np.array([4010, 4050])
nsimul = len(nsteps)  # Number of simulations to perform

tfin = 7776000  # Done : Verify that the value of tfin is EXACTLY the same as in the input file

dt = tfin / nsteps

energy = np.zeros(nsimul) # added 

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

paramstr2 = 'neps'
param2 = neps
nsimul2 = len(neps)

# ------------------------------------------- Simulations --------------------------------------------- #

if adaptative : # Simulations avec pas de temps adaptatif

    print("Adaptatif")
    
    outputs = []  # List to store output file names
    convergence_list = []
    for i in range(nsimul):
        output_file = f"{paramstr2}={param2[i]}.out"
        outputs.append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} {paramstr2}={param2[i]:.15g} output={output_file}"
        cmd = f"{executable} {input_filename} {paramstr}={param2[i]:.15g} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')

else : # Simulations avec pas de temps fixe

    print("Fixe")

    outputs = []  # List to store output file names
    convergence_list = []
    for i in range(nsimul):
        output_file = f"{paramstr}={param[i]}.out"
        outputs.append(output_file)
        cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
        cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
        print(cmd)
        subprocess.run(cmd, shell=True)
        print('Done.')    

error = np.zeros(nsimul)

#dref  = np.loadtxt(outputs[-1])
#xref = dref[-1,3]
#yref = dref[-1,4]
#Eref = dref[-1,5]

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]
    vx = data[-1, 1]  # final position, velocity, energy
    vy = data[-1, 2]
    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list.append(xx)
    # TODO compute the error for each simulation
    #error[i] = pow( (yref-yy)**2 + (xref-xx)**2 ,0.5)  
    #energy[i] = abs(Eref-En)

lw = 1.5
fs = 16

# ---------------------------------------- Zones Plots ---------------------------------------- #

def Trajectoire () :

    plt.figure()
    #Soleil = plt.Circle((38e7 * 5.972e24 / ( 5.972e24 + 7.348e22 ),0),1737100, color = 'yellow' , label = "Soleil")
    #Jupyter = plt.Circle((-38e7 * 7.348e22 / ( 5.972e24 + 7.348e22 ),0),6378100 , color = 'brown' , label = "Jupyter")    
    plt.plot(data[:, 3], data[:, 4], color = 'black' , label = '$n_{step} = $' + f"{nsteps[0]:.0f}")
    plt.xlabel('x [m]', fontsize=fs)
    plt.ylabel('y [m]', fontsize=fs)
    plt.legend()

def Energie () : # Energie en fonction du temps 

    plt.figure()
    plt.plot(data[:, 0], data[:, 5], color = 'black' , label = '$n_{step} = $' + f"{nsteps[0]:.0f}")
    plt.xlabel('t [s]', fontsize=fs)
    plt.ylabel('$E_{mec}$', fontsize=fs)
    plt.legend()

def PosFin_Conv ( norder = 1 ) : # convergeance sur la postion finale ( en x )

    plt.figure()
    plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
    plt.ticklabel_format(axis='y', style='scientific', scilimits = (-4,-4))
    plt.xlabel(f"$(\\Delta t)^{norder}$ [s]", fontsize=fs)
    plt.ylabel('$x_{final}$', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)

def dts () : # Pas de temps au cours du temps

    plt.figure()
    plt.plot( t, data[:,-1] , color = 'black' , label = '$\\epsilon = $' + f'{eps}') # à modifier 
    plt.xlabel('t', fontsize=fs)
    plt.ylabel('Pas de temps $dt$', fontsize=fs)    
    plt.legend()

def x() :

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(t, data[:, 3], color = 'black' , label = '$n_{step} = $' + f"{nsteps[0]:.0f}")
    ax.set_xlabel('t [s]', fontsize=fs)
    ax.set_ylabel('x [m]', fontsize=fs)
    plt.legend()


Trajectoire ()
Energie ()
PosFin_Conv ()
dts () 
x() 

plt.show()




















def energies (nord = 1) : 

    plt.figure()
    plt.loglog(dt[:-1],energy[:-1], 'r+-' , linewidth = lw)
    plt.xlabel(f"$\\Delta t $ [s]", fontsize = fs)
    plt.ylabel('error $E_{mec}$', fontsize=fs)
    plt.loglog(dt*nord, pow(dt,nord), color = 'black' ,linewidth = lw , label = f"$1/N^{nord}$" , linestyle = 'dashed')
    #plt.ticklabel_format(axis='y', style='scientific', scilimits = (5,5))
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.legend()

#energies() 

# uncomment the following if you want debug
#import pdb
#pbd.set_trace()

def pos_error (n) : 
    plt.figure()
    plt.loglog(dt[:-1], error[:-1], 'r+-', linewidth=lw)
    plt.loglog(dt, pow(dt,n), color = 'black' ,linewidth = lw , label = f"$1/N^{n}$" , linestyle = 'dashed')
    plt.xlabel('$\\Delta t$ [s]', fontsize=fs)
    plt.ylabel('final position error [m]', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.grid(True)
    plt.legend()

#pos_error (1)

"""
Si on n'a pas la solution analytique: on représente la quantite voulue
(ci-dessous v_y, modifier selon vos besoins)
en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 1 # Modify if needed = 1 (original)

plt.figure()
plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
plt.xlabel(f"$(\\Delta t)^{norder}$ [s]", fontsize=fs)
plt.ylabel('$x$ [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

#plt.show()
