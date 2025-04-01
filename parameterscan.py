import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/Matte/Desktop/EPFL/BA4/Physnum/Exo3/exe" #r"/Users/a-x-3/Desktop/Ex3_2024_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/Matte/Desktop/EPFL/BA4/Physnum/Exo3"#r"/Users/a-x-3/Desktop/Ex3_2024_student"
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
#eps    = Values[7,-1]
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
jupyter = Values[20,-1]
ma = Values[21,-1]
Rg = Values[22,-1]

# ---------------------------------------------------------------

nsteps = np.array([ 60 ])*1e3

eps = np.array( [1000 , 3000 ,  5000 , 10000 , 20000])

nsimul = len(nsteps)  # Number of simulations to perform ( pour un nombre de pas changeant )

dt = tfin / nsteps

energy = np.zeros(nsimul) # added 

paramstr = 'nsteps'  # Parameter name to scan
param = nsteps  # Parameter values to scan

paramstr2 = 'eps'
param2 = eps
nsimul2 = len(eps)

# ------------------------------------------- Simulations ---------------------------------------------

if jupyter :
 print("Avec Jupyter")
 if Rg :
     print("Référentiel RG")
 else : 
     print("Référentiel R'")
else :
 print("Sans Jupyter")



if adaptative : # Simulations avec pas de temps adaptatif

 print("Pas de temps adaptatif")

 outputs = []  # List to store output file names
 convergence_list = []
 for i in range(nsimul2):
     output_file = f"{paramstr2}={param2[i]}.out"
     outputs.append(output_file)
     cmd = f"{repertoire}{executable} {input_filename} {paramstr2}={param2[i]:.15g} output={output_file}"
     cmd = f"{executable} {input_filename} {paramstr2}={param2[i]:.15g} output={output_file}"
     print(cmd)
     subprocess.run(cmd, shell=True)
     print('Done.')

else : # Simulations avec pas de temps fixe

 print("Pas de temps fixe")

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

Eerr = np.array([]) # erreur sur l'énergie mécanique

jsteps_list = np.zeros(len(eps))

if adaptative : # schéma à pas de temps adaptatif 

 for i in range(nsimul2):  # Iterate through the results of all simulations
     data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
     t = data[:, 0]
     vx = data[-1, 1]  # final position, velocity, energy
     vy = data[-1, 2]
     xx = data[-1, 3]
     yy = data[-1, 4]
     En = data[-1, 5]
     jsteps = data[-1,6] # nombre de pas de temps total pour la simulation 

     convergence_list.append(xx)
     jsteps_list[i] = jsteps
     Eerr = np.append( Eerr , abs(data[1,5] - data[-1,5]) ) # erreur sur l'Emec : différence entre la valeur initiale et la valeur finale 

else : # schéma à pas de temps fixe 

 for i in range(nsimul):
     data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
     t = data[:, 0]
     vx = data[-1, 1]  # final position, velocity, energy
     vy = data[-1, 2]
     xx = data[-1, 3]
     yy = data[-1, 4]
     En = data[-1, 5]

     convergence_list.append(xx)
     Eerr = np.append( Eerr , abs(data[1,5] - data[-1,5]) ) # erreur sur l'Emec : différence entre la valeur initiale et la valeur finale

lw = 1.5
fs = 16

# ---------------------------------------- Zones Plots ---------------------------------------- #

# Vitesse max et min en y 
print( f" max(vy) = {max(abs(data[:,2]))}" )
print( f" min(vy) = {min(abs(data[:,2]))}" )
# Position max et min en x
print( f" max(x) = {max((data[:,3])):.3e}" )
print( f" min(x) = {min((data[:,3])):.3e}" )

def Trajectoire () :
 

    # Calcul des distances au Soleil pour chaque point de trajectoire
    r = np.sqrt(data[:,3]**2 + data[:,4]**2)

    # Distance maximale (aphélie)
    r_max = np.max(r)
    idx_max = np.argmax(r)  # Indice correspondant

    # Distance minimale (périhélie)
    r_min = np.min(r)
    idx_min = np.argmin(r)  # Indice correspondant

    # Coordonnées correspondantes
    x_max, y_max = data[idx_max, 3], data[idx_max, 4]
    x_min, y_min = data[idx_min, 3], data[idx_min, 4]

    # Affichage des résultats
    print(f"Aphélie: r_max = {r_max:.2e} m, coordonnées (x, y) = ({x_max:.2e}, {y_max:.2e})")
    print(f"Périhélie: r_min = {r_min:.2e} m, coordonnées (x, y) = ({x_min:.2e}, {y_min:.2e})")

    plt.figure()
    Soleil  = plt.scatter(  - a * mj / (mj + ms) , 0 , marker = 'o' , s=10,  color = 'gold' , label = 'Soleil', zorder=3 )

    if jupyter : # On affiche Jupyter

     referentiel = "$R'$"
     if Rg :
         referentiel = "$R_G$"
     plt.title(f"Référentiel {referentiel}",fontsize = fs - 2)
     Jupyter = plt.scatter( a * ms / ( ms + mj ) ,0,marker = 'o' , color = 'brown',  label = "Jupyter" )

    else : # On affiche la position minimale et maximale en x # -1.51e12 pk?

     rmax = plt.scatter( data[idx_max, 3], data[idx_max, 4], marker = '^' , s=10, color = 'red' , label = "$x_{max} =$" + f"{max((data[:,3])):.2e}" , zorder=3 )
     rmin = plt.scatter( data[idx_min, 3], data[idx_min, 4], marker = 'o' , s=10, color = 'red', label = "$x_{min} =$" + f"{min((data[:,3])):.2e}" , zorder=3)
     plt.plot([x_min, x_max], [y_min, y_max], 'c-', label="Axe Aphélie-Périhélie")
    #posinit = plt.scatter(data[0,3],data[0,4], marker = 'o' , color = 'grey' , label = "astéroide")
    plt.plot(data[:, 3], data[:, 4], color = 'black' ,  label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")
    plt.xlabel('x [m]', fontsize=fs)
    plt.ylabel('y [m]', fontsize=fs)
    plt.axis('equal')
    plt.grid()
    plt.legend()

def Energie () : # Energie en fonction du temps 

 plt.figure()
 plt.plot(data[:, 0], data[:, 5], color = 'black' , label = '$n_{step} = $' + f"{nsteps[-1]:.0f}")
 plt.xlabel('t [s]', fontsize=fs)
 plt.ylabel('$E_{mec}$', fontsize=fs)
 plt.legend()

def Emec_Err ( norder = 5 ) : # Erreur de l'Emec en fonction du temps ( fixe et sans jupyter : ordre 5 ; adaptatif : ordre ? )

 if adaptative : # Erreur avec pas de temps adaptatif

     plt.figure()
     plt.loglog( jsteps_list , Eerr,'k+-',linewidth = lw)
     plt.loglog( jsteps_list , (1/pow(jsteps_list,norder))*1e14 , color  = 'red' , label = f"$N^{norder}$" , linestyle = 'dashed')
     plt.xlabel("$j_{steps}$", fontsize=fs)
     plt.ylabel('$\\delta_{E_{mec}}$', fontsize=fs)
     plt.legend(fontsize = fs - 3)
     plt.plot()

 else : # Erreur avec pas de temps fixe  

     plt.figure()
     plt.loglog( dt,Eerr,'k+-',linewidth = lw)
     plt.loglog( dt , pow(dt,norder)/1e20 , color  = 'red' , label = f"$1/N^{norder}$" , linestyle = 'dashed')
     plt.xlabel(f"$\\Delta t$ [s]", fontsize=fs)
     plt.ylabel('$\\delta_{E_{mec}}$', fontsize=fs)
     plt.legend(fontsize = fs - 3)
     plt.plot()

def PosFin_Conv ( norder = 4 ) : # convergeance sur la postion finale ( en x )


 if adaptative : # convergeance en fonction de jsteps ( Runge Kutta adaptatif : ordre 1 ? )

     plt.figure()
     plt.plot(jsteps_list**norder, convergence_list , 'k+-' , linewidth = lw)
     plt.xlabel(f"$(jsteps)^{norder}$", fontsize=fs)
     plt.ylabel('$x_{final}$', fontsize=fs)
     plt.plot()

 else : # convergeance en fonction de dt ( Runge Kutta fixe : ordre 4)

     plt.figure()
     plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
     plt.xlabel(f"$(\\Delta t)^{norder}$ [s]", fontsize=fs)
     plt.ylabel('$x_{final}$', fontsize=fs)
     plt.xticks(fontsize=fs)
     plt.yticks(fontsize=fs)
     plt.grid(True)
     plt.plot()

def dts ( jstep ) : # Pas de temps au cours du temps et jsteps au cours du temps si adaptatif et jstep

 plt.figure()
 plt.plot( t, data[:,-1] , color = 'black' , label = '$\\epsilon = $' + f'{eps[-1]}') # à modifier 
 plt.xlabel('t [s] ', fontsize=fs)
 plt.ylabel('$dt$ [s]', fontsize=fs)    
 plt.legend()

 if adaptative and jstep :

     plt.figure()
     plt.plot(t,data[:,6], color = 'black' , label = '$\\epsilon = $' + f'{eps[-1]}' )
     plt.xlabel('t [s]', fontsize=fs)
     plt.ylabel('$j_{steps}$', fontsize=fs)    
     plt.legend()


def x() :

 fig, ax = plt.subplots(constrained_layout=True)
 ax.plot(t, data[:, 3], color = 'black' , label = '$n_{step} = $' + f"{nsteps[0]:.0f}")
 ax.set_xlabel('t [s]', fontsize=fs)
 ax.set_ylabel('x [m]', fontsize=fs)
 plt.legend()

def vy () :

 fig, ax = plt.subplots(constrained_layout=True)
 ax.plot(t, data[:, 2], color = 'black' , label = '$n_{step} = $' + f"{nsteps[0]:.0f}")
 ax.set_xlabel('t [s]', fontsize=fs)
 ax.set_ylabel('vy [ms$^{-1}$]', fontsize=fs)
 plt.legend()


Trajectoire ()
#Energie ()
#PosFin_Conv ()
#dts (True) 
#x()
#vy()
#Emec_Err ()

plt.show()

