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

nsteps = np.array([25 , 30 , 35 , 40 , 50 , 60 ])*1e3 # valeurs utilisées pour la convergeance 

#nsteps = np.array([40000])

#eps = np.array([])

eps = np.linspace(100,10e5,100)

#eps = np.array([100,200,300,400,700,800])

#eps = np.array([10e8,10e7,10e6,10e5,10e4,10e3,10e2,10])

#eps = np.array([10])

#nsteps = np.array([1000])

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
    convergence_list = np.array([])
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

conv_list_y = []

if adaptative : # schéma à pas de temps adaptatif 
    
    for i in range(nsimul2):  # Iterate through the results of all simulations
        #print(outputs[i])
        data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
        t = data[:, 0]
        vx = data[-1, 1]  # final position, velocity, energy
        vy = data[-1, 2]
        xx = data[-1, 3]
        yy = data[-1, 4]
        En = data[-1, 5]
        jsteps = data[-1,6] # nombre de pas de temps total pour la simulation 

        conv_list_y.append(yy)
        convergence_list = np.append(convergence_list,xx)
        jsteps_list[i] = jsteps
        Eerr = np.append( Eerr , abs(data[1,5] - data[-1,5]) ) # erreur sur l'Emec : différence entre la valeur initiale et la valeur finale 

else : # schéma à pas de temps fixe 

    for i in range(nsimul):
        data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
        #print(outputs[i])
        t = data[:, 0]
        vx = data[-1, 1]  # final position, velocity, energy
        vy = data[-1, 2]
        xx = data[-1, 3]
        yy = data[-1, 4]
        En = data[-1, 5]

        conv_list_y.append(yy)
        convergence_list.append(xx)
        Eerr = np.append( Eerr , abs(data[1,5] - data[-1,5]) ) # erreur sur l'Emec : différence entre la valeur initiale et la valeur finale

lw = 1.5
fs = 16

# ---------------------------------------- Zones Plots ---------------------------------------- #

print(" ")
# Vitesse max et min en y 
print( f" max(vy) = {max(abs(data[:,2]))}" )
print( f" min(vy) = {min(abs(data[:,2]))}" )
# Position max et min en x
print( f" max(x) = {max((data[:,3])):.3e}" )
print( f" min(x) = {min((data[:,3])):.3e}" )
print(" ")
# Erreur Emec
if adaptative :
    print( f" Emec err = {Eerr[0]} pour eps = {eps[0]} => jsteps = {jsteps_list[0]}" ) # On prend le premier epsilon (le plus précis si ordonné par odre croissant)
else :
    print( f" Emec err = {Eerr[-1]} pour nsteps = {nsteps[-1]}" ) # On prend le nombre de pas de temps max, donc (le plus précis si ordonné par ordre croissant)
print(" ")    

def Trajectoire () :

    plt.figure()
    lbl = '$n_{step} = $' + f"{nsteps[-1]:.0f}"
    if adaptative :
        lbl = '$\\epsilon = $' + f"{eps[-1]:.0f}"
        
    
    plt.plot(data[:, 3], data[:, 4], color = 'black' , label = lbl)
    Soleil  = plt.scatter(  - a * mj / (mj + ms) , 0 , marker = '*' , color = 'gold' , label = 'Soleil' )

    if jupyter : # On affiche Jupyter
        
        referentiel = "$R'$"
        if Rg :
            referentiel = "$R_G$"
        plt.title(f"Référentiel {referentiel}",fontsize = fs - 2)
        Jupyter = plt.scatter( a * ms / ( ms + mj ) ,0,marker = 'o' , color = 'brown' , label = "Jupyter" )
        
    else : # On affiche la position minimale et maximale en x

        print("")
        #xmax = plt.scatter( max(data[:,3]),-1.51e12, marker = 's' , color = 'red' , label = "$x_{max} =$" + f"{max((data[:,3])):.2e}" )
        #xmin = plt.scatter( min(data[:,3]),0, marker = '^' , color = 'red' , label = "$x_{min} =$" + f"{min((data[:,3])):.2e}" )
        
    #posinit = plt.scatter(data[0,3],data[0,4], marker = 'o' , color = 'grey' , label = "astéroide")
    
    plt.xlabel('x [m]', fontsize=fs)
    plt.ylabel('y [m]', fontsize=fs)
    plt.legend()

def Energie () : # Energie en fonction du temps

    lbl = '$n_{step} = $' + f"{nsteps[-1]:.0f}"
    if adaptative :
        lbl = '$\\epsilon = $' + f"{eps[-1]:.0f}"

    plt.figure()
    plt.plot(data[:, 0], data[:, 5], color = 'black' , label = lbl)
    plt.xlabel('t [s]', fontsize=fs)
    plt.ylabel('$E_{mec}$', fontsize=fs)
    plt.legend()

def Emec_Err ( norder = 5 ) : # Erreur de l'Emec en fonction du temps ( fixe et sans jupyter : ordre 5 ; adaptatif : ordre ? )

    if adaptative : # Erreur avec pas de temps adaptatif

        plt.figure()
        plt.loglog( 1/jsteps_list , Eerr,'k+-',linewidth = lw)
        plt.loglog( 1/jsteps_list , (1/pow(jsteps_list,norder))*1e14 , color  = 'red' , label = f"$1/N^{norder}$" , linestyle = 'dashed')
        plt.xlabel("$j_{steps}$", fontsize=fs)
        plt.ylabel('$\\delta_{E_{mec}}$', fontsize=fs)
        plt.legend(fontsize = fs - 3)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.plot()

    else : # Erreur avec pas de temps fixe  
        
        plt.figure()
        plt.loglog( dt,Eerr,'k+-',linewidth = lw)
        plt.loglog( dt , pow(dt,norder)/1e20 , color  = 'red' , label = f"$1/N^{norder}$" , linestyle = 'dashed')
        plt.xlabel(f"$\\Delta t$ [s]", fontsize=fs)
        plt.ylabel('$\\delta_{E_{mec}}$', fontsize=fs)
        plt.legend(fontsize = fs - 3)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.plot()

def Jsteps_graph() :

    plt.figure()
    plt.plot( eps , jsteps_list , 'k+-' )
    plt.xlabel( '$\\epsilon$' , fontsize = fs )
    plt.ylabel( '$j_{steps}$' , fontsize = fs )

def PosFin_Plot () :

    plt.figure()
    plt.plot(convergence_list,conv_list_y, 'k+-')
    plt.ylabel('$y_{final}$', fontsize=fs)
    plt.xlabel('$x_{final}$', fontsize=fs)

def PosFin_Conv ( norder = 5 ) : # convergeance sur la postion finale ( en x )

    if adaptative : # convergeance en fonction de jsteps ( Runge Kutta adaptatif : ordre 1 ? )

        plt.figure()
        plt.plot( pow(1/jsteps_list,norder) , convergence_list , 'k+-' , linewidth = lw ) 
        #plt.loglog(1/jsteps_list, convergence_list)
        plt.xlabel('$(1/j_{steps})$'+f"$^{norder}$", fontsize=fs)
        plt.ylabel('$x_{final}$ [m]', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)

    else : # convergeance en fonction de dt ( Runge Kutta fixe : ordre 4)

        plt.figure()
        plt.plot(pow(dt,norder), convergence_list, 'k+-', linewidth=lw)
        plt.xlabel(f"$(\\Delta t)^{norder}$ [s]", fontsize=fs)
        plt.ylabel('$x_{final}$ [m]', fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)

def dts ( jstep ) : # Pas de temps au cours du temps et jsteps au cours du temps si adaptatif et jstep = True

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

#Trajectoire ()
#Energie ()
PosFin_Conv ()
#PosFin_Plot ()
#Jsteps_graph ()
#dts (True) 
#x()
#vy()
#Emec_Err ()

plt.show()

def test():

    nsteps25000 = np.loadtxt('nsteps=25000.0.out')[-1,3]
    nsteps30000 = np.loadtxt('nsteps=30000.0.out')[-1,3]
    nsteps35000 = np.loadtxt('nsteps=35000.0.out')[-1,3]
    nsteps40000 = np.loadtxt('nsteps=40000.0.out')[-1,3]
    nsteps50000 = np.loadtxt('nsteps=50000.0.out')[-1,3]
    nsteps60000 = np.loadtxt('nsteps=60000.0.out')[-1,3]

    nsteps25000Err = abs (np.loadtxt('nsteps=25000.0.out')[0,5] - np.loadtxt('nsteps=25000.0.out')[-1,5])
    nsteps30000Err = abs (np.loadtxt('nsteps=30000.0.out')[0,5] - np.loadtxt('nsteps=30000.0.out')[-1,5])
    nsteps35000Err = abs (np.loadtxt('nsteps=35000.0.out')[0,5] - np.loadtxt('nsteps=35000.0.out')[-1,5])
    nsteps40000Err = abs (np.loadtxt('nsteps=40000.0.out')[0,5] - np.loadtxt('nsteps=40000.0.out')[-1,5])
    nsteps50000Err = abs (np.loadtxt('nsteps=50000.0.out')[0,5] - np.loadtxt('nsteps=50000.0.out')[-1,5])
    nsteps60000Err = abs (np.loadtxt('nsteps=60000.0.out')[0,5] - np.loadtxt('nsteps=60000.0.out')[-1,5])

    n = np.array([nsteps25000,nsteps30000,nsteps35000,nsteps40000,nsteps50000,nsteps60000])
    nerr = np.array([nsteps25000,nsteps30000,nsteps35000,nsteps40000,nsteps50000,nsteps60000])

    plt.figure()
    plt.loglog( pow(jsteps_list,4), convergence_list , 'k+-' , label = 'Adaptative' , linewidth = lw) # jsteps_list**norder)
    plt.loglog( pow(nsteps,4), n , 'r+-' , label = 'Fix' , linewidth = lw) # jsteps_list**norder)
    plt.xlabel('$(n_{steps})$', fontsize=fs)
    plt.ylabel('$x_{final}$', fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.legend()

    plt.figure()
    plt.loglog( pow(jsteps_list, 4), Eerr,'r+-',linewidth = lw , label = 'Adaptative')
    plt.loglog( pow(nsteps, 4) , nerr ,'k+-',linewidth = lw , label = 'Fix')
    #plt.loglog( jsteps_list , (1/pow(jsteps_list,norder))*1e14 , color  = 'red' , label = f"$N^{norder}$" , linestyle = 'dashed')
    plt.xlabel("$j_{steps}$", fontsize=fs)
    plt.ylabel('$\\delta_{E_{mec}}$', fontsize=fs)
    plt.legend(fontsize = fs - 3)
    plt.plot()

    plt.show()


#test()



