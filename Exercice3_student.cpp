#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
                          // Fichier .tpp car inclut fonctions template
#include <numeric>

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private:

// EngineEuler specific members
  unsigned int maxit; // nombre maximale d iterations
  double tol;         // tolerance methode iterative
  double alpha;       // parametre pour le scheme d'Euler

// définition des variables
double tfin;         // Temps final
unsigned int nsteps; // Nombre de pas de temps
double ml;           // Masse de la Lune
double mt;           // Masse de la Terre
double ma; 			 // Masse du satellite 
double dist;         // Distance Soleil_Jupyter
double Om;           // Vitesse de rotation du repère
double G_grav;       // Constante gravitationnelle
double xt;           // Position du Soleil
double xl;           // Position de Jupyter
double dist_s_t;     // Distance satellite-Terre
double dist_s_l;     // Distance satellite-Lune
double n ; // ordre de convergeance 
double jsteps ; 
double f ; 

  valarray<double> y0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> y  = std::valarray<double>(0.e0, 4); // Correctly initialized

  double t,dt;  // Temps courant pas de temps

  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  
  
  
  
  void printOut(bool write)
  {
    double Energy =  (y[0]*y[0]+y[1]*y[1])/2 - G_grav * mt / dist_s_t + G_grav * ml / dist_s_l - pow(Om,2) * ( pow(y[2],2) + pow(y[3],2) )/2 ;

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      *outputFile << t << " " << y[0] << " " << y[1] << " " \
      << y[2] << " " << y[3] << " " << Energy << " "<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

    void compute_f(valarray<double>& f) // Ne pas oublier de diviser par la masse de l'astéroide pour le cas sans jupyter 
    {
	  dist_s_l = sqrt( ( y[2] - xl ) * ( y[2] - xl )  +  y[3]*y[3] );
      dist_s_t = sqrt( ( y[2] - xt ) * ( y[2] - xt )  +  y[3]*y[3] );
      
      f[0]      =  - G_grav * mt * (y[2] - xt) / pow(dist_s_t,3) + G_grav * ml * (xl - y[2]) / pow(dist_s_l,3) + 2*Om*y[1] + pow(Om,2)*y[2] ; 
      f[1]      =  - G_grav * mt * y[3] / pow(dist_s_t,3) - G_grav * ml * y[3] / pow(dist_s_l,3) - 2*Om*y[0] + pow(Om,2)*y[3] ; 
      f[2]      = y[0] ; 
      f[3]      = y[1] ; 
      
    }

/********************************************** Runge Kutta ********************************************/ 



	std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double ti, double dt) 
	{
		std::valarray<double> k1, k2, k3, k4, ynew;
		
		k1 = dt * f(yold , ti) ; 
		k2 = dt * f(yold + k1/2 , ti + dt/2) ; 
		k3 = dt * f(yold + k2/2 , ti + dt/2) ; 
		k4 = dt * f(yold + k3, ti + dt) ; 
		
		ynew = yold + ( k1 + 2*k2 + 2*k3 + k4 )/6;
		return ynew;
    }

public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      tfin     = configFile.get<double>("tfin",tfin);	        // lire le temps final de simulation
      nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire le nombre de pas de temps
      y0[0]    = configFile.get<double>("vx0",y0[0]);  // vitesse initiale selon x	    
      y0[1]    = configFile.get<double>("vy0",y0[1]);  // vitesse initiale selon y       
      y0[2]    = configFile.get<double>("x0",y0[2]);   // position initiale selon x       
      y0[3]    = configFile.get<double>("y0",y0[3]);   // position initiale selon y	    
      G_grav   = configFile.get<double>("G_grav",G_grav);           
      ml       = configFile.get<double>("ml",ml);            
      mt       = configFile.get<double>("mt",mt);        
      dist     = configFile.get<double>("dist",dist);        
      sampling = configFile.get<unsigned int>("sampling",sampling);
      tol      = configFile.get<double>("tol", tol);
      maxit    = configFile.get<unsigned int>("maxit", maxit);
      alpha    = configFile.get<double>("alpha", alpha);
      /** DONE : calculer le time step **/
      dt       = tfin / nsteps; 

      
      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
      outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
    };


    // Destructeur virtuel
    virtual ~Engine()
    {
      outputFile->close();
      delete outputFile;
    };
      // Simulation complete
    void run()
    {
      /** TODO : Ajouter une def de eps , f = 0.999 , jsteps , n (ordre de convergence n = 4)**/ 
      
      xt = - dist * ml / (ml + mt) ;
      xl = dist * mt / ( mt + ml ) ; 
      Om = sqrt( G_grav * ( mt + ml ) / pow(dist,3) ) ;

      y0[2] = xl * ( ( 1 - ml/mt) / ( 1 + sqrt(ml/mt) ) ) ;

      t = 0.e0; // initialiser le temps
      y = y0;   // initialiser le position 
      last = 0; // initialise le parametre d'ecriture

      printOut(true); // ecrire la condition initiale

	  while ( t < tfin ) 
	  {
		  dt = min( dt , tfin-dt ) ; 
		  jsteps += 1 ; 
		  
		  y1 = RK4_do_onestep(y,t,dt) ; 
		  ytilde = RK4_do_onestep(y,t,dt/2) ; 
		  y2 = RK4_do_onestep(ytilde,t,dt/2) ;
		  
		  d = (y2-y1).max() ; 
		  
		  if ( d <= eps ) 
		  {
			  y = y2 ; 
			  t += dt ; 
			  dt = dt * pow( e/d , 1/(n+1) ) ; 
		  }
		  else 
		  {
			  while( d > eps )
			  {
				  dt = f * dt * pow( e/d , 1/(n+1) ) ; 
				  y1 = RK4_do_onestep(y,t,dt) ; 
				  ytilde = RK4_do_onestep(y,t,dt/2) ; 
				  y2 = RK4_do_onestep(ytilde,t,dt/2) ;
			  }
			  y = y2 ; 
			  t += dt ; 
		  }
		   
		  printOut(false);
	  }	

      for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
      {
        step();  // faire un pas de temps
        printOut(false); // ecrire le pas de temps actuel
      }
      printOut(true); // ecrire le dernier pas de temps

    };
   
};

// programme
int main(int argc, char* argv[])
{
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

  Engine* engine;

  // Create an instance of Engine instead of EngineEuler
  engine = new Engine(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}


