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
  double eps;         // tolerance methode iterative
  double alpha;       // parametre pour le scheme d'Euler

// définition des variables
double tfin;         // Temps final
unsigned int nsteps; // Nombre de pas de temps
double mj;           // Masse Jupyter
double ms;           // Masse du Soleil
double ma; 			 // Masse du satellite 
double a;         // Distance Soleil_Jupyter
double Om;           // Vitesse de rotation du repère
double G_grav;       // Constante gravitationnelle
double xs;           // Position du Soleil
double xj;           // Position de Jupyter
double dist_a_s;     // Distance asteroide Soleil
double dist_a_j;     // Distance asteroide Jupyter 
double n ; // ordre de convergeance 
double jsteps ; 
double f ; 
bool adaptative ; 
bool jupyter ; 
bool Rg ; // 

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
    double Energy =  ma*(y[0]*y[0]+y[1]*y[1])/2 - G_grav * ms / dist_a_s + G_grav * mj / dist_a_j - pow(Om,2) * ( pow(y[2],2) + pow(y[3],2) )/2 ;

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
		if ( jupyter and Rg ) // si jupyter et si on veut écrire dans le fichier ( Rg = 1 => ref non-tournant)
		{
			*outputFile << t << " " << y[0] << " " << y[1] << " " \
			<< y[2]*cos(Om*t) - y[3]*sin(Om*t) << " " << y[2]*sin(Om*t) + y[3]*cos(Om*t) << " " << Energy << " " << jsteps << " " << dt  << endl; // write output on file
			last = 1;
		}
		else // si sans jupyter ou si on veut le référentiel R tournant (Rg = 0)
		{
			*outputFile << t << " " << y[0] << " " << y[1] << " " \
			<< y[2] << " " << y[3] << " " << Energy << " " << jsteps << " " << dt  << endl; // write output on file
			last = 1;
		}
    }
    else
    {
      last++;
    }
  }

    valarray<double> compute_f(valarray<double> const & y) 
    { 

	  dist_a_j = sqrt( ( y[2] - xj ) * ( y[2] - xj )  +  y[3]*y[3] ); // distance asteroide jupyter
      dist_a_s = sqrt( ( y[2] - xs ) * ( y[2] - xs )  +  y[3]*y[3] ); // distance asteroide soleil 
      
      //cout << mj << endl ; 
      //cout << Om << endl ; 
      
      valarray<double> f = y ; 
      
      f[0]      =  G_grav * ms * (xs - y[2]) / pow(dist_a_s,3) + G_grav * mj  * (xj - y[2]) / pow(dist_a_j,3) + 2*Om*y[1] + pow(Om,2)*y[2] ; 
      f[1]      =  - G_grav * ms * y[3] / pow(dist_a_s,3) - G_grav * mj * y[3] / pow(dist_a_j,3) - 2*Om*y[0] + pow(Om,2)*y[3] ; 
      f[2]      = y[0] ; 
      f[3]      = y[1] ; 

      return f ; 
    }
    
	void print (valarray<double> const & k , short n )
	{ cout << "k" << n << " : (" <<  k[0] << ',' << k[1] << ',' << k[2] << ',' << k[3] << ')' << endl ; }

/********************************************** Runge Kutta ********************************************/ 



	std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double ti, double dti) 
	{
		std::valarray<double> k1, k2, k3, k4, ynew;
		
		k1 = dti * compute_f(yold) ; 
		k2 = dti * compute_f(yold + k1/2.) ; 
		k3 = dti * compute_f(yold + k2/2.) ; 
		k4 = dti * compute_f(yold + k3) ; 
		
		ynew = yold + ( k1 + 2.*k2 + 2.*k3 + k4 )/6. ;
		return ynew ;
    }

public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      adaptative = configFile.get<double>("adaptative",adaptative);
      f = configFile.get<double>("f",f);
      tfin     = configFile.get<double>("tfin",tfin);	        // lire le temps final de simulation
      nsteps   = configFile.get<unsigned int>("nsteps",nsteps); // lire le nombre de pas de temps
      y0[0]    = configFile.get<double>("vx0",y0[0]);  // vitesse initiale selon x	    
      y0[1]    = configFile.get<double>("vy0",y0[1]);  // vitesse initiale selon y       
      y0[2]    = configFile.get<double>("x0",y0[2]);   // position initiale selon x       
      y0[3]    = configFile.get<double>("y0",y0[3]);   // position initiale selon y	    
      G_grav   = configFile.get<double>("G_grav",G_grav);           
      mj       = configFile.get<double>("mj",mj);            
      ms       = configFile.get<double>("ms",ms);        
      a        = configFile.get<double>("a",a);        
      sampling = configFile.get<unsigned int>("sampling",sampling);
      eps      = configFile.get<double>("eps", eps);
      maxit    = configFile.get<unsigned int>("maxit", maxit);
      alpha    = configFile.get<double>("alpha", alpha);
      jupyter  = configFile.get<double>("jupyter", jupyter);
      /** DONE : calculer le time step **/
      dt       = tfin / nsteps; 
      ma = configFile.get<double>("ma", ma);
      Rg = configFile.get<double>("Rg", Rg); // 1 = dans le ref du centre de masse  ; 0 dans le ref R 

      
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
      /** TODO : Ajouter une def de eps , f = 0.999 , jsteps , n (ordre de convergence n = 4) , ajouter une condition dans le config pour avec pas de temps adaptatif et sans **/ 
      
      if ( jupyter )
      {
		xs = - a * mj / (mj + ms) ; // position du soleil 
		xj = a * ms / ( ms + mj ) ; // position de Jupyter 
		Om = sqrt( G_grav * ( ms + mj ) / pow(a,3) ) ;
		y0[2] = 2.*a + xs ; // position initiale en x
		y0[1] = y0[1] - y0[2]*Om; // vitesse initiale selon y 
		cout << "Jupyter" << endl ; 
		
		if (Rg)
        { cout << "Référentiel RG" << endl ; }
		else 
		{ cout << "Référentiel R'" << endl ; }
      }
      else 
	  {
		 xs = 0 ; 
		 Om = 0 ;  
		 y0[2] = 2.*a ; 
		 mj = 0 ; 
		 cout << "Sans Jupyter : " ; 
		 cout << "Om = " << Om ; 
		 cout << " et mj = " << mj << endl ; 
	  }
      
      jsteps = 0 ; // initialise le nombre de pas pour le schéma adaptatif 
      
      t = 0.e0; // initialiser le temps
      y = y0;   // initialiser le position 
      last = 0; // initialise le parametre d'ecriture

      printOut(true); // ecrire la condition initiale
      
      if ( adaptative == true ) // Si adaptative on fait avec la méthode de pas de temps adaptatif sinon sans (voir else)
	  {	
		cout << "ADAPTATIF" << endl ; 
		cout << "eps : " << eps << endl ; 
/****************************************************** Partie Avec Temps Adaptatif	**************************************************/ 

		n = 5 ;  // Ordre de convergeance de Runge Kutta ordre 4 ( 5 d'après nos simulations )

		valarray<double> y1 ; 
		valarray<double> ytilde ; 
		valarray<double> y2 ; 
		double d ; 
	  
		while ( t < tfin ) // on continue jusqu'à avoir atteint le temps final 
		{	
			dt = min( dt , abs(tfin-t) ) ; // on sélectionne le dernier pas de temps de sorte à ce que t + dt = tfin 
			++ jsteps ; // on ajoute une itération par cycle 
		  
			y1 = RK4_do_onestep(y,t,dt) ; // on fait un pas en entier
			/** On fait deux demi-pas **/ 
			ytilde = RK4_do_onestep(y,t,dt/2.) ; /** premier demi-pas **/
			y2 = RK4_do_onestep(ytilde,t,dt/2.) ; /** deuxième demi-pas **/
		  
			d = (abs(y2-y1)).max() ; // on calcule l'erreur entre les deux chemin (rouge et bleu)
		  
			if ( d <= eps ) // si la distance est inférieure à epsilon, on accepte le pas 
			{
				y = y2 ; // on pose la nouvelle position équivalente aux deux demi-pas de temps successifs 
				t+= dt ; // on augmente le temps du pas effectué
				
				if (d == 0) // évite les divisions pas 0 
				{ dt = 5.*dt ; } // on allonge le pas de temps d'un facteur 2
				else
				{ dt = dt * pow( eps/d , 1/(n+1) ) ; } // sinon on allonge le pas de temps normalement 
			}
			else // si d > epsilon
			{
				while( d > eps ) // tant que d > epsilon on réduit le pas de temps et on recommence
				{		
					dt = f * dt * pow( eps/d , 1/(n+1) ) ; // on réduit le pas de temps, f = 0.999 permet d'éviter les boucles infinies 
					
					y1 = RK4_do_onestep(y,t,dt) ; // on fait un pas en entier
					/** On fait deux demi-pas **/ 
					ytilde = RK4_do_onestep(y,t,dt/2.) ; /** premier demi-pas **/
					y2 = RK4_do_onestep(ytilde,t,dt/2.) ; /** deuxième demi-pas **/
					
					d = (abs(y2-y1)).max() ; // on calcule l'erreur entre les deux chemin (rouge et bleu)
				}
				y = y2 ; // on pose la nouvelle position équivalente aux deux demi-pas de temps successifs 
				t += dt ; // on augmente le temps 
			}
			printOut(false); // on écrit le résultat dans le fichier 
		}
		printOut(true);	// on écrit le dernier pas de temps 
      }
      
/**************************************************** Partie Sans dt Adaptatif *******************************************/      

      else // sans les temps adaptatifs
      {
		  
	  cout << "SANS TEMPS ADAPTATIF" << endl ; 	  
      for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
		{
			y = RK4_do_onestep(y,t,dt);  // faire un pas de temps
			t += dt ; 
			printOut(false); // ecrire le pas de temps actuel
		}
		printOut(true); // ecrire le dernier pas de temps
	  }
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


