// ISING 3D CON BORDO FIXED, TUTTI ALGORITMI INCLUSO WOLFF CON MEGASPIN, PER DATI DA CONFRONTARE

#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>
#include <fstream>
#include <cmath>


/* come ising 2d ma adesso tipo blume capel,punto di reticolo vale 0 se non c'è spin, se c'è ha p possibili valori, 4 è numero critico. 
si calcolano una volta all'inizio possibili valori di coseno di differenza angoli, così dopo basta prendere elemento che si
 vuole di questo vettore. bondenergy dà solo energia tra i e j dovuta a coseno, energia totale include pure un +qualcosa se spin non è 0 in un sito */
using namespace std;


const int p=2; // numero angoli
const int Lx=16+1;        // mettere L+1, che è numero totale siti lungo x. di questi L-1 sono variabili
const int Ly=6*(Lx-1);
const int Lz=Ly;        // SE SI CAMBIA QUESTO BISOGNA CAMBIARE NEXT-PREVIOUS E PURE GENERATORE RANDOM DI SITO DA 1 A LY
const int N=Lx*Ly*Lz;   const int sitivariabili=Lz*Ly*(Lx-2);
const int nmontecarlo=5E2;
const int steptermalizz=nmontecarlo/4;
const int stepdopotermalizz=nmontecarlo-steptermalizz;
const double numgiri=1;
const double pi= 4.*atan(1);
const double D= -.655;	// piùmeno 60 è abbastanza per rendere tutto vuoto, cioè 0, (se D negativo) o tutto nonzero (se D positivo). per ising D girato al contrario
      double beta= 0.387721735; // beta critico per D=1.02 è beta=0.5637963; beta=5 è praticamente T=0, andare oltre è peggio
// per ising D critico e' -0.655, beta crit=
const int r=0;

double expD=exp(D);
double expmD=1./expD;
double v[p];
double pesi[p];
double pesim[p];		// pesi di boltzmann per vari possibili valori di angolo e loro inversi, per non ricalcolarli

double pbond[2*p][2*p];
bool cluster[Lx][Ly][Lz];

int megaspin=2; // da 1 a p;

#include <iostream>
#include <stdio.h>     
#include <stdlib.h>    
#include <time.h>      
#include <random>
#include <fstream>
#include <cmath>


std::mt19937 rng(time(NULL));
std::uniform_int_distribution<int> random3(0,2);     // genera spin
std::uniform_int_distribution<int> randombit(1,2);      // genera bit per vuoto/pieno
std::uniform_real_distribution<double> randomreal(0.0, 1.0);    // genera reale per accettare/rifiutare mossa
std::uniform_int_distribution<int> randomx(1,Lx-2);             // genera coordinata x, non sui bordi (0 e Lx-1) 
std::uniform_int_distribution<int> randomy(0,Ly-1);             // genera coordinata y e si usa pure per z 

//std::uniform_int_distribution<int> cento(1,100); 

int next(int i)         //se si vuole avere pbc anche su x serve Lx=Ly oppure definire nuovi next/previous
    {
        if (i==Ly-1) return 0;
        else return i+1;
    }
int previous(int i)
    {
        if (i==0) return Ly-1;
        else return i-1;
    }

void sitodopo(int &i,int &j, int &k)

    {
        if (i<Lx-2) 
            
                i++;
        
        else 
            {
                i=1;
                if(j<Ly-1)
                    j++;
                else
                {   j=0;
                    if(k<Lz-1) k++;
                    else k=0;
                }
            }
        
    }// si usa per far crescere indici ijk in spinflip1 e spinflip2

void randfill(int a[Lx][Ly][Lz])
    {
        for(int k=0;k<Ly;k++)
        {
            for(int j=0;j<Ly;j++)
            {   
                a[0][j][k]=megaspin; a[Lx-1][j][k]=megaspin;    
                for(int i=1;i<Lx-1;i++)
                {   
                    a[i][j][k]=random3(rng);
                }
            }
        }
    }

double magnx(int m[Lx][Ly][Lz], int i)    // diversa da prima, ora si sommano spin fregandosene di megaspin, e alla fine si fa prodotto scalare con megaspin
    { 
        double mediax=0; double mediay=0;

       
            for(int j=0;j<Ly;j++)
            {
                for(int k=0;k<Lz;k++)
                {   
                    if(m[i][j][k]!=0) mediax+=4*(m[i][j][k]-1.5)*(megaspin-1.5);
                }
            }
            mediax/=(Ly*Lz);
            return mediax;
                
    }

double mediaperising(int m[Lx][Ly][Lz])      // uguale a mean a parte un fattore N
    {
        double mediax=0; double mediay=0;

       
            for(int j=0;j<Ly;j++)
            {
                for(int k=0;k<Lz;k++)
                {   
                   for(int i=0;i<Lx;i++) if(m[i][j][k]!=0) mediax+=4*(m[i][j][k]-1.5)*(megaspin-1.5);
                }
            }
        
    return mediax/(N);
    }
double mediabulk(int m[Lx][Ly][Lz])      // uguale a mean a parte un fattore N
    {
        double mediax=0; double mediay=0;

       
            for(int j=0;j<Ly;j++)
            {
                for(int k=0;k<Lz;k++)
                {   
                   for(int i=2;i<Lx-2;i++) if(m[i][j][k]!=0) mediax+=4*(m[i][j][k]-1.5)*(megaspin-1.5);
                }
            }
        
    return mediax/(Lz*Ly*(Lx-4));
    }
void printpiano(int m[Lx][Ly][Lz], int k)
    {
        for(int i=0;i<Lx;i++)
        {
            for(int j=0;j<Ly;j++)
            {
                cout << m[i][j][k];
            }
            cout << endl ;
        }
        cout << endl;
    }
void printcluster(int k)
    {
        for(int i=0;i<Lx;i++)
        {
            for(int j=0;j<Ly;j++)
            {
                cout <<" "<< cluster[i][j][k];
            }
            cout << endl ;
        }
        cout << endl;
    }//printa un piano di cluster di wolff


void spinflip2(int m[Lx][Ly][Lz], int i, int j, int k)      // ergodico
    {
        int angolonuovo=random3(rng);;
        double expmenodeltaE=1;
        
        if (angolonuovo!=m[i][j][k])
        {
            if(angolonuovo!=0) 
            {       // copincollato da spinflip1, potrebbe essere sbagliato
                expmenodeltaE*=expD;  
                if(i<Lx-1)  if (m[i+1][j][k]!=0) expmenodeltaE*= pesi[ abs(m[i+1][j][k]-angolonuovo)];      
                if(i>0)     if(m[i-1][j][k]!=0)  expmenodeltaE*= pesi[ abs(m[i-1][j][k]-angolonuovo)];
                if(m[i][next(j)][k]!=0) expmenodeltaE*= pesi[ abs(m[i][next(j)][k]-angolonuovo)];
                if(m[i][previous(j)][k]!=0) expmenodeltaE*= pesi[ abs(m[i][previous(j)][k]-angolonuovo)];
                if(m[i][j][next(k)]!=0) expmenodeltaE*= pesi[ abs(m[i][j][next(k)]-angolonuovo)];
                if(m[i][j][previous(k)]!=0) expmenodeltaE*= pesi[ abs(m[i][j][previous(k)]-angolonuovo)];
               
            }

            if(m[i][j][k]!=0)
            {   expmenodeltaE*=expmD; 
                if(i<Lx-1)  if (m[i+1][j][k]!=0) expmenodeltaE*= pesim[ abs(m[i+1][j][k]-m[i][j][k]) ];
                if(i>0)     if(m[i-1][j][k]!=0)  expmenodeltaE*= pesim[ abs(m[i-1][j][k]-m[i][j][k]) ];
                if(m[i][next(j)][k]!=0) expmenodeltaE*= pesim[ abs(m[i][next(j)][k]-m[i][j][k]) ]; 
                if(m[i][previous(j)][k]!=0) expmenodeltaE*= pesim[ abs(m[i][previous(j)][k]-m[i][j][k]) ];
                if(m[i][j][next(k)]!=0) expmenodeltaE*= pesim[ abs(m[i][j][next(k)]-m[i][j][k]) ]; 
                if(m[i][j][previous(k)]!=0) expmenodeltaE*= pesim[ abs(m[i][j][previous(k)]-m[i][j][k]) ];
           
            }
        

            //cout<< "scelto spin "<<i<<"," <<j<<endl;

      
        if(randomreal(rng)<expmenodeltaE) {m[i][j][k]=angolonuovo;}// printpiano(m,k); }//cout << "spin flippato in "<< angolonuovo << endl;}
           
        }
            //    cout << "si è flippato" <<peso<< endl;
    }

int flipangolo(int a)     // prende numero da 1 a p e dà numero da 1 a p
    { 
    return 3-a;
    }

void settapesi(double be)
{
   // vettore con possibili valori coseni
    for(int jj=0;jj<p;jj++)
    {
        v[jj]=cos(2.*pi*double(jj)/p);
       
    }

        // vettore con possibili pesi boltzmann

    for(int jj=0;jj<p;jj++)
    {
        pesi[jj]=exp(be*v[jj]);
        pesim[jj]=1/pesi[jj];
        
    }

    for(int jj=0;jj<2*p;jj++)
    {
        for(int kk=0;kk<2*p;kk++)
        {
            pbond[jj][kk]=1-exp(-2* be *cos(pi*double(jj)/p)* cos(pi*double(kk)/p) );   // servono per cluster update. è complementare di p secondo hasenbusch
            
        }
    }


}

void growCluster(int m[Lx][Ly][Lz], int i, int j, int k);
void tryAdd(int m[Lx][Ly][Lz], int i, int j, int k, int spinvecchio);

void oneMonteCarloStep(int m[Lx][Ly][Lz]) 
    {
    //cout<<"cluster di roba successa prima e' "<<endl; printcluster(1);
    // no cluster defined so clear the cluster array
    for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++)
    for (int k = 0; k < Lz; k++)
        cluster[i][j][k] = 0;

    // choose a random spin and grow a cluster
    int i = randomx(rng);
    int j = randomy(rng);
    int k = randomy(rng);
    //cout << "scelto sito "<< i<< j <<endl;
    //int r= 0;//randombit(rng)-1;       // questo dev'essere tra 0 e p-1

    //cout << "i e'" << i << " j e' " << j <<" r e' "<< r << endl;


   // printpiano(m,k); cout<< endl;
    if(m[i][j][k]!=0) growCluster(m, i, j, k);
    
    
   }

void growCluster(int m[Lx][Ly][Lz],int i, int j, int k) {
     
    // mark the spin as belonging to the cluster and flip it
    cluster[i][j][k] = 1;
    //cout<< "aggiunto spin " << i<< " "<<j<<endl;
    int spinvecchio=m[i][j][k];
   
    m[i][j][k]=flipangolo(m[i][j][k]);    
   
    // if the neighbor spin does not belong to the
    // cluster, then try to add it to the cluster
    if(i>0) if(!cluster[i-1][j][k])
        tryAdd(m, i-1, j,k, spinvecchio);
    if(i<Lx-1) if(!cluster[i+1][j][k])
        tryAdd(m, i+1, j,k, spinvecchio);
    if (!cluster[i][previous(j)][k])
        tryAdd(m, i,previous(j),k, spinvecchio);
    if (!cluster[i][next(j)][k])
        tryAdd(m, i, next(j),k, spinvecchio);
    if (!cluster[i][j][previous(k)])
        tryAdd(m,i, j,previous(k), spinvecchio);
    if (!cluster[i][j][next(k)])
        tryAdd(m, i, j,next(k), spinvecchio);
    
}

void megaflip(int m[Lx][Ly][Lz])
{
    int megaspinvecchio=megaspin;       // salva valore di megaspin prima di flipparlo
    megaspin=flipangolo(megaspin);    // lo flippa
    for(int j=0;j<Ly;j++)
        for(int k=0;k<Lz;k++)
        {
            cluster[Lx-1][j][k] = 1;
            cluster[0][j][k] = 1;
            m[0][j][k]=megaspin;
            m[Lx-1][j][k]=megaspin;
        }

       // cout<< "megaspin e' diventato "<<megaspin<<endl;

       for(int j=0;j<Ly;j++)
        for(int k=0;k<Lz;k++)
        {  
            if(!cluster[1][j][k])
            tryAdd(m, 1, j,k, megaspinvecchio);
            if(!cluster[Lx-2][j][k])
            tryAdd(m, Lx-2, j,k, megaspinvecchio);
        
        }

}


void tryAdd(int m[Lx][Ly][Lz], int i, int j, int k, int spinvecchio)
{
     if(m[i][j][k]!=0)
        if (randomreal(rng) <pbond[abs(2*(m[i][j][k]-1))][abs(2*(spinvecchio-1))])//pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*spinvecchio-r)])     //perchè r arriva fino a pi, altri angoli fino 2pi
            {if(i>0&&i<Lx-1)growCluster(m,i, j, k); else megaflip(m);}
        
}


void  tutto ()
{
    int m[Lx][Ly][Lz];

    double mediafinale=0;   // magnetizzazione complessiva
    double magnxfinale[Lx]={0};   // magnetizz in funzione della distanza dal bordo
   
    for(int iciclogrosso=0;iciclogrosso<numgiri;iciclogrosso++)
        {
        	double mediatemp=0;        
            double magnxtemp[Lx]={0};

            
            //double s4temp=0;
            randfill(m);
            //printpiano(m,1);
            int i=1; int j=0; int k=0;
   			
 //ofstream nomequalsiasi; nomequalsiasi.open("dati.txt");
            for(int mosse=0;mosse<nmontecarlo;mosse++)
            {
            
            //if(mosse%500==0) cout<<mean(lattice);
                for(int sweep=0;sweep<sitivariabili;sweep++)
                    {
                        spinflip2(m,i,j,k);
                        sitodopo(i,j,k);
                    }

                for(int rr=0;rr<6;rr++)
                {
                                        
                     oneMonteCarloStep(m); //printcluster(1); // questo e' cluster
                
                }
                
                    if(mosse>steptermalizz) 
                        {
                            for(int ii=1;ii<Lx/2+1;ii++) 
                                {
                                 magnxtemp[ii]+=(magnx(m,ii)+magnx(m,Lx-1-ii))/2;
                                }
                              
                            //  mediatemp+= mediabulk(m);
                            //  s4temp+=s4bulk(m);
                               
                               
                        }/*  }/* printpiano(m,1);*/
                     
            
               
        
            }       // finito un giro di montecarlo, roba che segue e' aggiungere a medie finali quelle di questo giro
            
            
            //mediatemp/=(stepdopotermalizz-1);
           // s4temp/=(stepdopotermalizz-1);
           // mediafinale+=mediatemp;
            //binder+=s4temp/pow(mediatemp,2);
                
            for(int u=1;u<Lx/2+1;u++) { magnxfinale[u]+=magnxtemp[u]/(stepdopotermalizz-1);}
        
    /*        ofstream nomequalsiasi;
    nomequalsiasi.open("dati.txt"); */
    /*for(int jj=0;jj<Lx;jj++)
        { 
      cout<<"{"<<jj <<" , "<<mediatemp[jj] <<" },"<< endl;//  nomequalsiasi <<"{"<<j <<" , "<<mediatemp[j] <<" },"<< endl;
        }*/
           
    //nomequalsiasi.close();

            //cout <<binder<<endl;
            //mediafinale+=binder;
        }
    
    
    for(int u=1;u<Lx/2+1;u++)cout<<magnxfinale[u]/numgiri<<" , ";
    //cout<<mediafinale/numgiri ;// <<endl;
   
   
}
    






int main()
{

    srand (time(NULL));
    settapesi(beta);
    tutto();

    /*for(int tt=0;tt<20;tt++)
        {
            double bet=.28+double(tt)/100;
            settapesi(bet);
           
           cout<<endl<<"{ "<<bet<<", "; tutto();cout<< "},";
            
        }
    */

	

/*    ofstream nomequalsiasi;
    nomequalsiasi.open("dati.txt");
    for(double i=1;i<8;i++)
        { beta+=.04; settapesi();
        nomequalsiasi <<"{"<<beta <<" , "<< tutto()<<" },"<< endl;
        }
    nomequalsiasi.close();
*/
   

return 0;



}