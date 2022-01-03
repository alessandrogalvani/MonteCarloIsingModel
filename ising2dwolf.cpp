#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <random>
#include <fstream>

//#include <funzioni.h>
using namespace std;

 // default_random_engine generator;
  // uniform_int_distribution<int> distribution(1,9);

double t=2.2691853;
double beta=1./t;
double expb=exp(2./t);
double expmb=exp(-2./t);


   	const int Lx=16+1;
   	const int Ly=6*(Lx-1);
    const int N=Lx*Ly;
    const int nmontecarlo=100000;
    const int steptermal=50000;
    const int stepdopoterm=nmontecarlo-steptermal;
	const double  numgiri=1;

static bool m[Lx][Ly];
static bool cluster[Lx][Ly];
int megaspin=1;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double randomreal(long *idum) // LA SUA VARIABILE DEV'ESSERE & coso, dove coso definito all'inizio coso=-time(NULL) SERVE IL MENO
{ 
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else {return temp;}
}

long int rng=-time(NULL);


int randint(int max){return 1+(int) (double(max)*rand()/(RAND_MAX+1.0)); }
int randombit(){return randint(2)-1;}
int randomx(){return randint(Lx-2);}
int randomy(){return randint(Ly)-1;}


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

void randfill()
{ for(int j=0;j<Ly;j++)
    {
    	m[0][j]=1; m[Lx-1][j]=1;
        for(int i=1;i<Lx-1;i++)
        {
        m[i][j]=randombit();
        };
    }
}

double mean()
{ double media=0;
 for(int i=0;i<Lx;i++)
    {
        for(int j=0;j<Ly;j++)
        {
        media+= (m[i][j]==megaspin? 1:-1) ;
        }
    }
    return media/N;
}

double meanx(int i)
{ double media=0;
 
        for(int j=0;j<Ly;j++)
        {
        media+= (m[i][j]==megaspin? 1:-1) ;
        }
    
    return media/Ly;
}

int bondenergy(bool i, bool j)
    {
            if (i==j) return -1; else return 1;
    }

int localenergy(int i, int j)
    {
      int e=0;
      e+=bondenergy(m[i][j],m[i][j+1<Ly?j+1:0]);
     if(i<Lx-1)   e+=bondenergy(m[i][j],m[i+1][j]);
     if(i>0)      e+=bondenergy(m[i][j],m[i-1][j]);
        e+=bondenergy(m[i][j],m[i][j>0?j-1:Ly-1]);
     return e;
    }

double bondweight(bool i, bool j)
    {
            if (i==j) return expmb; else return expb;
    }


double weight(int i, int j)
{
	double w=1.;
      w*=bondweight(m[i][j],m[i][previous(j)]);
     if(i<Lx-1)   w*=bondweight(m[i][j],m[i+1][j]);
     if(i>0)     w*=bondweight(m[i][j],m[i-1][j]);
        w*=bondweight(m[i][j],m[i][next(j)]);
     return w;
     //cout<<i<<" "<<j<<" "<<w;
}
/* void printbonds(bool h[][L], bool v[][L]) // prende matrici di bond e stampa disegnino di reticolo coi bond
    {
        for(int i=0;i<L;i++)
        {
            for(int j=0;j<L-1;j++)
            {
                cout << "O" << (h[i][j]?"-":" ");
            }
            cout <<"O" << endl;
            for(int k=0;k<L;k++)
            {   if(i<L-1)
                cout <<  (v[i][k]?"|":" ")<< " ";
            }
            cout << endl;
        }
    }
*/
void printlattice()
{
        for(int i=0;i<Lx;i++)
        {
            for(int j=0;j<Ly;j++)
            {
                cout << (m[i][j]==1?"_": "D");
            }
            cout << endl ;
        }
        cout << endl;
    }
void spinflip()
    {
        int i=randomx();
        int j=randomy();
              //  cout<< "vicini di "<<lattice[i][j]<<" sono " << lattice[i][j+1]<<lattice[i+1][j]<<lattice[i-1][j]<< lattice[i][j-1]<<endl;
        //double en=localenergy(lattice,i,j); // campo magnetico +(lattice[i][j]-.5)*800;
        double peso=weight(i,j); //exp(2.*en/t);
        if(randomreal(& rng)<peso)
         m[i][j]=1-m[i][j];
            //    cout << "si è flippato" <<peso<< endl;
    }

int ndaijk(int i, int j){return i+j*Ly;}
int idan(int n){return n%Ly;}
int jdan(int n){return n/Ly;}

double pbond=1.-exp(-2.*beta);

void wolff() // matrice Lx*Ly*Lz cluster serve per controllare che sito non sia gia' stato preso. vettore cluster serve per flippare tutti assieme alla fine. vettore sitidafare e' quello che si svuota gradualmente
{ 
  for (int i = 0; i < Lx; i++)
    {for (int j = 0; j < Ly; j++)
      cluster[i][j] = 0;
    }

  int i0 = randomx();
  int j0 = randomy();

 // int r= randomspin()-1;
 // per evitare di fare roba ricorsiva si fa un while su siti da fare. ogni terna che indica un sito viene convertita in un numero solo, che è scritto in una base in cui una cifra è in base Lx e le altre due in base Ly. poi le tre funzioni ida,jdan,kdan riconvertono questo numero in i,j,k


  vector<int> sitidafare;
  vector<int> siticluster;

  
    
     cluster[i0][j0]=1;
     sitidafare.push_back(ndaijk(i0, j0));
     siticluster.push_back(ndaijk(i0, j0));
    
  int aa=0;
  while(!sitidafare.empty())
  {   //prima si determina chi va in cluster, poi si flippa tutti
    int i=idan(sitidafare[0]);
    int j=jdan(sitidafare[0]);
    if( (i==0||i==Lx-1)&&aa==0 )
      { aa=1;
        megaspin=1-megaspin;
        for (int jj = 0; jj < Ly; jj++)
        {
            if(!cluster[0][jj])//if(i!=0||jj!=j||kk!=k)
            {
              cluster[0][jj]=1;
              sitidafare.push_back(ndaijk(0, jj)); 
              siticluster.push_back(ndaijk(0, jj));
            }
            if(!cluster[Lx-1][jj])//if(i!=Lx-1||jj!=j||kk!=k)
            {
              cluster[Lx-1][jj]=1;
              sitidafare.push_back(ndaijk(Lx-1, jj)); 
              siticluster.push_back(ndaijk(Lx-1, jj));
            }
          
         }
        
      }
      

    
      
        if(i>0&&!cluster[i-1][j])
          if (m[i-1][j]==m[i][j]&&randomreal(& rng) < pbond )
          { 
          cluster[i-1][j]=1;
          sitidafare.push_back(ndaijk(i-1, j));
          siticluster.push_back(ndaijk(i-1, j));
          }
        if(i<Lx-1&&!cluster[i+1][j])
          if (m[i+1][j]==m[i][j]&&randomreal(& rng) < pbond )
          {
          cluster[i+1][j]=1;
          //if(i!=Lx-1||cluster[0][j][k]==0)  // questo e' per evitare che si raggiunga bordo contemporaneamente da due
          sitidafare.push_back(ndaijk(i+1, j));
          siticluster.push_back(ndaijk(i+1, j));
          }
        if(i>0&&i<Lx-1&&!cluster[i][previous(j)])
          if (m[i][previous(j)]==m[i][j]&&randomreal(& rng) < pbond )
          {
          cluster[i][previous(j)]=1;
          sitidafare.push_back(ndaijk(i, previous(j))); 
          siticluster.push_back(ndaijk(i, previous(j)));
          }

        if(i>0&&i<Lx-1&&!cluster[i][next(j)])
          if (m[i][next(j)]==m[i][j]&&randomreal(& rng) < pbond )
          {
          cluster[i][next(j)]=1;
          sitidafare.push_back(ndaijk(i, next(j)));
          siticluster.push_back(ndaijk(i, next(j)));
          }        
/*
        if(i>0&&i<Lx-1&&!cluster[i][j][next(k)]&&m[i][j][next(k)]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][j][next(k)]-1)-r) ] )
          {
          cluster[i][j][next(k)]=1;
          sitidafare.push_back(ndaijk(i, j , next(k) ) ); 
          siticluster.push_back(ndaijk(i, j , next(k) ) );
          }  

        if(i>0&&i<Lx-1&&!cluster[i][j][previous(k)]&&m[i][j][previous(k)]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][j][previous(k)]-1)-r) ] )
          {
          cluster[i][j][previous(k)]=1;
          sitidafare.push_back(ndaijk(i, j,previous(k) )); 
          siticluster.push_back(ndaijk(i, j,previous(k) ));
          }
       */   
    
      
    sitidafare.erase(sitidafare.begin());
  }// fine del while
//cout<<siticluster.size()<<endl;
  for(vector<int>::iterator q=siticluster.begin();q!=siticluster.end();q++)
  { 
    int is=idan(*q);
    int js=jdan(*q);
   //cout<<is<<" "<<js<<" "<<ks<<endl;
    m[is][js]=1-m[is][js];
  }//cout<<"taglia "<<siticluster.size()<<endl;
  
}


double  tutto ()
{	

    double mediacicli=0;
    double magnx[Lx]={0};
            randfill();
         //   printlattice(lattice); cout<< endl;
            for(int mosse=0;mosse<nmontecarlo;mosse++)
            {
              wolff();
              //printlattice();
               //for(int ns=0;ns<12*N;ns++)
                //spinflip();
           	   
          // cout<<localmean(lattice)<<endl;
		         if(mosse>steptermal)
		         {  //printlattice(lattice);
		          mediacicli+= abs(mean());
		          for(int i=0;i<Lx;i++) magnx[i]+=meanx(i);
		         }
     		}
       //cout<<"media tot e' "<< mediacicli/(stepdopoterm-1)<<endl;
       for(int i=0;i<Lx;i++) cout<<magnx[i]/(stepdopoterm-1)<<" , ";
	cout<<endl<<"punto centrale "<< magnx[(Lx-1)/2]/(stepdopoterm-1);
       
}
int main()
{
   srand (time(NULL));
   tutto();
/*ofstream nomequalsiasi;

nomequalsiasi.open("dati.txt");
for(double t=1.4;t<3.1;t=t+.1)
      {nomequalsiasi <<t<< "  "<<tutto(t) << endl;}
   nomequalsiasi.close();
*/
return 0;

/*for(double t=.00008;t<3;t=t+.3)
      {cout <<"temperatura e' "<<t<< " magn "<<tutto(t) << endl;} */

}
