#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cmath>
#include <fstream>
#include <sys/resource.h>
#include <vector>

//#include <funzioni.h>
using namespace std;

 // default_random_engine generator;
  // uniform_int_distribution<int> distribution(1,9);

double beta=0.1496927;
double t=1./beta;//2.2691853;
double expb=exp(2.*beta);
double expmb=exp(-2.*beta);


   	const int Lx=12+1;
   	const int Ly=2*(Lx-1);
    const int N=Lx*Ly*Ly*Ly;
    const int nmontecarlo=200;
    const int steptermal=100;
    const int stepdopoterm=nmontecarlo-steptermal;
	const double  numgiri=1;

//static bool m[Lx][Ly][Ly][Ly];
int megaspin=1;



bool**** m;



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

long int rng=-2; //-time(NULL);


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
{
  for(int j=0;j<Ly;j++)
    {
      for (int k = 0; k < Ly; k++)
      {
        for (int l = 0; l < Ly; l++)
        {
          m[0][j][k][l] =1; m[Lx-1][j][k][l] =1;
          for(int i=1;i<Lx-1;i++)
          {
          m[i][j][k][l] =randombit();
          }
        }
      }
    }
}

double mean()
{
  double media=0;
  for(int i=0;i<Lx;i++)
    {
        for(int j=0;j<Ly;j++)
        {
          for (int k = 0; k < Ly; k++)
          {
            for (int l = 0; l < Ly; l++)
            {
              media+= (m[i][j][k][l]==megaspin? 1:-1) ;
            }
          }
        }
    }
    return media/N;
}

double meanx(int i)
{
  double media=0;

        for(int j=0;j<Ly;j++)
        {
          for (int k = 0; k < Ly; k++)
          {
            for (int l = 0; l < Ly; l++)
            {
              media+= (m[i][j][k][l] ==megaspin? 1:-1) ;
            }
          }
        }

    return media/(Ly*Ly*Ly);
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


double weight(int i, int j, int k, int l)
{
	double w=1.;

    if(i<Lx-1)   w*=bondweight(m[i][j][k][l],m[i+1][j][k][l]);
    if(i>0)     w*=bondweight(m[i][j][k][l],m[i-1][j][k][l]);
    w*=bondweight(m[i][j][k][l],m[i][next(j)][k][l]);
    w*=bondweight(m[i][j][k][l] ,m[i][previous(j)][k][l]);
    w*=bondweight(m[i][j][k][l] ,m[i][j][next(k)][l]);
    w*=bondweight(m[i][j][k][l] ,m[i][j][previous(k)][l]);
    w*=bondweight(m[i][j][k][l] ,m[i][j][k][next(l)]);
    w*=bondweight(m[i][j][k][l] ,m[i][j][k][previous(l)]);


     return w;
     //cout<<i<<" "<<j<<" "<<w;
}

void printlattice()
{
        for(int i=0;i<Lx;i++)
        {
            for(int j=0;j<Ly;j++)
            {
                cout << (m[i][j][0][0]==1?"_": "D");
            }
            cout << endl ;
        }
        cout << endl;
    }
void spinflip()
    {
        int i=randomx();
        int j=randomy();
        int k=randomy();
        int l=randomy();
              //  cout<< "vicini di "<<lattice[i][j]<<" sono " << lattice[i][j+1]<<lattice[i+1][j]<<lattice[i-1][j]<< lattice[i][j-1]<<endl;
        //double en=localenergy(lattice,i,j); // campo magnetico +(lattice[i][j]-.5)*800;
        double peso=weight(i,j,k,l); //exp(2.*en/t);
        if(randomreal(& rng)<peso)
         m[i][j][k][l] =1-m[i][j][k][l] ;
            //    cout << "si è flippato" <<peso<< endl;
    }

int ndaijk(int i, int j, int k, int l){return i+j*Lx + k*Lx*Ly + l*Lx*Ly*Ly;}
int idan(int n){return n%Lx;}
int jdan(int n){return (n%(Lx*Ly))/Lx ;}
int kdan(int n){return n%(Lx*Ly*Ly)/(Lx*Ly) ;}
int ldan(int n){return n/(Lx*Ly*Ly);}



double pbond=1.-exp(-2.*beta);

void wolffarray() // matrice Lx*Ly*Lz cluster serve per controllare che sito non sia gia' stato preso. vettore cluster serve per flippare tutti assieme alla fine. vettore sitidafare e' quello che si svuota gradualmente
{


  int i0 = randomx();
  int j0 = randomy();
  int k0 = randomy();
  int l0 = randomy();
  int n0=ndaijk(i0, j0,k0,l0);
  bool segno=m[i0][j0][k0][l0];

 // int r= randomspin()-1;
 // per evitare di fare roba ricorsiva si fa un while su siti da fare. ogni terna che indica un sito viene convertita in un numero solo, che è scritto in una base in cui una cifra è in base Lx e le altre due in base Ly. poi le tre funzioni ida,jdan,kdan riconvertono questo numero in i,j,k


  vector<int> sitidafare;
  vector<int> siticluster;

     sitidafare.push_back(n0);

  int aa=0;
  while(!sitidafare.empty())
  {   //prima si determina chi va in cluster, poi si flippa tutti
    int n=sitidafare[0];
    int i=idan(n);
    int j=jdan(n);
    int k=kdan(n);
    int l=ldan(n);
    if( (i==0||i==Lx-1)&&aa==0 )
      {
        aa=1;
        megaspin=1-megaspin;
        for (int j= 0; j < Ly; j++)// indice grosso che include j, k e l
        {
          for (int k = 0; k < Ly; k++)
          {
            for (int l = 0; l < Ly; l++)
            {
              m[0][j][l][k]=1-m[0][j][l][k];
              m[Lx-1][j][l][k]=1-m[Lx-1][j][l][k];

              sitidafare.push_back(ndaijk(0,j,l,k) );
              sitidafare.push_back(ndaijk(Lx-1,j,l,k) );
            }
          }

         }

      }

      else { if(i>0&&i<Lx-1) m[i][j][k][l]=1-m[i][j][k][l];if(m[i][j][k][l]==1-segno) cout<<"cazzo"; }



        if(i>0)
          if (m[i-1][j][k][l]==segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(n-1);
          }
        if(i<Lx-1 )
          if (m[i+1][j][k][l]==segno &&randomreal(& rng) < pbond )
          {
          //if(i!=Lx-1||cluster[0][j][k]==0)  // questo e' per evitare che si raggiunga bordo contemporaneamente da due
          sitidafare.push_back(n+1);
          }

        if(i>0&&i<Lx-1 )
          if (m[i][previous(j)][k][l]==segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(ndaijk(i, previous(j),k,l));
          }
        if(i>0&&i<Lx-1 )
          if (m[i][next(j)][k][l]==segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(ndaijk(i, next(j),k,l));
          }

        if(i>0&&i<Lx-1 )
          if (m[i][j][previous(k)][l]== segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(ndaijk(i, j,previous(k),l));
          }
          if(i>0&&i<Lx-1 )
          if (m[i][j][next(k)][l]== segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(ndaijk(i, j,next(k),l));
          }

          if(i>0&&i<Lx-1 )
          if (m[i][j][k][previous(l)]== segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(ndaijk(i, j,k,previous(l)));
          }
          if(i>0&&i<Lx-1 )
          if (m[i][j][k][next(l)]== segno &&randomreal(& rng) < pbond )
          {
          sitidafare.push_back(ndaijk(i, j,k,next(l)));
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


}


double  tutto ()
{

    double mediacicli=0;
    double magnx[Lx]={0};
            randfill();
         //   printlattice(lattice); cout<< endl;
            for(int mosse=0;mosse<nmontecarlo;mosse++)
            {
              wolffarray();

              for(int ns=0;ns<N;ns++) spinflip();
           	 // cout<<mosse<<" "<<time(NULL)<<" "<<mean()<<endl;
          // cout<<localmean(lattice)<<endl;
		         if(mosse>steptermal)
		         {  //printlattice(lattice);
		          mediacicli+= abs(mean());
		          for(int i=0;i<Lx;i++) magnx[i]+=meanx(i);
		         }
     		    }
       cout<<"media tot e' "<< mediacicli/(stepdopoterm-1)<<endl;
       for(int i=0;i<Lx;i++) cout<<magnx[i]/(stepdopoterm-1)<<" , ";
	cout<<endl<<"punto centrale "<< magnx[(Lx-1)/2]/(stepdopoterm-1);

}
int main()
{
   srand (time(NULL));

  m = new bool***[Lx];

  for(int i = 0; i < Lx; i++)
  {
    m[i] = new bool**[Ly];
  }

  for(int i = 0; i < Lx; i++)
  {
    for(int j=0; j<Ly;j++)
        m[i][j] = new bool*[Ly];
  }


  for(int i = 0; i < Lx; i++)
  {
    for(int j=0; j<Ly;j++)
      for(int k=0;k<Ly;k++)
        m[i][j][k] = new bool[Ly];
  }

   tutto();



for(int i = 0; i < Lx; i++)
{
  for(int j=0; j<Ly;j++)
    for(int k=0;k<Ly;k++)
        delete[] m[i][j][k];
}

for(int i = 0; i < Lx; i++)
{
  for(int j=0; j<Ly;j++)
      delete[] m[i][j] ;
}
for(int i = 0; i < Lx; i++)
{
  delete[] m[i];
}
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
