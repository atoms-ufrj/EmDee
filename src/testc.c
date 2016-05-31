#include "EmDee.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

typedef struct {
  int N, seed, Npassos, Nprop;
  double rho, L, Rc, Rs, Rc2, Temp, Press, Dt, InvL, Ws, SixWs, Ec, Dt_2;
  double *R, *V, *F;
} simpar;

#define M_PI 3.14159265358979323846

//--------------------------------------------------------------------------------------------------

double drand()   /* uniform distribution, (0..1] */
{
  return (rand() + 1.0)/(RAND_MAX + 1.0);
}

//--------------------------------------------------------------------------------------------------

double random_normal()  /* normal distribution, centered on 0, std dev 1 */
{
  return sqrt(-2*log(drand())) * cos(2*M_PI*drand());
}

//--------------------------------------------------------------------------------------------------

void read_data( simpar *par, char *filename )
{
  FILE *file;
  file = fopen(filename,"r");
  #define readline \
    if (fgets(line, sizeof(line), file)) \
    if (!fgets(line, sizeof(line), file)) { \
      printf("ERROR: could not read data."); \
      exit(0); \
    }
  double InvRc2, InvRc6, InvRc12;
  char line[256];
  if (file != NULL) {
    readline; sscanf( line, "%d", &par->N );
    readline; sscanf( line, "%lf", &par->Rc );
    readline; sscanf( line, "%lf", &par->Rs );
    readline; sscanf( line, "%d",  &par->seed );
    readline; sscanf( line, "%lf", &par->Dt );
    readline; sscanf( line, "%d",  &par->Npassos );
    readline; sscanf( line, "%d",  &par->Nprop );
    readline; sscanf( line, "%lf", &par->rho );
    readline; sscanf( line, "%lf", &par->Temp );
    readline; sscanf( line, "%lf", &par->Press );
  }
  #undef readline
  par->Rc2 = par->Rc*par->Rc;
  par->L = pow(par->N/par->rho,1.0/3.0);
  par->InvL = 1.0/par->L;
  InvRc2 = 1.0/par->Rc2;
  InvRc6 = InvRc2*InvRc2*InvRc2;
  InvRc12 = InvRc6*InvRc6;
  par->Ws = 2.0*InvRc12 - InvRc6;
  par->SixWs = 6.0*par->Ws;
  par->Ec = InvRc12 - InvRc6;
  par->Dt_2 = 0.5*par->Dt;
  int Nx3 = 3*par->N;
  par->R = malloc( Nx3*sizeof(double) );
  par->V = malloc( Nx3*sizeof(double) );
  par->F = malloc( Nx3*sizeof(double) );
}

//--------------------------------------------------------------------------------------------------

void create_configuration( simpar *par )
{
  double Vcm[3] = { 0.0, 0.0, 0.0 };
  int Nd = ceil(pow(par->N,1.0/3.0));
  int a[3];
  for (int i = 0; i < par->N; i++) {
    a[2] = i/(Nd*Nd);
    a[1] = (i - a[2]*Nd*Nd)/Nd;
    a[0] = i - a[1]*Nd - a[2]*Nd*Nd;
    for (int j = 0; j < 3; j++) {
      int k = 3*i + j;
      par->R[k] = (par->L/Nd)*a[j] + 0.5;
      par->V[k] = random_normal();
      Vcm[j] += par->V[k];
    }
  }
  for (int j = 0; j < 3; j++)
    Vcm[j] /= par->N;
  double factor = 0.0;
  for (int i = 0; i < par->N; i++)
    for (int j = 0; j < 3; j++) {
      int k = 3*i+j;
      par->V[k] -= Vcm[j];
      factor += par->V[k]*par->V[k];
    }
  factor = sqrt(par->Temp*(3*par->N - 3)/factor);
  for (int i = 0; i < 3*par->N; i++)
    par->V[i] *= factor;
}

//--------------------------------------------------------------------------------------------------

int main( int argc, char *argv[] )  {
  int threads;
  char *filename;
  if (argc == 2) {
    threads = 1;
    filename = argv[1];
  }
  else if (argc == 3) {
    threads = atoi(argv[1]);
    filename = argv[2];
  }
  else {
    printf("Usage: testc [number-of-threads] input-file\n");
    exit(EXIT_FAILURE);
  }
  simpar par;
  read_data( &par, filename );
  create_configuration( &par );
  tEmDee md;
  md_initialize( &md, threads, par.Rc, par.Rs, par.N, 1, NULL );
  tModel lj = pair_lj( 1.0, 1.0 );   
  md_set_pair( &md, 1, 1, &lj );
  md_compute_forces( &md, par.F, par.R, par.L );
  printf("%d %lf %lf\n", 0, md.Energy, md.Virial);
  clock_t start = clock();
  for (int passo = 1; passo <= par.Npassos; passo++) {
    if (passo % par.Nprop == 0) printf("%d %lf %lf\n", passo, md.Energy, md.Virial);
    for (int i = 0; i < 3*md.natoms; i++) {
      par.V[i] += par.Dt_2*par.F[i];
      par.R[i] += par.Dt*par.V[i];
    }
    md_compute_forces( &md, par.F, par.R, par.L );
    for (int i = 0; i < 3*md.natoms; i++)
      par.V[i] += par.Dt_2*par.F[i];
  }
  clock_t diff = clock() - start;
  int msec = diff * 1000 / CLOCKS_PER_SEC;
  printf("neighbor list builds = %d\n", md.builds);
  printf("pair time = %f s.\n", md.time);
  printf("excecution time = %f s.\n", msec/1000.0);
  return EXIT_SUCCESS;
}

