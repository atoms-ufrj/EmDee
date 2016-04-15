#include "EmDee.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "cblas.h"

#define println printf("%d\n",__LINE__)

#define nbcells 58

/* -------------------------------------------------------------------------------------------------
                                         LOCAL FUNCTIONS
------------------------------------------------------------------------------------------------- */

double max_appoach_sq( int N, double *R, double *R0 )
{
  double dx = R[0] - R0[0];
  double dy = R[1] - R0[1];
  double dz = R[2] - R0[2];
  double maximum = dx*dx + dy*dy + dz*dz;
  double next;
  int offset = 3;
  for (int i = 1; i < N; i++) {
    dx = R[offset] - R0[offset++];
    dy = R[offset] - R0[offset++];
    dz = R[offset] - R0[offset++];
    double value = dx*dx + dy*dy + dz*dz;
    if (value > maximum) {
      next = maximum;
      maximum = value;
    }
  }
  return maximum + 2*sqrt(maximum*next) + next;
}

// -------------------------------------------------------------------------------------------------

void make_cells( tEmDee *me, int M )
{
  int MM = M*M;
  me->mcells = M;
  me->ncells = M*MM;

  if (me->ncells > me->maxcells) {
    free( me->cell );
    me->cell = malloc( me->ncells*sizeof(tCell) );
    me->maxcells = me->ncells;
  }

  #define pbc( x ) if (x < 0) x += M; else if (x >= M) x -= M;
  for (int k = 0; k < nbcells; k++)
    for (int iz = 0; iz < M; iz++) {
      int jz = (iz + nb[k][2]) % M;
      pbc( jz );
      for (int iy = 0; iy < M; iy++) {
        int jy = (iy + nb[k][1]) % M;
        pbc( jy );
        for (int ix = 0; ix < M; ix++) {
          int jx = ix + nb[k][0];
          pbc( jx );
          me->cell[ix + iy*M + iz*MM].neighbor[k] = jx + jy*M + jz*MM;
        }
      }
    }
  #undef pbc
}

// -------------------------------------------------------------------------------------------------

void find_pairs( tEmDee *me, double L )
{
  int M = me->mcells;
  int MM = M*M;
  double invL = 1.0/L;
  double RcSq = me->xRcSq*invL*invL;

  // Distribute atoms over cells and save scaled coordinates:
  int *natoms = alloca( me->ncells*sizeof(int) );
  int *head = alloca( me->ncells*sizeof(int) );
  int *next = alloca( me->natoms*sizeof(int) );
  double *R = alloca( me->nx3*sizeof(double) );
  for (int icell = 0; icell < me->ncells; icell++) {
    head[icell] = -1;
    natoms[icell] = 0;
  }
  int offset = 0;
  for (int i = 0; i < me->natoms; i++) {
    double x = me->R[offset++]*invL;
    x -= floor(x);
    R[offset] = x;
    double y = me->R[offset++]*invL;
    y -= floor(y);
    R[offset] = y;
    double z = me->R[offset++]*invL;
    z -= floor(z);
    R[offset] = z;
    int icell = (int)(M*x) + M*(int)(M*y) + MM*(int)(M*z);
    next[i] = head[icell];
    head[icell] = i;
    natoms[icell]++;
  }

  // Safely allocate local arrays:
  int nmax = natoms[0];
  for (int icell = 1; icell < me->ncells; icell++)
    if (natoms[icell] > nmax)
      nmax = natoms[icell];
  int maxpairs = (nmax*((2*nbcells + 1)*nmax - 1))/2;
  int *atom = alloca( nmax*(nbcells + 1)*sizeof(int) );
  double *Rc = alloca( nmax*(nbcells + 1)*3*sizeof(double) );

  // Sweep all cells to search for neighbors:
  int npairs = 0;
  for (int icell = 0; icell < me->ncells; icell++) {
    int nlocal = natoms[icell];
    if (nlocal != 0) {

      if (me->maxpairs < npairs + maxpairs) {
        me->maxpairs = npairs + maxpairs + extraPairs;
        me->neighbor = realloc( me->neighbor, me->maxpairs*sizeof(int) );
        if (me->neighbor == NULL) {
          fprintf(stderr, "ERROR: not enough memory to allocate neighbor list.\n");
          exit(EXIT_FAILURE);
        }
      }

      // Build list of atoms in current cell and neighbor cells:
      int ntotal = 0;
      int n = 0;
      int j = head[icell];
      while (j != -1) {
        atom[ntotal++] = j;
        int jx3 = 3*j;
        Rc[n++] = R[jx3];
        Rc[n++] = R[jx3+1];
        Rc[n++] = R[jx3+2];
        j = next[j];
      }
      for (int k = 0; k < nbcells; k++) {
        j = head[me->cell[icell].neighbor[k]];
        while (j != -1) {
          atom[ntotal++] = j;
          int jx3 = 3*j;
          Rc[n++] = R[jx3];
          Rc[n++] = R[jx3+1];
          Rc[n++] = R[jx3+2];
          j = next[j];
        }
      }

      // Search for neighbors and add them to the list:
      for (int k = 0; k < nlocal; k++) {
        int i = atom[k];
        me->first[i] = npairs;
        int kx3 = 3*k;
        double Rix = Rc[kx3];
        double Riy = Rc[kx3+1];
        double Riz = Rc[kx3+2];
        for (int m = k+1; m < ntotal; m++) {
          int mx3 = 3*m;
          double dx = Rix - Rc[mx3];
          double dy = Riy - Rc[mx3+1];
          double dz = Riz - Rc[mx3+2];
          dx -= rint(dx);
          dy -= rint(dy);
          dz -= rint(dz);
          if (dx*dx + dy*dy + dz*dz <= RcSq)
            me->neighbor[npairs++] = atom[m];
        }
        me->last[i] = npairs - 1;
      }
    }
  }
  me->npairs = npairs;
}

/* -------------------------------------------------------------------------------------------------
                                         LIBRARY FUNCTIONS
------------------------------------------------------------------------------------------------- */

void md_initialize( tEmDee *me, double rc, double skin, int atoms, int *types )
{
  me->RcSq = rc*rc;
  me->xRc = rc + skin;
  me->xRcSq = me->xRc * me->xRc;
  me->skinSq = skin*skin;

  me->natoms = atoms;
  me->nx3 = 3*atoms;

  me->next = malloc( atoms*sizeof(int) );
  me->first = malloc( atoms*sizeof(int) );
  me->last = malloc( atoms*sizeof(int) );
  me->R = malloc( me->nx3*sizeof(double) );
  me->P = malloc( me->nx3*sizeof(double) );
  me->R0 = calloc( me->nx3, sizeof(double) );
  if (types) {
    me->type = malloc( atoms*sizeof(int) );
    for (int i = 0; i < atoms; i++)
      me->type[i] = types[i]-1;
  }
  else
    me->type = calloc( atoms, sizeof(int) );

  me->mcells = me->ncells = me->maxcells = 0;
  me->cell = malloc( 0 );

  me->maxpairs = extraPairs;


//  me->maxpairs = 4000000;

  me->neighbor = malloc( me->maxpairs );

  me->builds = 0;
}

// -------------------------------------------------------------------------------------------------

void md_upload( tEmDee *me, double *coords, double *momenta )
{
  if (coords)  cblas_dcopy( me->nx3, coords,  1, me->R, 1);
  if (momenta) cblas_dcopy( me->nx3, momenta, 1, me->P, 1);
//  if (coords)  memcpy( me->R, coords,  me->nx3*sizeof(double) );
//  if (momenta) memcpy( me->P, momenta, me->nx3*sizeof(double) );
}

// -------------------------------------------------------------------------------------------------

void md_download( tEmDee *me, double *coords, double *momenta, double *forces )
{
  if (coords)  cblas_dcopy( me->nx3, me->R, 1, coords,  1);
  if (momenta) cblas_dcopy( me->nx3, me->P, 1, momenta, 1);
  if (forces)  cblas_dcopy( me->nx3, me->F, 1, forces,  1);
//  if (coords)  memcpy( coords,  me->R, me->nx3*sizeof(double) );
//  if (momenta) memcpy( momenta, me->P, me->nx3*sizeof(double) );
//  if (forces)  memcpy( forces,  me->F, me->nx3*sizeof(double) );
}

// -------------------------------------------------------------------------------------------------

void md_change_coordinates( tEmDee *me, double a, double b )
{
  cblas_daxpby( me->nx3, b, me->P, 1, a, me->R, 1 );
}

// -------------------------------------------------------------------------------------------------

void md_change_momenta( tEmDee *me, double a, double b )
{
  cblas_daxpby( me->nx3, b, me->F, 1, a, me->P, 1 );
}

// -------------------------------------------------------------------------------------------------

void md_handle_neighbor_list( tEmDee *me, double L )
{
//  if (max_appoach_sq( me->natoms, me->R, me->R0 ) > me->skinSq) {
    int M = (int)floor(2*L/me->xRc);
//    if (M < 5) {
//      fprintf(stderr, "ERROR: simulation box is too small.\n");
//      exit(EXIT_FAILURE);
//    }
    if (M != me->mcells) make_cells( me, M );
    find_pairs( me, L );
//    memcpy( me->R0, me->R, me->nx3*sizeof(double) );
//    cblas_dcopy( me->nx3, me->R, 1, me->R0, 1 );
//    me->builds++;
//  }
}

// -------------------------------------------------------------------------------------------------

double md_kinetic_energy( tEmDee *me, double *mass )
{
  double ke;
  int offset = 0;
  for (int i = 0; i < me->natoms; i++) {
    double px = me->P[offset++];
    double py = me->P[offset++];
    double pz = me->P[offset++];
    ke += (px*px + py*py + pz*pz)/mass[me->type[i]];
  }
  return 0.5*ke;
}

// -------------------------------------------------------------------------------------------------

int count_pairs_directly( tEmDee *me, double L )
{
  // Scale coordinates and squared cut-off distance:
  double invL = 1.0/L;
  double R[3*me->natoms];
  for (int i = 0; i < 3*me->natoms; i++)
    R[i] = me->R[i]*invL;
  double RcSq = me->xRcSq*invL*invL;

  // Check all possible pairs:
  int npairs = 0;
  for (int i = 0; i < me->natoms-1; i++) {
    int offset = 3*i;
    double Rix = R[offset];
    double Riy = R[offset+1];
    double Riz = R[offset+2];
    for (int j = i+1; j < me->natoms; j++) {
      offset = 3*j;
      double dx = Rix - R[offset];
      double dy = Riy - R[offset+1];
      double dz = Riz - R[offset+2];
      dx -= round(dx);
      dy -= round(dy);
      dz -= round(dz);
      if (dx*dx + dy*dy + dz*dz <= RcSq)
        npairs++;
    }
  }
  return npairs;
}

// -------------------------------------------------------------------------------------------------

void compute_lennard_jones( tEmDee *me, double L )
{
/*  integer :: i, j, k*/
/*  real(8) :: r2, Rij(3), Ri(3), InvR2, InvR12, InvR6, Fij(3), Eij, Wij, InvL!, S*/
  double Epot = 0.0;
  double Virial = 0.0;
  memset( me->F, 0, me->nx3*sizeof(double) );
  md_handle_neighbor_list( me, L );
  for (int i = 0; i < me->natoms; i++) {
    int ix = 3*i;
    int iy = ix + 1;
    int iz = iy + 1;
    double Rix = me->R[ix];
    double Riy = me->R[iy];
    double Riz = me->R[iz];
    double Fix = 0.0;
    double Fiy = 0.0;
    double Fiz = 0.0;
    for (int k = me->first[i]; k < me->last[i]; k++) {
      int j = me->neighbor[k];
      int jx = 3*j;
      int jy = jx + 1;
      int jz = jy + 1;
      double dx = Rix - me->R[jx];
      double dy = Riy - me->R[jy];
      double dz = Riz - me->R[jz];
      dx -= round(dx);
      dy -= round(dy);
      dz -= round(dz);
      double r2 = dx*dx + dy*dy + dz*dz;
      if (r2 <= me->RcSq) {
        double invR2 = 1.0/r2;
        double invR6 = invR2*invR2*invR2;
        double invR12 = invR6*invR6;
        double Eij = invR12 - invR6;
        double Wij = invR12 + Eij;
        Epot += Eij;
        Virial += Wij;
        double Fa = Wij*invR2;
        double Fijx = Fa*dx;
        double Fijy = Fa*dy;
        double Fijz = Fa*dz;
        Fix += Fijx;
        Fiy += Fijy;
        Fiz += Fijz;
        me->F[jx] -= Fijx;
        me->F[jy] -= Fijy;
        me->F[jz] -= Fijz;
      }
      me->F[ix] += Fix;
      me->F[iy] += Fiy;
      me->F[iz] += Fiz;
    }
  }
  Epot *= 4.0;
  Virial *= 8.0;
  cblas_dscal( me->nx3, 24.0, me->F, 1 );
}

// -------------------------------------------------------------------------------------------------

void md_scale_momenta( tEmDee *me, double scaling )
{
  for (int i = 0; i < 3*me->natoms; i++)
    me->P[i] *= scaling;
}

// -------------------------------------------------------------------------------------------------

void md_translate_momenta( tEmDee *me, double force_factor )
{
  for (int i = 0; i < 3*me->natoms; i++)
    me->P[i] += force_factor*me->F[i];
}

// -------------------------------------------------------------------------------------------------

void md_scale_coords( tEmDee *me, double scaling )
{
  for (int i = 0; i < 3*me->natoms; i++)
    me->R[i] *= scaling;
}

// -------------------------------------------------------------------------------------------------

void md_translate_coords( tEmDee *me, double momentum_factor )
{
  for (int i = 0; i < 3*me->natoms; i++)
    me->R[i] += momentum_factor*me->P[i];
}

/* -------------------------------------------------------------------------------------------------
                               F O R T R A N    B I N D I N G S
------------------------------------------------------------------------------------------------- */

void md_allocate( void **me )
{
  *me = malloc(sizeof(tEmDee));
}

// -------------------------------------------------------------------------------------------------

void md_capture_first_and_last( tEmDee *me, void **first, void **last )
{
  *first = me->first;
  *last = me->last;
}

// -------------------------------------------------------------------------------------------------

int md_capture_neighbor( tEmDee *me, void **neighbor )
{
  *neighbor = me->neighbor;
  return me->npairs;
}

