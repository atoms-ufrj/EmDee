#include "EmDee.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#ifdef mkl
  #include "mkl_cblas.h"
#else
  #include "cblas.h"
#endif

#define println printf("%d\n",__LINE__)

/* -------------------------------------------------------------------------------------------------
                                         LOCAL FUNCTIONS
------------------------------------------------------------------------------------------------- */
/*
double max_appoach_sq( int N, double *R, double *R0 )
{
  double dx = R[0] - R0[0];
  double dy = R[1] - R0[1];
  double dz = R[2] - R0[2];
  double maximum = dx*dx + dy*dy + dz*dz;
  double next;
  int offset = 3;
  for (int i = 1; i < N; i++) {
    dx = R[offset] - R0[offset];
    offset++;
    dy = R[offset] - R0[offset];
    offset++;
    dz = R[offset] - R0[offset++];
    offset++;
    double value = dx*dx + dy*dy + dz*dz;
    if (value > maximum) {
      next = maximum;
      maximum = value;
    }
  }
  return maximum + 2*sqrt(maximum*next) + next;
}
*/
// -------------------------------------------------------------------------------------------------

//void make_cells( tEmDee *me, int M )
//{
//  int MM = M*M;
//  me->mcells = M;
//  me->ncells = M*MM;

//  if (me->ncells > me->maxcells) {
//    free( me->cell );
//    me->cell = malloc( me->ncells*sizeof(tCell) );
//    me->maxcells = me->ncells;
//  }

//  #define pbc( x ) if (x < 0) x += M; else if (x >= M) x -= M;
//  for (int k = 0; k < nbcells; k++) {
//    int kx = nb[k][0];
//    int ky = nb[k][1];
//    int kz = nb[k][2];
//    for (int iz = 0; iz < M; iz++) {
//      int jz = iz + kx;
//      pbc( jz );
//      for (int iy = 0; iy < M; iy++) {
//        int jy = iy + ky;
//        pbc( jy );
//        for (int ix = 0; ix < M; ix++) {
//          int jx = ix + kz;
//          pbc( jx );
//          me->cell[ix + iy*M + iz*MM].neighbor[k] = jx + jy*M + jz*MM;
//        }
//      }
//    }
//  }
//  #undef pbc
//}

// -------------------------------------------------------------------------------------------------

//void find_pairs( tEmDee *me, double L )
//{
//  int M = me->mcells;
//  int MM = M*M;
//  int ntypes = me->ntypes;
//  double invL = 1.0/L;
//  double xRcSq = me->xRcSq*invL*invL;

//  // Distribute atoms over cells and save scaled coordinates:
//  double *R = alloca( me->nx3*sizeof(double) );
//  int *next = alloca( me->natoms*sizeof(int) );
//  int *natoms = alloca( me->ncells*sizeof(int) );
//  int *head = alloca( me->ncells*sizeof(int) );
//  for (int icell = 0; icell < me->ncells; icell++) {
//    head[icell] = -1;
//    natoms[icell] = 0;
//  }
//  int offset = 0;
//  for (int i = 0; i < me->natoms; i++) {
//    double x = me->R[offset]*invL;
//    x -= floor(x);
//    R[offset++] = x;
//    double y = me->R[offset]*invL;
//    y -= floor(y);
//    R[offset++] = y;
//    double z = me->R[offset]*invL;
//    z -= floor(z);
//    R[offset++] = z;
//    int icell = (int)(M*x) + M*(int)(M*y) + MM*(int)(M*z);
//    next[i] = head[icell];
//    head[icell] = i;
//    natoms[icell]++;
//  }

//  // Safely allocate local arrays:
//  int nmax = natoms[0];
//  for (int icell = 1; icell < me->ncells; icell++)
//    if (natoms[icell] > nmax)
//      nmax = natoms[icell];
//  int maxpairs = (nmax*((2*nbcells + 1)*nmax - 1))/2;
//  int *atom = alloca( nmax*(nbcells + 1)*sizeof(int) );

//  // Sweep all cells to search for neighbors:
//  int npairs = 0;
//  for (int icell = 0; icell < me->ncells; icell++) {
//    int nlocal = natoms[icell];
//    if (nlocal != 0) {

//      if (me->maxpairs < npairs + maxpairs) {
//        me->maxpairs = npairs + maxpairs;
//        me->neighbor = realloc( me->neighbor, me->maxpairs*sizeof(int) );
//        if (me->neighbor == NULL) {
//          fprintf(stderr, "ERROR: not enough memory to allocate neighbor list.\n");
//          exit(EXIT_FAILURE);
//        }
//      }

//      // Build list of atoms in current cell and neighbor cells:
//      int ntotal = 0;
//      int j = head[icell];
//      while (j != -1) {
//        atom[ntotal++] = j;
//        j = next[j];
//      }
//      for (int k = 0; k < nbcells; k++) {
//        j = head[me->cell[icell].neighbor[k]];
//        while (j != -1) {
//          atom[ntotal++] = j;
//          j = next[j];
//        }
//      }

//      // Search for neighbors and add them to the list:
//      for (int k = 0; k < nlocal; k++) {
//        int i = atom[k];
//        me->first[i] = npairs;
//        int ix3 = 3*i;
//        double Rix = R[ix3];
//        double Riy = R[++ix3];
//        double Riz = R[++ix3];
//        int nxitype = ntypes*me->type[i];
//        for (int m = k+1; m < ntotal; m++) {
//          int j = atom[m];
//          if (me->pairType[nxitype + me->type[j]].model) {
//            int jx3 = 3*j;
//            double dx = Rix - R[jx3];
//            double dy = Riy - R[++jx3];
//            double dz = Riz - R[++jx3];
//            dx -= round(dx);
//            dy -= round(dy);
//            dz -= round(dz);
//            if (dx*dx + dy*dy + dz*dz <= xRcSq)
//                me->neighbor[npairs++] = j;
//          }
//        }
//        me->last[i] = npairs - 1;
//      }
//    }
//  }
//  me->npairs = npairs;
//}

// -------------------------------------------------------------------------------------------------

void find_pairs_brute_force( tEmDee *me, double L )
{
  int ntypes = me->ntypes;
  double invL = 1.0/L;
  double xRcSq = me->xRcSq*invL*invL;

  // Distribute atoms over cells and save scaled coordinates:
  int offset = 0;
  double *R = alloca( me->nx3*sizeof(double) );
  for (int i = 0; i < me->natoms; i++) {
    double x = me->R[offset]*invL;
    R[offset++] = x;
    double y = me->R[offset]*invL;
    R[offset++] = y;
    double z = me->R[offset]*invL;
    R[offset++] = z;
  }

  int npairs = 0;

  // Search for neighbors and add them to the list:
  for (int i = 0; i < me->natoms-1; i++) {
    me->first[i] = npairs;
    int ix3 = 3*i;
    double Rix = R[ix3];
    double Riy = R[++ix3];
    double Riz = R[++ix3];
    int nxitype = ntypes*me->type[i];
    for (int j = i+1; j < me->natoms; j++)
      if (me->pairType[nxitype + me->type[j]].model) {
        int jx3 = 3*j;
        double dx = Rix - R[jx3];
        double dy = Riy - R[++jx3];
        double dz = Riz - R[++jx3];
        dx -= round(dx);
        dy -= round(dy);
        dz -= round(dz);
        if (dx*dx + dy*dy + dz*dz <= xRcSq)
npairs++;
//            me->neighbor[npairs++] = j;
      }
    me->last[i] = npairs - 1;
  }
  me->npairs = npairs;
}

// -------------------------------------------------------------------------------------------------

//void handle_neighbor_list( tEmDee *me, double L )
//{
////printf("%f %f\n",max_appoach_sq( me->natoms, me->R, me->R0 ), me->skinSq);
//  if (max_appoach_sq( me->natoms, me->R, me->R0 ) > me->skinSq) {
//    int M = floor(ndiv*L/me->xRc);
//    if (M < 5) {
//      fprintf(stderr, "ERROR: simulation box is too small.\n");
//      exit(EXIT_FAILURE);
//    }
//    if (M != me->mcells) make_cells( me, M );
////printf("find pairs\n");
//    find_pairs( me, L );
//    cblas_dcopy( me->nx3, me->R, 1, me->R0, 1 );
////    memcpy( me->R0, me->R, me->nx3*sizeof(double) );
//    me->builds++;
//  }
//}

// -------------------------------------------------------------------------------------------------

void configure_atom_types( tEmDee *me, int *types )
{
  int min, max;
  if (types) {
    me->type = malloc( me->natoms*sizeof(int) );
    min = max = types[0];
    for (int i = 0; i < me->natoms; i++) {
      int itype = types[i];
      if (itype < min) min = itype;
      if (itype > max) max = itype;
      me->type[i] = itype-1;
    }
  }
  else {
    min = max = 1;
    me->type = calloc( me->natoms, sizeof(int) );
  }

  if (min != 1) {
    fprintf(stderr, "ERROR: type identifiers must be positive integers.\n");
    exit(EXIT_FAILURE);
  }
  me->ntypes = max;
  me->pairType = malloc( max*max*sizeof(tPairType) );
  for (int i = 0; i < max*max; i++)
    me->pairType[i].model = NONE;
}

// -------------------------------------------------------------------------------------------------

tModelOutput lennard_jones( double sr2, double eps4 )
{
  double sr6 = sr2*sr2*sr2;
  double sr12 = sr6*sr6;
  tModelOutput lj;
  lj.Eij = eps4*(sr12 - sr6);
  lj.Wij = 6.0*(eps4*sr12 + lj.Eij);
  return lj;
}

// -------------------------------------------------------------------------------------------------

tModelOutput force_shifting( tModelOutput model, double r2, double rc, double Ec, double Fs )
{
  double r = sqrt(r2);
  tModelOutput result;
  result.Wij = model.Wij - Fs*r;
  result.Eij = model.Eij - Ec + Fs*(r - rc);
  return result;
}

/* -------------------------------------------------------------------------------------------------
                                         LIBRARY FUNCTIONS
------------------------------------------------------------------------------------------------- */

//void md_initialize( tEmDee *me, double rc, double skin, int atoms, int types, 
//                    int *type_index, double *mass )
//{
//  me->Rc = rc;
//  me->RcSq = rc*rc;
//  me->xRc = rc + skin;
//  me->xRcSq = (rc + skin)*(rc + skin);
//  me->skinSq = skin*skin;

//  me->natoms = atoms;
//  me->nx3 = 3*atoms;

//  me->first = malloc( atoms*sizeof(int) );
//  me->last = malloc( atoms*sizeof(int) );
//  me->R = malloc( me->nx3*sizeof(double) );
//  me->P = malloc( me->nx3*sizeof(double) );
//  me->F = malloc( me->nx3*sizeof(double) );
//  me->R0 = calloc( me->nx3, sizeof(double) );

//  me->ntypes = types;
//  me->type = malloc( me->natoms*sizeof(int) );
//  if (type_index)
//    for (int i = 0; i < me->natoms; i++)
//      me->type[i] = type_index[i] - 1;
//  else
//    memset( me->type, 0, me->natoms*sizeof(int) );

//  me->pairType = malloc( types*types*sizeof(tPairType) );
//  for (int i = 0; i < types*types; i++)
//    me->pairType[i].model = NONE;

//  me->invmass = malloc( me->nx3*sizeof(double) );
//  for (int i = 0; i < me->natoms; i++) {
//    double invmass = 1.0/mass[me->type[i]];
//    me->invmass[3*i] = me->invmass[3*i+1] = me->invmass[3*i+2] = invmass;
//  }

//  me->mcells = me->ncells = me->maxcells = me->maxpairs = me->builds = 0;
//  me->cell = malloc( 0 );
//  me->neighbor = malloc( 0 );
//}

//--------------------------------------------------------------------------------------------------

//void md_set_lj( tEmDee *me, int i, int j, double sigma, double epsilon )
//{
//  tPairType *a = &me->pairType[i-1 + me->ntypes*(j-1)];
//  tPairType *b = &me->pairType[j-1 + me->ntypes*(i-1)];
//  a->model = b->model = LJ;
//  a->p1 = b->p1 = sigma*sigma;
//  a->p2 = b->p2 = 4.0*epsilon;
//}

//--------------------------------------------------------------------------------------------------

void md_set_shifted_force_lj( tEmDee *me, int i, int j, double sigma, double epsilon )
{
  tPairType *a = &me->pairType[i-1 + me->ntypes*(j-1)];
  tPairType *b = &me->pairType[j-1 + me->ntypes*(i-1)];
  a->model = b->model = SHIFTED_FORCE_LJ;
  double eps4 = 4.0*epsilon;
  double sigbyr6 = pow(sigma/me->Rc,6);
  a->p1 = b->p1 = pow(sigma,6);
  a->p2 = b->p2 = eps4;
  a->p3 = b->p3 = eps4*sigbyr6*(sigbyr6 - 1.0);
  a->p4 = b->p4 = 6.0*eps4*sigbyr6*(2.0*sigbyr6 - 1.0)/me->Rc;
}

// -------------------------------------------------------------------------------------------------

//void md_upload( tEmDee *me, double *coords, double *momenta )
//{
//  if (coords)  cblas_dcopy( me->nx3, coords,  1, me->R, 1);
//  if (momenta) cblas_dcopy( me->nx3, momenta, 1, me->P, 1);
////  if (coords)  memcpy( me->R, coords,  me->nx3*sizeof(double) );
////  if (momenta) memcpy( me->P, momenta, me->nx3*sizeof(double) );
//}

// -------------------------------------------------------------------------------------------------

//void md_download( tEmDee *me, double *coords, double *momenta, double *forces )
//{
//  if (coords)  cblas_dcopy( me->nx3, me->R, 1, coords,  1);
//  if (momenta) cblas_dcopy( me->nx3, me->P, 1, momenta, 1);
//  if (forces)  cblas_dcopy( me->nx3, me->F, 1, forces,  1);
////  if (coords)  memcpy( coords,  me->R, me->nx3*sizeof(double) );
////  if (momenta) memcpy( momenta, me->P, me->nx3*sizeof(double) );
////  if (forces)  memcpy( forces,  me->F, me->nx3*sizeof(double) );
//}

// -------------------------------------------------------------------------------------------------

//void md_change_coordinates( tEmDee *me, double a, double b )
//{
//double *R = me->R;
//double *invmass = me->invmass;
//double *P = me->P;
//for (int i = 0; i < me->nx3; i++) {
//  R[i] = a*R[i] + b*invmass[i]*P[i];
//}
//return;

//  int N = me->nx3;
//  cblas_dgbmv( CblasRowMajor, CblasNoTrans, N, N, 0, 0, b, me->invmass, 1, me->P, 1, a, me->R, 1 );

//}

// -------------------------------------------------------------------------------------------------

//void md_change_momenta( tEmDee *me, double a, double b )
//{
//double *P = me->P;
//double *F = me->F;
//for (int i = 0; i < me->nx3; i++)
//  P[i] = a*P[i] + b*F[i];
//return;

//  cblas_daxpby( me->nx3, b, me->F, 1, a, me->P, 1 );
//}

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

void md_compute_brute_force( tEmDee *me, double L )
{
  double invL = 1.0/L;
  for (int i = 0; i < me->nx3; i++) {
    me->F[i] = 0.0;
  }
  me->Energy = me->Virial = 0.0;
int npairs = 0;
  for (int i = 0; i < me->natoms-1; i++) {
    int ix = 3*i;
    int iy = ix + 1;
    int iz = iy + 1;
    double Rix = me->R[ix];
    double Riy = me->R[iy];
    double Riz = me->R[iz];
    double Fix = 0.0;
    double Fiy = 0.0;
    double Fiz = 0.0;
    for (int j = i+1; j < me->natoms; j++) {
      int jx = 3*j;
      int jy = jx + 1;
      int jz = jy + 1;
      double dx = Rix - me->R[jx];
      double dy = Riy - me->R[jy];
      double dz = Riz - me->R[jz];
      dx -= L*rint(invL*dx);
      dy -= L*rint(invL*dy);
      dz -= L*rint(invL*dz);
      double r2 = dx*dx + dy*dy + dz*dz;

      if (r2 <= me->RcSq) {
npairs++;
        double invR2 = 1.0/r2;
        double InvR6 = invR2*invR2*invR2;
        double InvR12 = InvR6*InvR6;
        double Eij = 4.0*(InvR12 - InvR6);
        double Wij = 24.0*(2.0*InvR12 - InvR6);
        me->Energy += Eij;
        me->Virial += Wij;
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
    }
    me->F[ix] += Fix;
    me->F[iy] += Fiy;
    me->F[iz] += Fiz;
  }
  me->Virial /= 3.0;
printf("real pairs C = %d\n",npairs);
//exit(0);

//for (int i = 0; i < me->nx3; i++)
//  printf("%f\n",me->F[i]);
//exit(0);
}

// -------------------------------------------------------------------------------------------------

//void md_compute_forces( tEmDee *me, double L )
//{
//  int ix, iy, iz, jx, jy, jz, itype;
//  double Rix, Riy, Riz, Fix, Fiy, Fiz, Fijx, Fijy, Fijz, dx, dy, dz;

//  handle_neighbor_list( me, L );
////  md_compute_brute_force( me, L );
////  printf("%f %f %f\n",me->F[0],me->F[1],me->F[2]);

////  static void *label[] = { NULL, &&LJ, &&SHIFTED_FORCE_LJ };
//  double *restrict F = alloca( me->nx3*sizeof(double) );
////  double *F = me->F;
//  double Epot = 0.0;
//  double Virial = 0.0;
//  double invL = 1.0/L;
//  double invL2 = invL*invL;
//  double RcSq = me->RcSq*invL2;
//  double *restrict R = alloca( me->nx3*sizeof(double) );
//  cblas_dcopy( me->nx3, me->R, 1, R, 1 );
//  cblas_dscal( me->nx3, invL, R, 1 );
//  memset( F, 0, me->nx3*sizeof(double) );
//  int ntypes = me->ntypes;
////int npairs = 0;
//  for (int i = 0; i < me->natoms; i++) {
//    ix = 3*i;
//    iy = ix + 1;
//    iz = iy + 1;
//    itype = me->type[i];
//    Rix = R[ix];
//    Riy = R[iy];
//    Riz = R[iz];
//    Fix = 0.0;
//    Fiy = 0.0;
//    Fiz = 0.0;
//    for (int k = me->first[i]; k <= me->last[i]; k++) {
//      int j = me->neighbor[k];
//      jx = 3*j;
//      jy = jx + 1;
//      jz = jy + 1;
//      dx = Rix - R[jx];
//      dy = Riy - R[jy];
//      dz = Riz - R[jz];
//      dx -= round(dx);
//      dy -= round(dy);
//      dz -= round(dz);
//      double r2 = dx*dx + dy*dy + dz*dz;
//      if (r2 <= RcSq) {
////npairs++;
//        tPairType *pair = &me->pairType[itype + ntypes*me->type[j]];
//        double invR2 = invL2/r2;

////        double InvR6 = invR2*invR2*invR2;
////        double InvR12 = InvR6*InvR6;
////        tModelOutput result;
////        result.Eij = InvR12 - InvR6;
////        result.Wij = InvR12 + result.Eij;



////printf("params = %f %f\n",pair->p1,pair->p2);
////printf("%f %f\n",4.0*result.Eij,24.0*result.Wij);
//        tModelOutput result;
//        switch (pair->model) {
//          case LJ:
//            result = lennard_jones( invR2*pair->p1, pair->p2 );
//            break;
//          case SHIFTED_FORCE_LJ:
//            result = lennard_jones( invR2*pair->p1, pair->p2 );
//            result = force_shifting( result, r2, me->Rc, pair->p3, pair->p4 );
//            break;
//        }

////printf("%f %f\n",result.Eij,result.Wij);
////exit(0);

////        goto *label[pair->model];
////        LJ:
////          result = lennard_jones( invR2, pair->p1, pair->p2 );
////          goto END;
////        SHIFTED_FORCE_LJ:
////          result = lennard_jones( invR2, pair->p1, pair->p2 );
////          result = force_shifting( result, r2, me->Rc, pair->p3, pair->p4 );
////        END:
//        Epot += result.Eij;
//        Virial += result.Wij;
//        double Fa = result.Wij*invR2*L;
//        Fijx = Fa*dx;
//        Fijy = Fa*dy;
//        Fijz = Fa*dz;
//        Fix += Fijx;
//        Fiy += Fijy;
//        Fiz += Fijz;
//        F[jx] -= Fijx;
//        F[jy] -= Fijy;
//        F[jz] -= Fijz;
//      }
//    }
//    F[ix] += Fix;
//    F[iy] += Fiy;
//    F[iz] += Fiz;
//  }
//  me->Energy = Epot;
//  me->Virial = Virial/3.0;
////  cblas_dcopy( me->nx3, F, 1, me->F, 1 );
//  memcpy( me->F, F, me->nx3*sizeof(double) );
////  me->Energy = 4.0*Epot;
////  me->Virial = 8.0*Virial;
////  cblas_dscal( me->nx3, 24.0, F, 1 );
////printf("real pairs = %d\n",npairs);
////printf("%f %f %f\n",me->F[0],me->F[1],me->F[2]);
////exit(0);
//}

