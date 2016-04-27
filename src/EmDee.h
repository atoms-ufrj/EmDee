typedef struct {
  int model;
  double p1;
  double p2;
  double p3;
  double p4;
} tPairType;

#define ndiv 2
#define nbcells 62
static const int nb[nbcells][3] = { { 1, 0, 0}, { 2, 0, 0}, {-2, 1, 0}, {-1, 1, 0}, { 0, 1, 0},
                                    { 1, 1, 0}, { 2, 1, 0}, {-2, 2, 0}, {-1, 2, 0}, { 0, 2, 0},
                                    { 1, 2, 0}, { 2, 2, 0}, {-2,-2, 1}, {-1,-2, 1}, { 0,-2, 1},
                                    { 1,-2, 1}, { 2,-2, 1}, {-2,-1, 1}, {-1,-1, 1}, { 0,-1, 1},
                                    { 1,-1, 1}, { 2,-1, 1}, {-2, 0, 1}, {-1, 0, 1}, { 0, 0, 1},
                                    { 1, 0, 1}, { 2, 0, 1}, {-2, 1, 1}, {-1, 1, 1}, { 0, 1, 1},
                                    { 1, 1, 1}, { 2, 1, 1}, {-2, 2, 1}, {-1, 2, 1}, { 0, 2, 1},
                                    { 1, 2, 1}, { 2, 2, 1}, {-2,-2, 2}, {-1,-2, 2}, { 0,-2, 2},
                                    { 1,-2, 2}, { 2,-2, 2}, {-2,-1, 2}, {-1,-1, 2}, { 0,-1, 2},
                                    { 1,-1, 2}, { 2,-1, 2}, {-2, 0, 2}, {-1, 0, 2}, { 0, 0, 2},
                                    { 1, 0, 2}, { 2, 0, 2}, {-2, 1, 2}, {-1, 1, 2}, { 0, 1, 2},
                                    { 1, 1, 2}, { 2, 1, 2}, {-2, 2, 2}, {-1, 2, 2}, { 0, 2, 2},
                                    { 1, 2, 2}, { 2, 2, 2} };

typedef struct {
  int neighbor[nbcells];
} tCell;

typedef struct {
  int nitems;
  int count;
  int *first;
  int *last;
  int *item;
} tList;

typedef struct {
  int builds;        // Number of neighbor-list builds

  tList neighbor;
  tList exclude;

  double time;

  int natoms;        // Number of atoms
  int nx3;           // Three times the number of atoms
  int mcells;        // Number of cells at each dimension
  int ncells;        // Total number of cells
  int maxcells;      // Maximum number of cells

  double Rc;         // Cut-off distance
  double RcSq;       // Cut-off distance squared
  double xRc;        // Extended cutoff distance (including skin)
  double xRcSq;      // Extended cutoff distance squared
  double skinSq;     // Square of the neighbor list skin width

  tCell *cell;

  int *type;         // Atom types
  double *R0;        // Atom positions at list building
  double *R;         // Pointer to dynamic atom positions
  double *F;

  int ntypes;
  tPairType *pairType;

  double Energy;
  double Virial;

} tEmDee;

void md_initialize( tEmDee *me, double rc, double skin, int atoms, int types,
                    int *type_index, double *coords, double *forces );
void md_set_pair( tEmDee *me, int i, int j, tPairType model );
void md_compute_forces( tEmDee *me, double L );


tPairType lennard_jones( double sigma, double epsilon );
tPairType shifted_force_lennard_jones( double sigma, double epsilon, double rc );

