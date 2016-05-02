typedef struct {
  int model;
  double p1;
  double p2;
  double p3;
  double p4;
} tModel;

typedef struct {
  int i;
  int j;
  tModel model;
} tBond;

typedef struct {
  int i;
  int j;
  tModel model;
} tAngle;

typedef struct {
  int neighbor[62];
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
  tModel *pairType;

  int nbonds;
  tBond *bond;

  double Energy;
  double Virial;

} tEmDee;

void md_initialize( tEmDee *me, double rc, double skin, int atoms, int types,
                    int *type_index, double *coords, double *forces );
void md_set_pair( tEmDee *me, int i, int j, tModel model );
void md_compute_forces( tEmDee *me, double L );


tModel lennard_jones( double sigma, double epsilon );
tModel shifted_force_lennard_jones( double sigma, double epsilon, double rc );

