typedef struct {
  void *data;
  void *params;
} EmDee_Model;

typedef struct {

  int builds;            // Number of neighbor-list builds
  double time;           // Total time taken in force calculations
  double Rc;             // Cut-off distance
  double RcSq;           // Cut-off distance squared
  double xRc;            // Extended cutoff distance (including skin)
  double xRcSq;          // Extended cutoff distance squared
  double skinSq;         // Square of the neighbor list skin width

  int mcells;            // Number of cells at each dimension
  int ncells;            // Total number of cells
  int maxcells;          // Maximum number of cells
  int maxatoms;
  int maxpairs;          // Maximum number of pairs containing all atoms of a cell
  void *cell;            // Array containing all neighbor cells of each cell

  int natoms;            // Number of atoms in the system
  int *type;             // The type of each atom
  double *R0;            // The position of each atom at the latest neighbor list building
  double *charge;        // Pointer to the electric charge of each atom
  int chargeFlag;

  int ntypes;            // Number of atom types
  void *pairParams;      // Model parameters of each type of atom pair
  void *pairData;        // Model data of each type of atom pair

  void *bond;            // List of bonds
  void *angle;           // List of angles
  void *dihedral;        // List of dihedrals

  double Energy;         // Total potential energy of the system
  double Virial;         // Total internal virial of the system

  int nthreads;          // Number of parallel openmp threads
  void *cellAtom;        // List of atoms belonging to each cell
  void *threadCell;      // List of cells to be dealt with in each parallel thread
  void *neighbor;        // Pointer to neighbor lists
  void *excluded;        // List of pairs excluded from the neighbor lists

} tEmDee;

tEmDee EmDee_system( int threads, double rc, double skin, int N, int *types );

void EmDee_set_charges( tEmDee *md, double *charges );

void EmDee_set_pair( tEmDee *md, int itype, int jtype, EmDee_Model *model );

void EmDee_apply_mixing_rules( tEmDee *md );

void EmDee_add_bond( tEmDee *md, int i, int j, EmDee_Model *model );

void EmDee_add_angle( tEmDee *md, int i, int j, int k, EmDee_Model *model );

void EmDee_add_dihedral( tEmDee *md, int i, int j, int k, int l, EmDee_Model *model );

void EmDee_exclude_pair( tEmDee *md, int i, int j );

void EmDee_compute_forces( tEmDee *md, double *forces, double *coords, double L );

EmDee_Model EmDee_pair_lj( double sigma, double epsilon );

EmDee_Model EmDee_pair_lj_sf( double sigma, double epsilon, double cutoff );

EmDee_Model EmDee_pair_lj_coul_sf( double sigma, double epsilon, double cutoff );

EmDee_Model EmDee_bond_harmonic( double k, double r0 );

EmDee_Model EmDee_bond_morse( double D, double alpha, double r0 );

EmDee_Model EmDee_angle_harmonic( double k, double theta0 );

EmDee_Model EmDee_dihedral_harmonic( double k, double phi0 );

