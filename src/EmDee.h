typedef struct {
  int neighbor[58];
} tCell;

typedef struct {
  int builds;        // Number of neighbor-list builds
  int *first;        // First neighbor of each atom
  int *last;         // Last neighbor of each atom 
  int *neighbor;     // List of neighbors

  double time, checktime; // APAGAR DEPOIS DE TESTAR

  int natoms;        // Number of atoms
  int nx3;           // Three times the number of atoms
  int npairs;        // Number of neighbor pairs
  int maxpairs;      // Maximum number of neighbor pairs
  int mcells;        // Number of cells at each dimension
  int ncells;        // Total number of cells
  int maxcells;      // Maximum number of cells

  double RcSq;       // Cut-off distance squared
  double xRc;        // Extended cutoff distance (including skin)
  double xRcSq;      // Extended cutoff distance squared
  double skinSq;     // Square of the neighbor list skin width

  tCell *cell;
  int *next;         // Next atom in a linked list of atoms

  int *type;         // Atom types
  double *R0;        // Atom positions at list building
  double *R;         // Pointer to dynamic atom positions
  double *P;         // Pointer to dynamic atom momenta
  double *F;

  double Energy;
  double Virial;
} tEmDee;

const int minMcell = 5;
const int extraPairs = 1000;

const int nb[58][3] = { { 0, 0,-1}, { 0,-1, 0}, {-1, 0, 0}, { 0,-1,-1}, {-1, 0,-1}, { 0, 1,-1},
                        { 1, 0,-1}, { 1,-1, 0}, {-1,-1, 0}, {-1,-1,-1}, {-1, 1,-1}, { 1, 1,-1},
                        { 1,-1,-1}, { 0, 0,-2}, { 0,-2, 0}, {-2, 0, 0}, { 1, 0,-2}, { 2,-1, 0},
                        { 0, 2,-1}, { 0,-1,-2}, {-1, 0,-2}, { 0, 1,-2}, { 0,-2,-1}, {-2, 0,-1},
                        { 2, 0,-1}, {-1,-2, 0}, { 1,-2, 0}, {-2,-1, 0}, {-1,-1,-2}, {-1, 1,-2},
                        { 1, 1,-2}, {-1,-2,-1}, { 1,-2,-1}, {-2,-1,-1}, { 2,-1,-1}, { 1,-1,-2},
                        {-2, 1,-1}, { 2, 1,-1}, { 1, 2,-1}, {-1, 2,-1}, { 0,-2,-2}, {-2, 0,-2},
                        {-2,-2, 0}, { 2,-2, 0}, { 0, 2,-2}, {-1,-2,-2}, { 1,-2,-2}, {-2,-1,-2},
                        { 2,-1,-2}, { 2, 0,-2}, {-2, 1,-2}, { 2, 1,-2}, {-1, 2,-2}, { 1, 2,-2},
                        {-2,-2,-1}, { 2,-2,-1}, {-2, 2,-1}, { 2, 2,-1} };

void md_initialize( tEmDee *me, double rc, double skin, int atoms, int *types );
void md_upload( tEmDee *me, double *coords, double *momenta );
void md_change_coordinates( tEmDee *me, double a, double b );
void md_change_momenta( tEmDee *me, double a, double b );



void md_compute_lennard_jones( tEmDee *me, double L );
void md_download( tEmDee *me, double *coords, double *momenta, double *forces );


void md_handle_neighbor_list( tEmDee *me, double L );
void md_scale_momenta( tEmDee *me, double scaling );
void md_translate_momenta( tEmDee *me, double force_factor );
void md_scale_coords( tEmDee *me, double scaling );
void md_translate_coords( tEmDee *me, double momentum_factor );

/* -------------------------------------------------------------------------------------------------
                               F O R T R A N    B I N D I N G S
------------------------------------------------------------------------------------------------- */
void md_allocate( void **me );
void md_capture_first_and_last( tEmDee *me, void **first, void **last );
int md_capture_neighbor( tEmDee *me, void **neighbor );

