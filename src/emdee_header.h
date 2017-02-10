typedef struct {
  int    builds;         // Number of neighbor-list builds
  double pairTime;       // Time taken in force calculations
  double totalTime;      // Total time since initialization
  struct {
    double Potential;    // Total potential energy of the system
    double Dispersion;   // Dispersion (vdW) part of the potential energy
    double Coulomb;      // Electrostatic part of the potential energy
    double Fourier;      // Reciprocal part of the electrostatic potential
    double Kinetic;      // Total kinetic energy of the system
    double KinPart[3];   // Kinetic energy at each dimension
    double Rotational;   // Rotational kinetic energy of the system
    double RotPart[3];   // Rotational kinetic energy around each principal axis
    _Bool  Compute;      // Flag to activate/deactivate energy computations
    _Bool  UpToDate;     // Flag to attest whether energies have been computed
  } Energy;
  double Virial;         // Total internal virial of the system
  double BodyVirial;     // Rigid body contribution to the internal virial
  int    DOF;            // Total number of degrees of freedom
  int    RDOF;           // Number of rotational degrees of freedom
  void*  Data;           // Pointer to system data
  struct {
    _Bool translate;     // Flag to activate/deactivate translations
    _Bool rotate;        // Flag to activate/deactivate rotations
    int   rotationMode;  // Algorithm used for free rotation of rigid bodies
  } Options;
} tEmDee;

tEmDee EmDee_system( int threads, int layers, double rc, double skin, int N, int* types, 
                     double* masses, int* bodies );

void EmDee_switch_model_layer( tEmDee* md, int layer );

void EmDee_set_pair_model( tEmDee md, int itype, int jtype, void* model, double kCoul );

void EmDee_set_pair_multimodel( tEmDee md, int itype, int jtype, void* model[], double kCoul[] );

void EmDee_set_kspace_model( tEmDee md, void* model );

void EmDee_set_coul_model( tEmDee md, void* model );

void EmDee_set_coul_multimodel( tEmDee md, void* model[] );

void EmDee_ignore_pair( tEmDee md, int i, int j );

void EmDee_add_bond( tEmDee md, int i, int j, void* model );

void EmDee_add_angle( tEmDee md, int i, int j, int k, void* model );

void EmDee_add_dihedral( tEmDee md, int i, int j, int k, int l, void* model );

void EmDee_upload( tEmDee* md, char *option, double* address );

void EmDee_download( tEmDee md, char *option, double* address );

void EmDee_random_momenta( tEmDee* md, double kT, _Bool adjust, int seed );

//void EmDee_save_state( tEmDee md, int rigid );
//void EmDee_restore_state( tEmDee md );

void EmDee_boost( tEmDee* md, double lambda, double alpha, double dt );

void EmDee_move( tEmDee* md, double lambda, double alpha, double dt );

void* EmDee_pair_none();

void* EmDee_coul_none();

void* EmDee_bond_none();

void* EmDee_angle_none();

void* EmDee_dihedral_none();

