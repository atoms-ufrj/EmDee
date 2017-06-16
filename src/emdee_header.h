typedef struct {
  int    Builds;             // Number of neighbor-list builds
  struct {
    double Pair;             // Time taken in force calculations
    double FastPair;
    double Neighbor;
    double Total;            // Total time since initialization
  } Time;
  struct {
    double Potential;        // Total potential energy of the system
    double Dispersion;       // Dispersion (vdW) part of the potential energy
    double Coulomb;          // Electrostatic part of the potential energy
    double Fourier;          // Reciprocal part of the electrostatic potential
    double Kinetic;          // Total kinetic energy of the system
    double TransPart[3];     // Translational kinetic energy at each dimension
    double Rotational;       // Rotational kinetic energy of the system
    double RotPart[3];       // Rotational kinetic energy around each principal axis
    double *LayerPotential;  // Vector with multilayer potential energy components
    double *LayerDispersion; // Vector with multilayer dispersion energy components
    double *LayerCoulomb;    // Vector with multilayer coulombic energy components
    double *LayerFourier;    // Vector with multilayer reciprocal energy components
    _Bool  Compute;          // Flag to activate/deactivate energy computations
    _Bool  UpToDate;         // Flag to attest whether energies have been computed
  } Energy;
  double Virial;             // Total internal virial of the system
  double BodyVirial;         // Rigid body contribution to the internal virial
  int    DoF;                // Total number of degrees of freedom
  int    RotDoF;             // Number of rotational degrees of freedom
  void*  Data;               // Pointer to system data
  struct {
    _Bool  Translate;        // Flag to activate/deactivate translations
    _Bool  Rotate;           // Flag to activate/deactivate rotations
    int    RotationMode;     // Algorithm used for free rotation of rigid bodies
    _Bool  AutoForceCompute; // Flag to activate/deactivate automatic force computations
    _Bool  AutoBodyUpdate;   // Flag to activate/deactivate automatic rigid body update
  } Options;
} tEmDee;

tEmDee EmDee_system( int threads, int layers, double rc, double skin, int N, int* types, 
                     double* masses, int* bodies );

void EmDee_share_phase_space( tEmDee mdkeep, tEmDee* mdlose );

void EmDee_switch_model_layer( tEmDee md, int layer );

void EmDee_set_pair_model( tEmDee md, int itype, int jtype, void* model, double kCoul );

void EmDee_set_pair_multimodel( tEmDee md, int itype, int jtype, void* model[], double kCoul[] );

void EmDee_set_kspace_model( tEmDee md, void* model );

void EmDee_set_coul_model( tEmDee md, void* model );

void EmDee_set_coul_multimodel( tEmDee md, void* model[] );

void EmDee_ignore_pair( tEmDee md, int i, int j );

void EmDee_add_bond( tEmDee md, int i, int j, void* model );

void EmDee_add_angle( tEmDee md, int i, int j, int k, void* model );

void EmDee_add_dihedral( tEmDee md, int i, int j, int k, int l, void* model );

void EmDee_set_respa( tEmDee md, double InRc, double ExRc, int Npair, int Nbond );

void EmDee_upload( tEmDee* md, char *option, double* address );

void EmDee_download( tEmDee md, char *option, double* address );

void EmDee_random_momenta( tEmDee* md, double kT, _Bool adjust, int seed );

void EmDee_boost( tEmDee* md, double lambda, double alpha, double dt );

void EmDee_displace( tEmDee* md, double lambda, double alpha, double dt );

void EmDee_compute_forces( tEmDee* md );

void EmDee_advance( tEmDee* md, double alpha_R, double alpha_P, double dt );

void EmDee_rdf( tEmDee md, int bins, int pairs, int itype[], int jtype[], double count[] );

void* EmDee_shifted_force( void* model );

void* EmDee_smoothed( void* model, double Rm );

void* EmDee_pair_none();

void* EmDee_coul_none();

void* EmDee_bond_none();

void* EmDee_angle_none();

void* EmDee_dihedral_none();

