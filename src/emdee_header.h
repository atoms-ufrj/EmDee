typedef struct {
  int    Builds;             // Number of neighbor-list builds
  struct {
    double Pair;             // Time taken in force calculations
    double Motion;
    double Neighbor;
    double Total;            // Total time since initialization
  } Time;
  struct {
    double Potential;        // Total potential energy of the system
    double Dispersion;       // Dispersion (vdW) part of the potential energy
    double Coulomb;          // Electrostatic part of the potential energy
    double Bond;
    double Angle;
    double Dihedral;
    double ShadowPotential;
    _Bool  UpToDate;         // Flag to attest whether energies have been computed
  } Energy;
  struct {
    double Kinetic;          // Total kinetic energy of the system
    double TransPart[3];     // Translational kinetic energy at each dimension
    double Rotational;       // Rotational kinetic energy of the system
    double RotPart[3];       // Rotational kinetic energy around each principal axis
    double ShadowKinetic;
    double ShadowRotational;
  } Kinetic;
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
    _Bool  Compute;          // Flag to activate/deactivate energy computations
  } Options;
} tEmDee;

tEmDee EmDee_system( int threads, int layers, double rc, double skin, int N, int* types,
                     double* masses, int* bodies );

void* EmDee_memory_address( tEmDee md, char *option );

void EmDee_share_phase_space( tEmDee mdkeep, tEmDee* mdlose );

void EmDee_layer_based_parameters( tEmDee md, double InternalRc, int* Apply, int* Bonded );

void EmDee_set_pair_model( tEmDee md, int itype, int jtype, void* model, double kCoul );

void EmDee_set_pair_multimodel( tEmDee md, int itype, int jtype, void* model[], double kCoul[] );

void EmDee_set_kspace_model( tEmDee md, void* model );

void EmDee_set_coul_model( tEmDee md, void* model );

void EmDee_set_coul_multimodel( tEmDee md, void* model[] );

void EmDee_ignore_pair( tEmDee md, int i, int j );

void EmDee_add_bond( tEmDee md, int i, int j, void* model );

void EmDee_add_angle( tEmDee md, int i, int j, int k, void* model );

void EmDee_add_dihedral( tEmDee md, int i, int j, int k, int l, void* model );

void EmDee_download( tEmDee md, char *option, double* address );

void EmDee_upload( tEmDee* md, char *option, double* address );

void EmDee_switch_model_layer( tEmDee* md, int layer );

void EmDee_random_momenta( tEmDee* md, double kT, _Bool adjust, int seed );

void EmDee_boost( tEmDee* md, double lambda, double alpha, double dt );

void EmDee_displace( tEmDee* md, double lambda, double alpha, double dt );

void EmDee_verlet_step( tEmDee* md, double dt );

void EmDee_compute_forces( tEmDee* md );

void EmDee_rdf( tEmDee md, int bins, double Rc, int pairs, int itype[], int jtype[], double count[] );

void* EmDee_shifted_force( void* model );

void* EmDee_smoothed( void* model, double skin );

void* EmDee_shifted_smoothed( void* model, double skin );

void* EmDee_square_smoothed( void* model, double skin );

void* EmDee_shifted_square_smoothed( void* model, double skin );

void* EmDee_pair_none();

void* EmDee_coul_none();

void* EmDee_bond_none();

void* EmDee_angle_none();

void* EmDee_dihedral_none();
