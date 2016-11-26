typedef struct {
  int    builds;         // Number of neighbor-list builds
  double pairTime;       // Time taken in force calculations
  double totalTime;      // Total time since initialization
  double Potential;      // Total potential energy of the system
  double Kinetic;        // Total kinetic energy of the system
  double Rotational;     // Rotational kinetic energy of the system
  double Virial;         // Total internal virial of the system
  int    DOF;            // Total number of degrees of freedom
  int    RDOF;           // Number of rotational degrees of freedom
  void*  Data;           // Pointer to system data
  struct {
    int translate;       // Flag to activate/deactivate translations
    int rotate;          // Flag to activate/deactivate rotations
    int rotationMode;    // Algorithm used for free rotation of rigid bodies
  } Options;
} tEmDee;

tEmDee EmDee_system( int threads, int layers, double rc, double skin, int N, int* types, double* masses );
void EmDee_switch_model_layer( tEmDee* md, int layer );
void EmDee_set_charges( tEmDee md, double* charges );
void EmDee_set_pair_type( tEmDee md, int itype, int jtype, void* model );
void EmDee_ignore_pair( tEmDee md, int i, int j );
void EmDee_add_bond( tEmDee md, int i, int j, void* model );
void EmDee_add_angle( tEmDee md, int i, int j, int k, void* model );
void EmDee_add_dihedral( tEmDee md, int i, int j, int k, int l, void* model );
void EmDee_add_rigid_body( tEmDee md, int N, int* indexes );
void EmDee_upload( tEmDee* md, double* Lbox, double* coords, double* momenta, double* forces );
void EmDee_download( tEmDee md, double* Lbox, double* coords, double* momenta, double* forces );
void EmDee_random_momenta( tEmDee* md, double kT, int adjust, int seed );
//void EmDee_save_state( tEmDee md, int rigid );
//void EmDee_restore_state( tEmDee md );
void EmDee_boost( tEmDee* md, double lambda, double alpha, double dt );
void EmDee_move( tEmDee* md, double lambda, double alpha, double dt );
void EmDee_group_energy( tEmDee md, int na, double* atoms, int ne, double* energies );
void* EmDee_pair_none();
void* EmDee_bond_none();
void* EmDee_angle_none();
void* EmDee_dihedral_none();
