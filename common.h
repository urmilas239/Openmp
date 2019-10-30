#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;



//
// particle data structure
//
typedef struct 
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
  double index;
} particle_t;


typedef struct 
{
	int bin_id;
	int num_particles;
	int particle_offset;
	int neighbor_id[8] ={-1, -1,-1,-1,-1,-1,-1,-1};
}particle_bin_mapping;


typedef struct 
{
	int neighbor_bin_id;
	int num_particles;
	int particle_offset;
}neighbor_bin_mapping;


//Added 

 void set_bin_count(int n);
 int compute_bin_index_from_xy(double x, double y);
 std::vector<std::vector<int> > initialize_neighbor_bins();
std::vector<std::vector<int> > initialize_bin_vector();
//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles1( int n, particle_t *p, std::vector<std::vector<int> > &bin_map);
void init_particles( int n, particle_t *p);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );
void move1( particle_t &p,  std::vector<std::vector<int> > &bin_map);
void bin_particles(int n, particle_t *p ,  std::vector<std::vector<int> > &bin_map);

//MPI Functions
void form_particles_array_for_MPI(std::vector<int>  &bin_ids,
    std::vector<int>  &boundary_bin_ids,
    std::vector<std::vector<int> > &bin_map,
    std::vector<std::vector<int> > &neighbor_bins,
    particle_t *particles_to_send, 
    particle_bin_mapping *pbm,
    neighbor_bin_mapping *n_bins,
    particle_t *particles,
    int &partition_size_per_rank);

std::vector<std::vector<int> > get_boundary_bins(std::vector<std::vector<int> > process_bins, std::vector<std::vector<int> > neighbor_bins);
std::vector<std::vector<int> > assign_bins_to_processes_mpi(int num_of_processes, std::vector<std::vector<int> > bin_map);




//MPI alternate logic functions
std::vector<int > assign_bins_to_current_process_mpi(int num_of_processes, int current_process, std::vector<std::vector<int> > bin_map, std::vector<int> bin_process_map);
std::vector<int> get_boundary_bins_for_curr_process(std::vector<int > process_bins, std::vector<std::vector<int> > neighbor_bins);

//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
