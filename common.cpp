#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include "common.h"

double size;
int num_of_bins_x;
int num_of_bins_y;
int total_bin_count;
int bin_density;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005



//
//  timer
//
double read_timer()
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

std::vector<std::vector<int>> initialize_bin_vector()
{
    std::vector<std::vector<int>> bin_map(total_bin_count, std::vector<int>(bin_density+2)); 
    return bin_map;
}


std::vector<std::vector<int>> initialize_neighbor_bins()
{
     std::vector<std::vector<int>> neighbor_bins(total_bin_count, std::vector<int>(8)); 
     int index;
     int last_row = num_of_bins_y*(num_of_bins_y-1);

     //1st row 
     
     
     for(int i=0; i<num_of_bins_y;i++)
     { 
        index =0;




       if(i%num_of_bins_y != 0)
       {
        //not first column
          neighbor_bins[i][index] =  i-1;
          index++;

          if(i > num_of_bins_y)
          {
            //not first row
            neighbor_bins[i][index] =  i-1-num_of_bins_y;
            index++;
          }

          if(i>= last_row)
          {
            neighbor_bins[i][index] =  i-1+num_of_bins_y;
            index++;
          }
          
         
       }


       if((i+1)%num_of_bins_y != 0)
       {
        //not last column
          neighbor_bins[i][index] =  i+1;
          index++;

          if(i > num_of_bins_y)
          {
            neighbor_bins[i][index] =  i+1-num_of_bins_y;
            index++;
          }

          if(i>= last_row)
          {
            neighbor_bins[i][index] =  i+1+num_of_bins_y;
            index++;
          } 
          
       }


       if(i > num_of_bins_y)
          {
            neighbor_bins[i][index] =  i-num_of_bins_y;
            index++;
          }

          if(i>= last_row)
          {
            neighbor_bins[i][index] =  i+num_of_bins_y;
            index++;
          } 
      return neighbor_bins;
        
     }

     

    

     return neighbor_bins;
}

/*
 * Set number of bins
 */
 void set_bin_count(int n)
 {
    num_of_bins_x = ceil(size/min_r);
    num_of_bins_y = num_of_bins_x;
    total_bin_count = num_of_bins_x*num_of_bins_y;
    bin_density = n/total_bin_count;
 }


 int compute_bin_index_from_xy(int x, int y, int bin_size)
 {
    int bin_index = -1;
    bin_index = (floor(x/bin_size)*num_of_bins_y)+floor(y /bin_size);
    return bin_index;

 }

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p ,  std::vector<std::vector<int>> &bin_map)
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];
        
        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //TODO: Assign bin for the particle
        //bin_map.push_back() vvvv
        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
	if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
