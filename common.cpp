#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <vector>
#include <iostream>
#include "common.h"

double size;
int num_of_bins_x;
int num_of_bins_y;
int total_bin_count;
int bin_density;
double bin_size;

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

std::vector<std::vector<int> > initialize_bin_vector()
{
    std::vector<std::vector<int> > bin_map(total_bin_count, std::vector<int>()); 
    return bin_map;
}


std::vector<std::vector<int> > initialize_neighbor_bins()
{
    std::cout << ":::IN neighbor_bins::: " << std::endl;
     std::vector<std::vector<int> > neighbor_bins(total_bin_count, std::vector<int>(8, -1)); 

     std::cout << ":::neighbor_bins.size()::: " << neighbor_bins.size()<<std::endl;
     int index;
     int last_row = num_of_bins_y*(num_of_bins_y-1);
     std::cout << ":::last_row::: " << last_row<<std::endl;

    

     
     for(int i=0; i<total_bin_count;i++)
     { 
        index =0;

      /*std::cout << ":::i  :: " <<i<<std::endl;
      std::cout << ":::i mod num_of_bins_y :: " << (i%num_of_bins_y ) <<std::endl;
      std::cout << ":::(i+1) mod num_of_bins_y :: " << ((i+1)%num_of_bins_y) <<std::endl;*/

       if(i%num_of_bins_y != 0)
       {
        //not first column
          neighbor_bins[i][index] =  i-1;
          index++;

          if(i >= num_of_bins_y)
          {
            //not first row
            neighbor_bins[i][index] =  i-1-num_of_bins_y;
            index++;
          }

          if(i < last_row)
          {
            // Not last row
            neighbor_bins[i][index] =  i-1+num_of_bins_y;
            index++;
          }
          
         
       }


       if((i+1)%num_of_bins_y != 0)
       {
        //not last column
          neighbor_bins[i][index] =  i+1;
          index++;

          if(i >= num_of_bins_y)
          {
            //not first row
            neighbor_bins[i][index] =  i+1-num_of_bins_y;
            index++;
          }

          if(i < last_row)
          {
            neighbor_bins[i][index] =  i+1+num_of_bins_y;
            index++;
          } 
          
       }


       if(i >= num_of_bins_y)
          {
            neighbor_bins[i][index] =  i-num_of_bins_y;
            index++;
          }

          if(i < last_row)
          {
            neighbor_bins[i][index] =  i+num_of_bins_y;
            index++;
          } 
        
     }

    

      /*
      std::vector<int>  my_neighbors;
        for(int i = 0; i < neighbor_bins.size(); i++)
        {
            std::cout << "Current bin::: " << i << ":: My neighbors:::" << std::endl;
            my_neighbors = neighbor_bins.at(i);
            if(!my_neighbors.empty())
            {
                for(int j = 0; j < my_neighbors.size(); j++)
                {
                    std::cout << "::: " << my_neighbors.at(j) << std::endl;
                }
            }
            

        }
    */

     return neighbor_bins;
}

/*
 * Set number of bins
 */
 void set_bin_count(int n)
 {

    std::cout<< "n::" <<n<<std::endl;
    std::cout<< "size::" <<size<<std::endl;

    bin_size = cutoff;
    std::cout<< "bin_size::" <<bin_size<<std::endl;
    num_of_bins_x = ceil(size/bin_size);
    std::cout<< "num_of_bins_x::" <<num_of_bins_x<<std::endl;
    num_of_bins_y = num_of_bins_x;
    std::cout<< "num_of_bins_y::" <<num_of_bins_y<<std::endl;
    total_bin_count = num_of_bins_x*num_of_bins_y;
    std::cout<< "total_bin_count::" <<total_bin_count<<std::endl;
    bin_density = n/total_bin_count;
    std::cout<< "bin_density::" <<bin_density<<std::endl;
    
 }


 int compute_bin_index_from_xy(double x, double y)
 {
    
    int bin_index = (floor(x/bin_size)*num_of_bins_y)+floor(y /bin_size);
    return bin_index;

 }


 void bin_particles(int n, particle_t *p ,  std::vector<std::vector<int> > &bin_map)
 {

    int bin_index;
    std::vector<int> particle_list;
    for( int i = 0; i < n; i++ ) 
    {
        //particle_list.clear();
         bin_index = compute_bin_index_from_xy( p[i].x, p[i].y );
         std::vector<int> particle_list = bin_map.at(bin_index);
            particle_list.push_back(i);
            bin_map.insert(bin_map.begin() + bin_index, particle_list);

         
         
    }
 }

//
//  Initialize the particle positions and velocities
//
void init_particles1( int n, particle_t *p ,  std::vector<std::vector<int> > &bin_map)
{
    srand48( time( NULL ) );
        
    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    std::vector<int> particle_list;
    int bin_index = -1;
    for( int i = 0; i < n; i++ ) 
    {
        particle_list.clear();
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
        std::cout<< "p[i].x:: " << p[i].x << ",:: p[i].y::" << p[i].y << std::endl;
        //TODO: Assign bin for the particle
         bin_index = compute_bin_index_from_xy( p[i].x, p[i].y );
         particle_list = bin_map.at(bin_index);
         particle_list.push_back(i);
         bin_map.insert(bin_map.begin() + bin_index, particle_list);

         std::cout<< "init_particles1 ::: bin_index::: " << bin_index << std::endl;



        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
        p[i].index = i;
    }
    free( shuffle );
}





void init_particles( int n, particle_t *p)
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

        //p[i].index = i;
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




void move1( particle_t &p, std::vector<std::vector<int> > &bin_map)
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //

    //Remove the current particle from its old location
    int old_bin_index = compute_bin_index_from_xy( p.x, p.y);
    //std::cout<< "old_bin_index::" <<old_bin_index<<std::endl;
    

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

     
     //Insert the particle into its new bin based on its new position cordinates
     int bin_index = compute_bin_index_from_xy( p.x, p.y);
     //std::cout<< "bin_index::" <<bin_index<<std::endl;
     if(old_bin_index != bin_index)
     {
        std::vector<int> particle_list = bin_map.at(old_bin_index);
        std::remove(particle_list.begin(), particle_list.end(), p.index); 
        bin_map.insert(bin_map.begin() + old_bin_index, particle_list);


        particle_list.clear();
      particle_list = bin_map.at(bin_index);
      particle_list.push_back(p.index);
      bin_map.insert(bin_map.begin() + bin_index, particle_list);
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
